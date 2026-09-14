"""Can planar matchgate contraction reproduce signed two-adder coefficients?

PREDICTIONS BEFORE EXECUTION: complement parity makes local carry factors
ternary matchgates; explicit gadgets realize every tested deletion signature;
FKT with signed calibration equals independent tiny perfect-matching sums;
the two-adder word map matches actual Cuccaro/SWAP gate replay on all inputs;
selected signed coefficients agree exactly for every width3 bit permutation.
Controls: naive Pfaffian orientation on K4, magnitude-only cancellation and
accepting nonplanar K3,3 must fail. No width/rank advantage is predicted.

Fix width3 and the same20 mask fixtures; vary only the intermediate bit
permutation (all6). Full circuit has13qubits,8192basis labels. No unitary,
tensor SVD or large memory pilot. Tiny matching enumeration is limited to
local gadgets or4nodes. Scalar simulation uses the known planar FKT method.

Run: uv run --with networkx==3.5 python -m experiments.experiment_planar_carry
"""
from fractions import Fraction
from itertools import permutations

import networkx as nx

from circuits import Circuit, ripple_adder
from lab.harness import Experiment
from lab.planar_carry import (_Gadgets, NonPlanarCarry, carry_factors, hadamard,
                              pfaffian, planar_perfect_matching,
                              two_addition_coefficient)
from walsh import classical_permutation


def direct_matching(graph):
    if not graph:
        return Fraction(1)
    first = next(iter(graph))
    result = Fraction(0)
    for neighbor in list(graph[first]):
        smaller = graph.copy()
        smaller.remove_nodes_from((first, neighbor))
        result += graph[first][neighbor].get("weight", 1)*direct_matching(smaller)
    return result


def actual_circuit(perm):
    w = len(perm)
    qc = Circuit(3*w+4)
    adder, regs = ripple_adder(w)
    def add(a_start, carry, z):
        mapping = {regs["c0"]: carry, regs["z"]: z}
        mapping.update({q: a_start+i for i, q in enumerate(regs["a"])})
        mapping.update({q: w+i for i, q in enumerate(regs["b"])})
        for op in adder.logical:
            getattr(qc, op[0])(*(mapping[q] for q in op[1:]))
    add(0, 3*w, 3*w+2)
    contents = list(range(w))
    for i, source in enumerate(perm):
        j = contents.index(source)
        if i != j:
            qc.swap(w+i, w+j)
            contents[i], contents[j] = contents[j], contents[i]
    add(2*w, 3*w+1, 3*w+3)
    return qc


def expected_image(x, perm):
    w, mask = len(perm), (1 << len(perm))-1
    a, b, c = x & mask, (x >> w) & mask, (x >> (2*w)) & mask
    first = a+b+((x >> (3*w)) & 1)
    shuffled = sum(((first >> source) & 1) << i for i, source in enumerate(perm))
    second = shuffled+c+((x >> (3*w+1)) & 1)
    result = (x & ~(mask << w)) | ((second & mask) << w)
    return result ^ ((first >> w) << (3*w+2)) ^ ((second >> w) << (3*w+3))


def run():
    exp = Experiment("planar_carry", doc=__doc__)
    exp.predict("P1", "carry factors transform to pure parity and all local gadget signatures agree")
    exp.predict("P2", "signed planar FKT equals independent tiny perfect-matching sums")
    exp.predict("P3", "actual Cuccaro/SWAP circuits equal the full-space word map")
    exp.predict("P4", "all selected exact signed word coefficients equal full-space gate replay")
    exp.must_fail("C1", "unoriented K4 Pfaffian must disagree with its three perfect matchings")
    exp.must_fail("C2", "replacing negative edge weights by magnitudes must destroy exact cancellation")
    exp.must_fail("C3", "ordinary nonplanar K3,3 must be rejected")
    fixtures = []
    for a in (0, 1):
        for b in (0, 1):
            k, r = carry_factors(a, b, a, b)
            fixtures.extend((hadamard(k), hadamard(r)))
    fixtures.extend(([7,0,0,0,0,0,0,0], [0,0,0,0,0,0,0,-3],
                     [0,0,0,1,0,2,-1,0], [0,1,-2,0,3,0,0,0], [0]*8))
    gadget_errors = 0
    for values in fixtures:
        parity = {i.bit_count() % 2 for i, v in enumerate(values) if v}
        gadget_errors += int(len(parity) > 1)
        g = _Gadgets()
        ports = g.signature(values)
        for bits, expected in enumerate(values):
            graph = g.graph.copy()
            graph.remove_nodes_from(ports[j] for j in range(3) if bits & (1 << j))
            actual = g.scale*direct_matching(graph)
            gadget_errors += int(actual != expected)
    exp.check("P1", gadget_errors == 0, f"{len(fixtures)*8} deletion signatures; errors={gadget_errors}")

    matching_errors = 0
    for edge_weight in (-2, -1, 1, 2):
        graph = nx.complete_graph(4)
        nx.set_edge_attributes(graph, 1, "weight")
        graph[0][1]["weight"] = edge_weight
        expected = direct_matching(graph)
        actual = planar_perfect_matching(graph)
        matching_errors += int(actual != expected)
    exp.check("P2", matching_errors == 0, f"four signed K4 fixtures; errors={matching_errors}")
    graph = nx.complete_graph(4)
    naive = [[0 if i == j else 1 if i < j else -1 for j in range(4)] for i in range(4)]
    exp.fail_check("C1", pfaffian(naive) != direct_matching(graph),
                   f"naive Pfaffian={pfaffian(naive)}, matching sum={direct_matching(graph)}")
    graph[0][1]["weight"] = -2
    signed = planar_perfect_matching(graph)
    graph[0][1]["weight"] = 2
    absolute = planar_perfect_matching(graph)
    exp.fail_check("C2", signed == 0 and absolute == 4, f"signed={signed}, magnitudes={absolute}")
    rejected = False
    try:
        planar_perfect_matching(nx.complete_bipartite_graph(3, 3))
    except NonPlanarCarry:
        rejected = True
    exp.fail_check("C3", rejected, "K3,3 rejects before Pfaffian extraction")

    w = 3
    masks = [(4+(i & 3), 4+((i >> 2) & 3), 4+((3*i+1) & 3),
              4+((i >> 1) & 3), i & 1, (i >> 1) & 1) for i in range(16)]
    masks += [(0,0,0,0,0,0), (1,1,1,1,1,1), (2,2,2,2,0,0), (3,2,2,2,1,0)]
    rows = []
    map_errors = coefficient_errors = negatives = nonzero = 0
    for perm in permutations(range(w)):
        qc = actual_circuit(perm)
        images = classical_permutation(qc)
        bad = sum(int(y) != expected_image(x, perm) for x, y in enumerate(images))
        map_errors += bad
        for alpha, beta, gamma, delta, lam, mu in masks:
            query = alpha | (beta << w) | (gamma << (2*w)) | (lam << (3*w)) | (mu << (3*w+1))
            obs = delta << w
            numerator = sum(-1 if ((x & query).bit_count() + (int(y) & obs).bit_count()) % 2 else 1
                            for x, y in enumerate(images))
            expected = Fraction(numerator, len(images))
            result = two_addition_coefficient(perm, alpha, beta, gamma, delta, carry_masks=(lam, mu))
            coefficient_errors += int(result.coefficient != expected)
            negatives += int(expected < 0)
            nonzero += int(expected != 0)
            rows.append(dict(permutation=list(perm), masks=[alpha,beta,gamma,delta,lam,mu],
                             expected=str(expected), actual=str(result.coefficient),
                             graph_nodes=result.graph_nodes, graph_edges=result.graph_edges,
                             contraction_edges=result.contraction_edges, full_map_errors=bad))
        exp.log(f"permutation={perm}; map errors={bad}; cumulative coefficient errors={coefficient_errors}")
    exp.check("P3", map_errors == 0, f"six actual13-qubit maps,8192inputs each; errors={map_errors}")
    exp.check("P4", coefficient_errors == 0 and nonzero > 6,
              f"{len(rows)} coefficients; errors={coefficient_errors}; nonzero={nonzero}, negative={negatives}")
    exp.finish(report_path="out/planar_carry_report.json", rows=rows,
               metadata=dict(width=w, networkx=nx.__version__, exact="integer/Fraction throughout",
                             scope="known planar matchgate/FKT contraction, no novelty or width-separation assertion"))


if __name__ == "__main__":
    run()
