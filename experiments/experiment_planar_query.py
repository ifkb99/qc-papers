"""Does the explicit planar family retain a nonzero arithmetic query?

PREDICTIONS BEFORE EXECUTION: rotated-zigzag permutations at w=8*k*k,
k=1,2,3 are planar. All four masks equal 2**w-1; both incoming carries
are uniformly averaged without a character. Complementing transformed
carry bits reduces Q to (-1)**w * 2**(2-2*w) times the weighted perfect-
matching sum of the architecture with its four path endpoints removed.
Carry edges have weight 2 and interface edges weight 1. For even w the
all-carry matching proves Q >= 2**(-w). At k=1 compare with independently
enumerated cyclic autocorrelation and the generic C86 matchgate method.

Controls: omitting the global sign must fail at odd width 3; all-zero
query masks must defeat the inference that architecture width alone makes
every query nontrivial. No performance advantage or scalar hardness is
predicted. A few exact coefficients do not exclude a further closed form.

Run: uv run --with networkx==3.5 python -m experiments.experiment_planar_query
"""
from fractions import Fraction

import networkx as nx

from lab.harness import Experiment
from lab.planar_carry import planar_perfect_matching, two_addition_coefficient


def rotated_zigzag(k):
    if not isinstance(k, int) or not 1 <= k <= 3:
        raise ValueError("bounded geometry fixture requires k=1,2,3")
    w = 8*k*k
    z = [i//2 if i % 2 == 0 else w-1-i//2 for i in range(w)]
    inverse = {value: i for i, value in enumerate(z)}
    return tuple(inverse[(value+k) % w] for value in z)


def carry_architecture(perm):
    w = len(perm)
    if w < 2 or sorted(perm) != list(range(w)):
        raise ValueError("expected a bit permutation of width at least two")
    graph = nx.Graph()
    graph.add_nodes_from(range(2*w))
    for start in (0, w):
        graph.add_weighted_edges_from((start+i, start+i+1, 2) for i in range(w-1))
    graph.add_weighted_edges_from((source, w+i, 1) for i, source in enumerate(perm))
    return graph


def all_bits_dimer(perm):
    """Restricted all-ones query, using the independently derived reduction."""
    w = len(perm)
    graph = carry_architecture(perm)
    graph.remove_nodes_from((0, w-1, w, 2*w-1))
    matching_sum = planar_perfect_matching(graph)
    coefficient = (-1)**w * matching_sum / (1 << (2*w-2))
    return coefficient, matching_sum, graph


def autocorrelation_query(perm, mask):
    """Independent integer word calculation, restricted to at most eight bits."""
    w, size = len(perm), 1 << len(perm)
    if w > 8:
        raise ValueError("independent autocorrelation exceeds its enumeration budget")
    chars = [(-1)**((mask & x).bit_count()) for x in range(size)]
    correlations = [sum(chars[x]*chars[(x+s) % size] for x in range(size))
                    for s in range(size)]
    twice_g = [correlations[s]+correlations[(s+1) % size] for s in range(size)]
    numerator = 0
    for s in range(size):
        shuffled = sum(((s >> source) & 1) << i for i, source in enumerate(perm))
        numerator += twice_g[s]*twice_g[shuffled]
    return (-1)**mask.bit_count() * Fraction(numerator, 4*size**3)


def run():
    exp = Experiment("planar_query", doc=__doc__)
    exp.predict("P1", "explicit family is planar and direct dimer respects the positive even-width bound")
    exp.predict("P2", "k=1 direct dimer equals independent autocorrelation and generic matchgate contraction")
    exp.must_fail("C1", "discarding the global sign disagrees on a nonzero odd-width query")
    exp.must_fail("C2", "large architecture alone cannot make the all-zero query nontrivial")
    rows = []
    for k in (1, 2, 3):
        perm = rotated_zigzag(k)
        w = len(perm)
        architecture = carry_architecture(perm)
        value, z, reduced = all_bits_dimer(perm)
        planar = nx.check_planarity(architecture)[0]
        positive_bound = value >= Fraction(1, 1 << w)
        rows.append(dict(k=k, width=w, permutation=list(perm), planar=planar,
                         architecture_nodes=len(architecture),
                         reduced_nodes=len(reduced), reduced_edges=reduced.number_of_edges(),
                         matching_sum=str(z), coefficient=str(value),
                         positive_bound=positive_bound,
                         numerator_bits=value.numerator.bit_length(),
                         denominator_bits=value.denominator.bit_length()))
        exp.log(f"k={k}, w={w}, direct graph={len(reduced)} vertices, Q={value}")
    exp.check("P1", all(row["planar"] and row["positive_bound"] for row in rows),
              "three exact rows; no empirical width or memory claim")

    perm = rotated_zigzag(1)
    mask = (1 << len(perm))-1
    direct = Fraction(rows[0]["coefficient"])
    reference = autocorrelation_query(perm, mask)
    generic = two_addition_coefficient(perm, mask, mask, mask, mask)
    exp.check("P2", direct == reference == generic.coefficient,
              f"dimer={direct}; independent={reference}; generic={generic.coefficient}")

    odd_perm = (0, 1, 2)
    odd, _, _ = all_bits_dimer(odd_perm)
    odd_reference = autocorrelation_query(odd_perm, 7)
    exp.fail_check("C1", odd == odd_reference and odd < 0 and abs(odd) != odd_reference,
                   f"signed={odd}; removing global sign={abs(odd)}")
    trivial = two_addition_coefficient(perm, 0, 0, 0, 0).coefficient
    exp.fail_check("C2", trivial == 1 and reference not in (0, 1),
                   f"same architecture: all-zero query={trivial}, all-bit query={reference}")
    exp.finish(report_path="out/planar_query_report.json", rows=rows,
               metadata=dict(exact="integer/Fraction", networkx=nx.__version__,
                             independent_width=8, generic_gadget_nodes=generic.graph_nodes,
                             global_sign_control=str(odd),
                             scope="restricted direct dimer reduction; no practical advantage or hardness claim"))


if __name__ == "__main__":
    run()
