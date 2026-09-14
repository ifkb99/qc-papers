"""Can paired sparse elimination reduce memory for the explicit carry query?

PREDICTIONS BEFORE EXECUTION: all-carry matching pairs give nonzero leading
even Pfaffians under any pair ordering because weights are positive and
each prefix/complement is matched. Natural pair ordering has bandwidth
at most 4*k+2. Sparse Schur elimination equals existing dense FKT and
generic pair-min-fill, with a bounded live entry count. Hypothesis: at
k=3 its traced peak is below the existing dense implementation; no win
over min-fill or known planar nested dissection is assumed.

Fixed family k=1,2,3, w=8*k*k, all-bit query and uniform initial carries.
Fresh subprocess measurements include graph/order/orientation, exact
arithmetic and extraction; Python imports/startup are outside the trace.
Controls reject a nonmatching pair order and an incorrect allowed-edge
deletion. The z^2 coefficient is independently counted by even-gap tests.
This is a known sparse-elimination application, not a novelty claim.

Run: uv run --with networkx==3.5 python -m experiments.experiment_carry_matching
"""
from fractions import Fraction
import gc
import json
import subprocess
import sys
import time
import tracemalloc

import networkx as nx

from experiments.experiment_planar_query import carry_architecture, rotated_zigzag
from lab.harness import Experiment
from lab.planar_carry import _orientation, planar_perfect_matching


def graph_for(k):
    perm = rotated_zigzag(k)
    w = len(perm)
    graph = carry_architecture(perm)
    graph.remove_nodes_from((0, w-1, w, 2*w-1))
    return w, perm, graph


def reference_pairs(w):
    return [(offset+i, offset+i+1)
            for i in range(1, w-1, 2) for offset in (0, w)]


def pair_min_fill(graph, pairs):
    """Known greedy fill heuristic on the graph contracted by matching pairs."""
    owner = {v: i for i, pair in enumerate(pairs) for v in pair}
    quotient = nx.Graph()
    quotient.add_nodes_from(range(len(pairs)))
    quotient.add_edges_from((owner[u], owner[v]) for u, v in graph.edges
                            if owner[u] != owner[v])
    order = []
    while quotient:
        def score(v):
            neighbors = sorted(quotient[v])
            missing = sum(not quotient.has_edge(a, b)
                          for j, a in enumerate(neighbors) for b in neighbors[j+1:])
            return missing, len(neighbors), v
        chosen = min(quotient, key=score)
        neighbors = list(quotient[chosen])
        quotient.add_edges_from((a, b) for j, a in enumerate(neighbors)
                                for b in neighbors[j+1:])
        quotient.remove_node(chosen)
        order.append(pairs[chosen])
    return order


def paired_sparse_matching(graph, pairs):
    """Known exact skew Schur elimination; positive matched prefixes certify pivots."""
    order = [v for pair in pairs for v in pair]
    if len(order) != len(graph) or len(set(order)) != len(order) or set(order) != set(graph):
        raise ValueError("pairs must cover every graph vertex exactly once")
    if any(not graph.has_edge(*pair) for pair in pairs):
        raise ValueError("each prescribed pair must be a graph edge")
    if any(data.get("weight", 1) <= 0 for _, _, data in graph.edges(data=True)):
        raise ValueError("positive weights required for certified nonzero pivots")
    if not nx.check_planarity(graph)[0]:
        raise ValueError("planar graph required")
    position = {v: i for i, v in enumerate(order)}
    orientations = _orientation(graph)
    rows = [dict() for _ in order]
    for (u, v), sign in orientations.items():
        i, j = position[u], position[v]
        value = sign*Fraction(graph[u][v].get("weight", 1))
        if i > j:
            i, j, value = j, i, -value
        rows[i][j] = value
    calibration = 1
    for i in range(0, len(order), 2):
        calibration *= 1 if rows[i][i+1] > 0 else -1
    initial = sum(map(len, rows))
    bandwidth = max((j-i for i, row in enumerate(rows) for j in row), default=0)
    peak, updates, largest_front = initial, 0, 0
    max_fraction_bits = 0
    value = Fraction(1)
    for pivot in range(0, len(order), 2):
        a, b = rows[pivot], rows[pivot+1]
        p = a.get(pivot+1, 0)
        if not p:
            raise ArithmeticError("certified matched-prefix pivot unexpectedly vanished")
        value *= p
        neighbors = sorted((set(a) | set(b)) - {pivot+1})
        largest_front = max(largest_front, len(neighbors))
        for offset, i in enumerate(neighbors):
            ai, bi = a.get(i, 0), b.get(i, 0)
            for j in neighbors[offset+1:]:
                correction = (bi*a.get(j, 0) - ai*b.get(j, 0))/p
                if not correction:
                    continue
                entry = rows[i].get(j, 0)+correction
                updates += 1
                if entry:
                    rows[i][j] = entry
                    max_fraction_bits = max(max_fraction_bits,
                                            entry.numerator.bit_length()+entry.denominator.bit_length())
                else:
                    rows[i].pop(j, None)
        # Count the transient updated front before releasing the pivot rows.
        peak = max(peak, sum(map(len, rows)))
        rows[pivot], rows[pivot+1] = {}, {}
    return calibration*value, dict(initial_entries=initial, peak_entries=peak,
                                    bandwidth=bandwidth, largest_front=largest_front,
                                    updates=updates, max_updated_fraction_bits=max_fraction_bits)


def measure_case(k, mode):
    gc.collect()
    tracemalloc.start()
    start = time.perf_counter()
    w, _, graph = graph_for(k)
    stats = {}
    if mode == "dense":
        z = planar_perfect_matching(graph)
    else:
        pairs = reference_pairs(w)
        if mode == "min_fill":
            pairs = pair_min_fill(graph, pairs)
        elif mode != "banded":
            raise ValueError("unknown case")
        z, stats = paired_sparse_matching(graph, pairs)
    coefficient = z/(1 << (2*w-2))
    elapsed = time.perf_counter()-start
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return dict(k=k, width=w, mode=mode, matching_sum=str(z), coefficient=str(coefficient),
                elapsed_seconds=elapsed, traced_peak_bytes=peak, **stats)


def two_cross_count(graph, w):
    """Independent interface-subset count: remaining path gaps must be even."""
    cross = [(u, v-w) for u, v in graph.edges if u < w <= v]
    def valid(a, b):
        lo, hi = sorted((a, b))
        return (lo-1) % 2 == 0 and (hi-lo-1) % 2 == 0 and (w-2-hi) % 2 == 0
    return sum(valid(a, c) and valid(b, d)
               for offset, (a, b) in enumerate(cross) for c, d in cross[offset+1:])


def run():
    exp = Experiment("carry_matching", doc=__doc__)
    exp.predict("P1", "dense, banded, and generic pair-min-fill give identical exact scalars")
    exp.predict("P2", "natural paired bandwidth and live sparse-entry bounds hold without pivot search")
    exp.predict("P3", "at k=3, natural sparse traced peak is below the existing dense evaluator")
    exp.predict("P4", "two-interface coefficient equals binomial(4k^2-2k,2)")
    exp.must_fail("C1", "an arbitrary nonmatching pair order must be rejected")
    exp.must_fail("C2", "deleting an allowed interface edge must change the matching sum")
    rows, counts = [], []
    for k in (1, 2, 3):
        for mode in ("dense", "banded", "min_fill"):
            child = subprocess.run([sys.executable, "-m", "experiments.experiment_carry_matching",
                                    "--case", str(k), mode], capture_output=True, text=True,
                                   check=True, timeout=60)
            row = json.loads(child.stdout)
            rows.append(row)
            exp.log(f"k={k} {mode}: peak={row['traced_peak_bytes']} bytes, "
                    f"seconds={row['elapsed_seconds']:.4f}, Q={row['coefficient']}")
        w, _, graph = graph_for(k)
        count = two_cross_count(graph, w)
        t = 4*k*k-2*k
        counts.append(dict(k=k, observed=count, expected=t*(t-1)//2))
    exp.check("P1", all(len({row["coefficient"] for row in rows if row["k"] == k}) == 1
                        for k in (1, 2, 3)), "three modes on each of three fixed-family instances")
    banded = [row for row in rows if row["mode"] == "banded"]
    exp.check("P2", all(row["bandwidth"] <= 4*row["k"]+2 and
                         row["peak_entries"] <= row["initial_entries"]+(row["bandwidth"]+2)**2
                         for row in banded), "exact matched pivots and structural upper-entry bound")
    high = {row["mode"]: row for row in rows if row["k"] == 3}
    exp.check("P3", high["banded"]["traced_peak_bytes"] < high["dense"]["traced_peak_bytes"],
              f"dense={high['dense']['traced_peak_bytes']}, banded={high['banded']['traced_peak_bytes']}, "
              f"min-fill={high['min_fill']['traced_peak_bytes']}; no best-known-method claim")
    exp.check("P4", all(row["observed"] == row["expected"] for row in counts), str(counts))
    w, _, graph = graph_for(1)
    pairs = reference_pairs(w)
    wrong = pairs.copy()
    wrong[0], wrong[1] = (pairs[0][0], pairs[1][0]), (pairs[0][1], pairs[1][1])
    rejected = False
    try:
        paired_sparse_matching(graph, wrong)
    except ValueError:
        rejected = True
    exp.fail_check("C1", rejected, "nonadjacent first pair cannot certify a pivot")
    changed = graph.copy()
    changed.remove_edge(1, 11)
    original = planar_perfect_matching(graph)
    altered = planar_perfect_matching(changed)
    exp.fail_check("C2", original != altered, f"original={original}, incorrect deletion={altered}")
    exp.finish(report_path="out/carry_matching_report.json", rows=rows,
               metadata=dict(two_cross_coefficients=counts, networkx=nx.__version__,
                             exact="integer/Fraction", measurement="fresh subprocess per case; Python traced allocations",
                             excluded="Python/import startup and process-wide native/RSS memory",
                             scope="known sparse elimination on selected arithmetic scalar; no novelty or optimality claim"))


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "--case":
        print(json.dumps(measure_case(int(sys.argv[2]), sys.argv[3])))
    else:
        run()
