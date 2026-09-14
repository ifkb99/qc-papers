"""Experimental exact carry contraction via known planar matchgate machinery.

Requires NetworkX (the bounded experiment pins3.5). No general CNOT or
nonplanar simulation claim. The API takes an explicit two-adder word map,
not an arbitrary Circuit certificate. All sums and signs use Fractions.
"""
from dataclasses import dataclass
from fractions import Fraction
from operator import index

import networkx as nx


class NonPlanarCarry(ValueError):
    pass


def pfaffian(matrix):
    """Exact skew elimination, including simultaneous row/column pivots."""
    a = [[Fraction(x) for x in row] for row in matrix]
    n = len(a)
    if any(len(row) != n for row in a):
        raise ValueError("expected square matrix")
    if n % 2:
        return Fraction(0)
    result = Fraction(1)
    for k in range(0, n, 2):
        pivot = next((j for j in range(k+1, n) if a[k][j]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != k+1:
            a[k+1], a[pivot] = a[pivot], a[k+1]
            for row in a:
                row[k+1], row[pivot] = row[pivot], row[k+1]
            result = -result
        p = a[k][k+1]
        result *= p
        for i in range(k+2, n):
            for j in range(i+1, n):
                a[i][j] += (a[k+1][i]*a[k][j] - a[k][i]*a[k+1][j])/p
                a[j][i] = -a[i][j]
    return result


def _orientation(graph):
    """Solve odd face-orientation constraints in each planar component."""
    edges = sorted(tuple(sorted(e)) for e in graph.edges)
    positions = {e: i for i, e in enumerate(edges)}
    equations = []
    for nodes in nx.connected_components(graph):
        component = graph.subgraph(nodes)
        planar, embedding = nx.check_planarity(component)
        if not planar:
            raise NonPlanarCarry("weighted matching graph is nonplanar")
        visited, faces = set(), []
        for u, v in embedding.edges():
            if (u, v) not in visited:
                faces.append(embedding.traverse_face(u, v, visited))
        # Any face may serve as the outer face. Bridges occur twice with
        # opposite directions; their variable cancels but their constant does not.
        for face in faces[1:]:
            mask, rhs = 0, 1
            for u, v in zip(face, face[1:]+face[:1]):
                mask ^= 1 << positions[tuple(sorted((u, v))) ]
                rhs ^= int(u > v)
            equations.append((mask, rhs))
    basis = {}
    for mask, rhs in equations:
        while mask:
            p = (mask & -mask).bit_length()-1
            if p not in basis:
                basis[p] = (mask, rhs)
                break
            row, value = basis[p]
            mask, rhs = mask ^ row, rhs ^ value
        else:
            if rhs:
                raise ArithmeticError("inconsistent planar face orientation")
    solution = 0
    for p in sorted(basis, reverse=True):
        row, rhs = basis[p]
        if rhs ^ ((row & solution).bit_count() & 1):
            solution |= 1 << p
    return {e: 1 if (solution >> i) & 1 else -1 for e, i in positions.items()}


def planar_perfect_matching(graph, *, max_nodes=1024):
    """Known FKT algorithm; calibrate its global sign with one matching.

    Signed weights are retained. Planarity is checked on nonzero edges;
    an ordinary crossing is never replaced by a fermionic swap.
    """
    if len(graph) > max_nodes:
        raise ValueError("matching graph exceeds node budget")
    g = nx.Graph()
    g.add_nodes_from(graph.nodes)
    for u, v, data in graph.edges(data=True):
        if u == v:
            raise ValueError("self-loop unsupported")
        weight = Fraction(data.get("weight", 1))
        if weight:
            g.add_edge(u, v, weight=weight)
    # Relabeling permits arbitrary hashable input node labels.
    g = nx.convert_node_labels_to_integers(g)
    if not nx.check_planarity(g)[0]:
        raise NonPlanarCarry("weighted matching graph is nonplanar")
    n = len(g)
    if n % 2:
        return Fraction(0)
    matching = nx.max_weight_matching(g, maxcardinality=True, weight=None)
    if len(matching)*2 != n:
        return Fraction(0)
    orientation = _orientation(g)
    a = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    for (u, v), sign in orientation.items():
        a[u][v] = sign*g[u][v]["weight"]
        a[v][u] = -a[u][v]
    pairs = sorted(tuple(sorted(e)) for e in matching)
    sequence = [q for pair in pairs for q in pair]
    inversions = sum(sequence[i] > sequence[j] for i in range(n) for j in range(i+1, n))
    calibration = -1 if inversions % 2 else 1
    for pair in pairs:
        calibration *= orientation[pair]
    return calibration*pfaffian(a)


def hadamard(values):
    """Unnormalized local Walsh transform; each joined edge later costs1/2."""
    return [sum(value*(-1 if (mask & bits).bit_count() % 2 else 1)
                for bits, value in enumerate(values)) for mask in range(len(values))]


class _Gadgets:
    def __init__(self):
        self.graph = nx.Graph()
        self.scale = Fraction(1)

    def node(self):
        n = len(self.graph)
        self.graph.add_node(n)
        return n

    def edge(self, a, b, weight=1):
        if weight:
            self.graph.add_edge(a, b, weight=Fraction(weight))

    def signature(self, values):
        """Realize pure-parity unary/ternary signatures as planar gadgets.

        A bit1 means deletion of that external vertex. For three even ports,
        a center plus the external triangle realizes all pair-deletion values.
        Adding a leaf flips one external bit and gives odd parity.
        """
        f = [Fraction(x) for x in values]
        arity = len(f).bit_length()-1
        if len(f) not in (2, 8):
            raise ValueError("only unary or ternary signatures supported")
        parities = {i.bit_count() % 2 for i, x in enumerate(f) if x}
        if len(parities) > 1:
            raise ValueError("signature has mixed parity")
        if not parities:
            self.scale = Fraction(0)
            return [self.node() for _ in range(arity)]
        if arity == 1:
            port = self.node()
            if f[0]:
                self.edge(port, self.node(), f[0])
            else:
                self.scale *= f[1]
            return [port]
        odd = next(iter(parities)) == 1
        if odd:
            ports = self.signature([f[i ^ 1] for i in range(8)])
            leaf = self.node()
            self.edge(ports[0], leaf)
            ports[0] = leaf
            return ports
        ports = [self.node() for _ in range(3)]
        pivot = next((j for j in range(3) if f[7 ^ (1 << j)]), None)
        if pivot is None:
            # Only f000 survives: three private forced edges, including zeros.
            for j, port in enumerate(ports):
                self.edge(port, self.node(), f[0] if j == 0 else 1)
            return ports
        center = self.node()
        for j, port in enumerate(ports):
            self.edge(center, port, f[7 ^ (1 << j)])
        other = [j for j in range(3) if j != pivot]
        self.edge(ports[other[0]], ports[other[1]], f[0]/f[7 ^ (1 << pivot)])
        return ports


def carry_factors(alpha, beta, gamma, delta):
    """Unnormalized K(p,q,s), R(r,t,s), little-endian tensor indices."""
    k, r = [0]*8, [0]*8
    for p in (0, 1):
        for a in (0, 1):
            for b in (0, 1):
                total = a+b+p
                k[p | ((total >> 1) << 1) | ((total & 1) << 2)] += (-1)**(alpha*a+beta*b)
    for carry in (0, 1):
        for s in (0, 1):
            for c in (0, 1):
                total = s+c+carry
                r[carry | ((total >> 1) << 1) | (s << 2)] += (-1)**(gamma*c+delta*(total & 1))
    return k, r


@dataclass(frozen=True)
class PlanarCarryResult:
    coefficient: Fraction
    graph_nodes: int
    graph_edges: int
    contraction_edges: int


def two_addition_coefficient(permutation, alpha, beta, gamma, delta, *,
                             carry_masks=(0, 0), max_width=64):
    """Exact E[(-1)^(alpha.A+beta.B+gamma.C+delta.Y+lambda.c+mu.d)].

    S=(A+B+c) mod2^w, Y=(P(S)+C+d) mod2^w, with every input uniform.
    P(S)_i=S_(permutation[i]); independent c,d remain unchanged. Two
    unqueried output-carry XOR bits average out and are allowed physically.
    Only the explicit word map is assumed; circuit recognition is additional.
    Refuses nonplanar nonzero matching graphs; no treewidth claim is made.
    """
    perm = tuple(index(q) for q in permutation)
    w = len(perm)
    if not 1 <= w <= index(max_width) or sorted(perm) != list(range(w)):
        raise ValueError("expected bounded bit permutation")
    masks = tuple(index(q) for q in (alpha, beta, gamma, delta))
    if any(q < 0 or q.bit_length() > w for q in masks):
        raise ValueError("word mask outside width")
    if len(carry_masks) != 2 or any(q not in (0, 1) for q in carry_masks):
        raise ValueError("expected two carry query bits")
    gadgets, first, second = _Gadgets(), [], []
    for i in range(w):
        k, r = carry_factors(*((q >> i) & 1 for q in masks))
        first.append(gadgets.signature(hadamard(k)))
        second.append(gadgets.signature(hadamard(r)))
    joined = 0
    for chain, carry_mask in zip((first, second), carry_masks):
        left = gadgets.signature(hadamard([1, (-1)**carry_mask]))[0]
        right = gadgets.signature(hadamard([1, 1]))[0]
        gadgets.edge(left, chain[0][0])
        gadgets.edge(chain[-1][1], right)
        joined += 2
        for i in range(w-1):
            gadgets.edge(chain[i][1], chain[i+1][0])
            joined += 1
    for i, source in enumerate(perm):
        gadgets.edge(first[source][2], second[i][2])
        joined += 1
    graph = gadgets.graph
    value = (gadgets.scale*planar_perfect_matching(graph)
             if gadgets.scale else Fraction(0))
    # 2^joined undoes H on both ends of every edge; 2^(3w+2)
    # averages over A,B,C and both independent initial carry bits.
    value /= 1 << (joined+3*w+2)
    return PlanarCarryResult(value, len(graph), graph.number_of_edges(), joined)
