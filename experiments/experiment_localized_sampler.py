"""Does finite basis support remove the exponential early-prefix vector?

PREDICTIONS WRITTEN BEFORE MEASUREMENT (TODO 19).
P1: Every joint (k,y) and marginal agrees with the existing SparseOrbitPrefix
    and spectral effects, for r=6,t=8,theta=pi/2 with only s=1..7 changing.
P2: Enumerating the component proposal and its acceptance reproduces each
    early conditional distribution; mean cost over k is <=5(d+1), not a
    pointwise bound. All progression sampler branches implement gcd lifting.
P3: Actual sampled outputs at t=63 with varying split, and a separately
    supplied large abstract period, require neither orbit nor prefix arrays.
C1: Deleted interference fails the normalized conditional distribution test.
C2: Uniform final eigenphase weights fail the marginal output test.
C3: A near-cancelled phase violates a proposed pointwise 5(d+1) bound.

Known order and orbit indices are inputs. No order-discovery or compiled
large-modulus simulation is claimed. Dense references are capped BEFORE
allocation at t<=8,r<=8. A sampled path's probability is joint, not marginal.
Run with Python 3.12, numpy<2.5, OPENBLAS_NUM_THREADS=1, as in the other probes.
"""
from __future__ import annotations
import math
import time
from datetime import datetime, timezone
import numpy as np
from lab import Experiment
from lab.localized import LocalizedOrbitDefect
from lab.fourier_sampling import interval_path, progression_sample
from lab.prefix import SparseOrbitPrefix
from lab.spectral import single_defect_effects


def mixer(theta=np.pi/2):
    c, s = np.cos(theta/2), np.sin(theta/2)
    return np.array([[c, -1j*s], [-1j*s, c]])


def tiny_reference(defect, width, split):
    r = defect.period
    if not 0 <= width <= 8 or not 1 <= r <= 8 or not 0 <= split <= width:
        raise ValueError("reference dimensions exceed pre-allocation budget")
    V = np.eye(r, dtype=complex)
    V[np.ix_(defect.indices, defect.indices)] = defect.unitary
    def column(l):
        v = V[:, l % r]
        ids = np.flatnonzero(v != 0)
        return ids, v[ids]
    prefix = SparseOrbitPrefix(r, split, column, defect.D)
    F = np.exp(-2j*np.pi*np.arange(r)[:, None]*np.arange(r)[None, :]/r)/np.sqrt(r)
    E = single_defect_effects(r, width, split, F.conj().T @ V @ F)
    return prefix, E.sum(axis=(1, 2)).real/r


class ReplayBits:
    """Traverse each positive-probability proposal path without Monte Carlo."""
    def __init__(self, output, lift):
        self.output, self.lift, self.step = output, lift, 0
    def random(self):
        bit = (self.output >> self.step) & 1
        self.step += 1
        return np.nextafter(1., 0.) if bit else 0.
    def integers(self, high):
        assert 0 <= self.lift < high
        return self.lift


def progression_proposal(split, start, stride, count, q, Q):
    L = 1 << split
    if split > 8:
        raise ValueError("proposal enumeration cap")
    g = math.gcd(stride, L)
    T, beta = L//g, stride//g
    p = np.zeros(L)
    branch_error = 0
    for w in range(T):
        path = interval_path(T.bit_length()-1, count, stride*q, Q, output=w)
        if path["path_probability"] == 0:
            continue
        residue = (pow(beta, -1, T)*w) % T if T > 1 else 0
        for lift in range(g):
            generated = progression_sample(split, start, stride, count, q, Q, ReplayBits(w, lift))
            expected = residue+T*lift
            branch_error = max(branch_error, abs(generated["output"]-expected))
            p[generated["output"]] += path["path_probability"]/g
    return p, branch_error


def envelope_audit(defect, width, split, prefix):
    L, H, Q = 1 << split, 1 << (width-split), 1 << width
    worst = incoherent_error = branch_error = mean = 0.
    for k in range(defect.period):
        data = defect._phase_data(split, k)
        if data["norm"] == 0:
            continue
        mean += (data["norm"]/defect.period)*defect.m*data["T"]/data["norm"]
        # One fixed q per k; the proof covers all q, while this independently
        # checks feedback, proposal normalization and rejection for each row.
        q = min(1, H-1)
        proposal = data["weights"][0]/data["T"]*np.array([
            interval_path(split, L, q*defect.period-k*Q, Q*defect.period,
                          output=z)["path_probability"] for z in range(L)])
        for ip, (p, n) in enumerate(zip(defect.indices, data["counts"])):
            if data["weights"][ip+1] == 0:
                continue
            component, error = progression_proposal(split, p, defect.period, n, q, Q)
            branch_error = max(branch_error, error)
            proposal += data["weights"][ip+1]/data["T"]*component
        accepted = np.zeros(L)
        incoherent = np.zeros(L)
        for z in range(L):
            b = defect._components(width, data, q+H*z)
            mag, diagonal = float(abs(b.sum())**2), float(np.vdot(b, b).real)
            accepted[z] = proposal[z]*mag/(defect.m*diagonal) if diagonal > 0 else 0.
            incoherent[z] = diagonal
        row = prefix.phase_row(k)
        reference = np.abs(np.fft.fft(row*np.exp(-2j*np.pi*q*np.arange(L)/Q)))**2
        reference /= float(np.vdot(row, row).real)*L
        expected_acceptance = data["norm"]/(defect.m*data["T"])
        worst = max(worst, abs(proposal.sum()-1), abs(accepted.sum()-expected_acceptance),
                    float(np.max(np.abs(accepted/expected_acceptance-reference))))
        incoherent /= incoherent.sum()
        incoherent_error = max(incoherent_error, float(np.max(np.abs(incoherent-reference))))
    return dict(error=worst, branch_error=int(branch_error), mean_early_attempts=mean,
                bound=5*defect.m, deleted_interference_error=incoherent_error)


def main():
    exp = Experiment("localized_sampler", doc=__doc__)
    exp.predict("P1", "localized joint/marginal distributions equal existing prefix and spectral routes")
    exp.predict("P2", "normalized proposal/rejection and enumerated progression sampling paths agree")
    exp.predict("P3", "wide actual samples have fixed matrix payload and no exponential arrays")
    exp.must_fail("C1", "deleting interference changes normalized early distributions")
    exp.must_fail("C2", "uniform final eigenphases change the marginal")
    exp.must_fail("C3", "conditional expected rejection count has no uniform 5m bound")
    started, rows = time.perf_counter(), []
    max_interference = max_uniform = 0.
    defect = LocalizedOrbitDefect(6, (0, 5), mixer())
    for split in range(1, 8):
        width = 8
        prefix, spectral = tiny_reference(defect, width, split)
        actual = np.array([[defect.forced_joint(width, split, k, y)["joint_latent_output_probability"]
                            for y in range(1 << width)] for k in range(6)])
        expected = np.array([[prefix.forced_joint(width, k, y)["joint_latent_output_probability"]
                              for y in range(1 << width)] for k in range(6)])
        weights = np.array([defect.phase_probability(split, k) for k in range(6)])
        error = float(np.max(np.abs(actual-expected)))
        spectral_error = float(np.max(np.abs(actual.sum(axis=0)-spectral)))
        exp.check("P1", max(error, spectral_error, abs(actual.sum()-1)) < 2e-11,
                  f"s={split}: joint={error:.2e}, spectral={spectral_error:.2e}")
        audit = envelope_audit(defect, width, split, prefix)
        exp.check("P2", audit["error"] < 2e-11 and audit["branch_error"] == 0
                  and audit["mean_early_attempts"] <= audit["bound"]+1e-12,
                  f"s={split}: {audit}")
        max_interference = max(max_interference, audit["deleted_interference_error"])
        max_uniform = max(max_uniform, float(np.max(np.abs((actual/weights[:, None]).mean(axis=0)-spectral))))
        rows.append(dict(kind="fixed_sweep", width=width, split=split,
                         joint_error=error, spectral_error=spectral_error, envelope=audit))

    # Separate fixed edge cases, not presented as a one-parameter scaling law.
    rng = np.random.default_rng(49019)
    z = rng.normal(size=(3, 3))+1j*rng.normal(size=(3, 3))
    complex_block = np.linalg.qr(z)[0]
    for edge, width, split in ((defect, 5, 0), (defect, 5, 5),
                               (LocalizedOrbitDefect(5, (0, 2, 4), complex_block), 5, 3),
                               (LocalizedOrbitDefect(1, (), np.empty((0, 0))), 0, 0)):
        prefix, spectral = tiny_reference(edge, width, split)
        marginal = np.array([sum(edge.forced_joint(width, split, k, y)["joint_latent_output_probability"]
                                for k in range(edge.period)) for y in range(1 << width)])
        error = float(np.max(np.abs(marginal-spectral)))
        exp.check("P1", error < 2e-11, f"edge r={edge.period}, t={width}, s={split}: {error:.2e}")
        rows.append(dict(kind="edge", period=edge.period, width=width, split=split, error=error))

    for period, split in ((6, 2), (6, 16), (6, 32), (6, 62), (1_000_000_007, 32)):
        wide = LocalizedOrbitDefect(period, (0, period-1), mixer())
        samples = [wide.sample(63, split, rng) for _ in range(16)]
        exp.check("P3", all(0 <= x["output"] < 1 << 63
                            and x["joint_latent_output_probability"] > 0 for x in samples)
                  and wide.stats()["matrix_payload_bytes"] == 64,
                  f"r={period}, t=63, s={split}: sampled 16 paths, matrix=64 bytes")
        rows.append(dict(kind="wide", period=period, width=63, split=split, samples=samples,
                         stats=wide.stats(), abstract_known_period=True))

    angle = np.pi/4-1e-4
    c, s = np.cos(angle), np.sin(angle)
    near = LocalizedOrbitDefect(3, (0, 1), np.array([[c, s], [-s, c]]))
    cost = near.forced_joint(3, 0, 0, 0)["expected_early_attempts"]
    exp.fail_check("C1", max_interference > .01, f"normalized deleted-interference error={max_interference:.6g}")
    exp.fail_check("C2", max_uniform > .001, f"uniform-weight marginal error={max_uniform:.6g}")
    exp.fail_check("C3", cost > 1000*5*near.m, f"near-zero phase conditional expected attempts={cost:.6g}")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    exp.finish(report_path=f"out/localized_sampler_{stamp}.json", rows=rows,
               metadata=dict(numpy=np.__version__, elapsed_seconds=time.perf_counter()-started,
                             reference_cap="t<=8,r<=8; spectral entries<=16384",
                             precision="float64, no rare-output or bit-complexity guarantee"))


if __name__ == "__main__":
    main()
