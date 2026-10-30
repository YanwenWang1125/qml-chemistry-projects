"""E2.0 step 2: self-test of COLOR_17 (supplementary plan, E2.0 step 2).

Pass: generators commute pairwise, GF(2) rank 16 (k = 1); no nontrivial logical of weight <= 4
and one of weight 5 (d = 5); the stated logicals X17, Z17 commute with all generators and
anticommute with each other.
Run:  python test_color17.py
"""
import numpy as np

from ftrs.codes import COLOR_17, _matrix, css_distance_check, gf2_rank

c = COLOR_17
n = c.n_data
hx, hz = _matrix(c.x_supports, n), _matrix(c.z_supports, n)
ok = True


def report(name, cond):
    global ok
    ok &= bool(cond)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")


report("X and Z generators commute pairwise (Hx Hz^T = 0 mod 2)", not (hx @ hz.T % 2).any())
rank = gf2_rank(hx) + gf2_rank(hz)
report(f"GF(2) rank {gf2_rank(hx)} + {gf2_rank(hz)} = {rank}, k = {n - rank}", rank == 16)
lx, lz = _matrix(c.logical_x, n)[0], _matrix(c.logical_z, n)[0]
report("logical X and Z commute with the Z and X generators",
       not (hz @ lx % 2).any() and not (hx @ lz % 2).any())
report("logical X and Z anticommute", int(lx @ lz) % 2 == 1)
x17 = np.zeros(n, dtype=np.uint8)
x17[16] = 1
print("note: the paper's X_17 anticommutes with Z generators",
      [f"g{i + 1}" for i in np.nonzero(hz @ x17 % 2)[0]], "(so it is not a logical of Eq. (1))")
d = css_distance_check(c, 5)
print("lowest-weight logicals found (0-indexed qubits):", d)
report("distance 5: no logical of weight <= 4, one of weight 5 for both types",
       d["X"] and d["Z"] and d["X"][0] == 5 and d["Z"][0] == 5)
report(f"stabilizer group via Code: {len(c._sx)} X-type x {len(c._sz)} Z-type elements", len(c._sx) == 256)
print("\nALL PASS" if ok else "\nSOME CHECKS FAILED")
