#!/usr/bin/env python3
"""
test_fingerprint_compare.py — 对比 Cython 指纹计算与 Octave 原版

用同一个真实结构，分别用 Cython 和 Python(Octave逻辑) 计算，
验证 order, fing, atom_fing 是否一致。
"""

import sys
import os
import numpy as np
sys.path.insert(0, '/home/feng/uspex_tnt/uspex_gp/FunctionFolder/sys')
import uspex_fast_core as ufc


def octave_makeMatrices(lattice, coordinates, numIons, atomType, Rmax=12.0):
    """Octave makeMatrices.m 的 Python 忠实翻译"""
    species = len(numIons)
    natom = sum(numIons)
    N = np.zeros(natom)
    atomType1 = np.zeros(natom)
    typ_i = np.zeros(natom, dtype=int)
    typ_j = np.array([], dtype=int)
    Nfull = 0
    dist_matrix = np.zeros((0, natom))
    V = abs(np.linalg.det(lattice))

    c = 0
    for i in range(species):
        for j in range(numIons[i]):
            N[c] = numIons[i]
            atomType1[c] = atomType[i]
            typ_i[c] = i + 1  # 1-based
            c += 1

    signum = np.array([[1,1,1],[-1,1,1],[1,-1,1],[1,1,-1],[-1,-1,1],[1,-1,-1],[-1,1,-1],[-1,-1,-1]])
    condition = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,0,1],[1,1,1]])

    vect = np.zeros((13, 3))
    vect[0] = lattice[0]; vect[1] = lattice[1]; vect[2] = lattice[2]
    vect[3] = vect[0]+vect[1]; vect[4] = vect[0]-vect[1]
    vect[5] = vect[0]+vect[2]; vect[6] = vect[0]-vect[2]
    vect[7] = vect[2]+vect[1]; vect[8] = vect[2]-vect[1]
    vect[9] = vect[0]+vect[1]+vect[2]; vect[10] = vect[0]+vect[1]-vect[2]
    vect[11] = vect[0]-vect[1]+vect[2]; vect[12] = -vect[0]+vect[1]+vect[2]
    abs_vect = np.sqrt(np.sum(vect**2, axis=1))

    lengthX = int(np.ceil((Rmax + max(abs_vect)) / min(abs_vect))) + 1

    for i in range(lengthX + 1):
        quit_x = True
        for j in range(lengthX + 1):
            quit_y = True
            for k in range(lengthX + 1):
                quit_z = True
                for quad in range(8):
                    if condition[quad,0]*(i==0) + condition[quad,1]*(j==0) + condition[quad,2]*(k==0) != 0:
                        continue
                    for cc in range(natom):
                        distances = np.zeros(natom)
                        marker = 0
                        for bc in range(natom):
                            x = coordinates[cc,0] + signum[quad,0]*i - coordinates[bc,0]
                            y = coordinates[cc,1] + signum[quad,1]*j - coordinates[bc,1]
                            z = coordinates[cc,2] + signum[quad,2]*k - coordinates[bc,2]
                            Rij = (x*lattice[0,0]+y*lattice[1,0]+z*lattice[2,0])**2
                            Rij += (x*lattice[0,1]+y*lattice[1,1]+z*lattice[2,1])**2
                            Rij += (x*lattice[0,2]+y*lattice[1,2]+z*lattice[2,2])**2
                            if Rij < Rmax**2:
                                quit_z = quit_y = quit_x = False
                                if marker == 0 and Nfull >= natom:
                                    N = np.append(N, N[cc])
                                marker = 1
                                distances[bc] = np.sqrt(Rij)
                        if marker:
                            if dist_matrix.shape[0] == 0:
                                dist_matrix = distances.reshape(1, -1)
                            else:
                                dist_matrix = np.vstack([dist_matrix, distances])
                            typ_j = np.append(typ_j, typ_i[cc])
                            Nfull_dist = dist_matrix.shape[0]
                if quit_z:
                    break
            if quit_y:
                break
        if quit_x:
            break

    return N, V, dist_matrix, typ_i, typ_j


def octave_fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons,
                             Rmax=12.0, sigma=0.05, delta=0.08):
    """Octave fingerprint_calc.m 的 Python 忠实翻译"""
    from scipy.special import erf
    species = len(numIons)
    Ncell = sum(numIons)
    Nfull = dist_matrix.shape[0]
    normaliser = 1.0

    if Ncell == 0:
        return np.array([]), np.array([]), np.array([])

    numBins = int(round(Rmax / delta))
    fing = np.zeros((species * species, numBins))
    fing[:, 0] = -1.0
    order = np.zeros(Ncell)
    atom_fing = np.zeros((Ncell, species, numBins))
    sigm = sigma / np.sqrt(2 * np.log(2))

    for bins in range(2, numBins):  # 0-based, Octave 2:numBins
        for i in range(Ncell):
            for j in range(Nfull):
                d = dist_matrix[j, i]
                if d > 0 and abs(d - delta * (bins - 0.5)) < 4 * sigm:
                    R0 = d
                    # interval(2) = upper
                    val_upper = delta * bins - R0
                    if abs(val_upper / np.sqrt(2)) <= 5 * sigm:
                        upper = val_upper / (np.sqrt(2) * sigm)
                    else:
                        upper = 5 * np.sign(val_upper)
                    # interval(1) = lower
                    val_lower = delta * (bins - 1) - R0
                    if abs(val_lower / np.sqrt(2)) <= 5 * sigm:
                        lower = val_lower / (np.sqrt(2) * sigm)
                    else:
                        lower = 5 * np.sign(val_lower)

                    delt = 0.5 * (erf(upper) - erf(lower))
                    atom_fing[i, typ_j[j] - 1, bins] += delt / (Ni[j] * R0**2)
                    fing[(typ_i[i] - 1) * species + (typ_j[j] - 1), bins] += delt / R0**2

            atom_fing[i, :, bins] = V * atom_fing[i, :, bins] / (4 * np.pi * delta) - normaliser
            for j in range(species):
                weight = numIons[j] / sum(numIons)
                order[i] += weight * delta * atom_fing[i, j, bins]**2 / (V / Ncell)**(1.0/3.0)

        for i in range(species):
            for j in range(species):
                fing[(i-1+1)*species + j, bins] = V * fing[(i)*species + j, bins] / (4 * np.pi * numIons[i] * numIons[j] * delta) - normaliser

    order = np.sqrt(order)
    return order, fing, atom_fing


def test_compare():
    print("=== 指纹计算对比测试 ===")

    # 构造测试结构 (8个分子质心)
    np.random.seed(42)
    lattice = np.eye(3) * 15.0 + np.random.randn(3, 3) * 0.5
    coords = np.random.rand(8, 3)
    numIons = np.array([4, 4], dtype=np.int32)  # 2 种分子
    atomType = np.array([1, 2], dtype=np.int32)

    Rmax = 12.0; sigma = 0.05; delta = 0.08

    # Octave 逻辑
    print("  运行 Octave 逻辑 (Python 翻译)...")
    Ni_oct, V_oct, dist_oct, ti_oct, tj_oct = octave_makeMatrices(
        lattice, coords, numIons, atomType, Rmax)
    order_oct, fing_oct, af_oct = octave_fingerprint_calc(
        Ni_oct, V_oct, dist_oct, ti_oct, tj_oct, numIons, Rmax, sigma, delta)

    # Cython
    print("  运行 Cython...")
    result = ufc.compute_all(coords, lattice, numIons, atomType, Rmax, sigma, delta, 3)
    order_cy = result['order']
    fing_cy = result['fing']
    af_cy = result['atom_fing']

    # 对比
    print(f"\n  order: Octave shape={order_oct.shape}, Cython shape={order_cy.shape}")
    if order_oct.shape == order_cy.shape:
        diff = np.abs(order_oct - order_cy)
        print(f"  order max diff: {diff.max():.6e}")
        print(f"  order: Octave={order_oct[:4]}")
        print(f"  order: Cython={order_cy[:4]}")

    print(f"\n  fing: Octave shape={fing_oct.shape}, Cython shape={fing_cy.shape}")
    if fing_oct.shape == fing_cy.shape:
        diff = np.abs(fing_oct - fing_cy)
        print(f"  fing max diff: {diff.max():.6e}")
        # 看非零列
        nonzero_cols = np.where(np.any(np.abs(fing_oct) > 1e-10, axis=0))[0]
        if len(nonzero_cols) > 0:
            c = nonzero_cols[0]
            print(f"  fing col {c}: Octave={fing_oct[:, c]}")
            print(f"  fing col {c}: Cython={fing_cy[:, c]}")

    print(f"\n  atom_fing: Octave shape={af_oct.shape}, Cython shape={af_cy.shape}")
    if af_oct.shape == af_cy.shape:
        diff = np.abs(af_oct - af_cy)
        print(f"  atom_fing max diff: {diff.max():.6e}")

    print(f"\n  V: Octave={V_oct:.6f}, Cython={result['V']:.6f}")
    print(f"  n_pairs (Cython): {result['n_pairs']}")
    print(f"  Nfull (Octave): {dist_oct.shape[0]}")


if __name__ == '__main__':
    test_compare()
