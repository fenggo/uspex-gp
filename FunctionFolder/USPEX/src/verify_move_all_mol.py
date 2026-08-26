#!/usr/bin/env python3
"""验证 move_all_mol_Mutation 向量化版本与原版数学等价"""
import numpy as np

def move_all_mol_original(ZMATRIX_centers, LATTICE, order, max_sigma):
    """原版逐分子循环"""
    N = len(ZMATRIX_centers)
    new_Coord = ZMATRIX_centers.copy()
    new_Coord = new_Coord @ np.linalg.inv(LATTICE)  # Cartesian → fractional

    if len(LATTICE.shape) == 1 or LATTICE.shape[0] == 6:
        lat = LATTICE if LATTICE.shape[0] == 3 else LATTICE[:3].reshape(3, 3)
    else:
        lat = LATTICE
    # 简化: 直接用对角线
    temp_potLat = np.array([np.linalg.norm(lat[0]),
                             np.linalg.norm(lat[1]),
                             np.linalg.norm(lat[2])])

    ranking = np.argsort(order)
    r1 = order[ranking[0]]
    rN = order[ranking[-1]]

    if rN > r1:
        for i in range(N):
            rI = order[ranking[i]]
            koef = (rN - rI) / (rN - r1)
            deviat_dist = np.random.randn(3) * max_sigma * koef
            new_Coord[ranking[i], 0] += deviat_dist[0] / temp_potLat[0]
            new_Coord[ranking[i], 1] += deviat_dist[1] / temp_potLat[1]
            new_Coord[ranking[i], 2] += deviat_dist[2] / temp_potLat[2]
            new_Coord[ranking[i], 0] -= np.floor(new_Coord[ranking[i], 0])
            new_Coord[ranking[i], 1] -= np.floor(new_Coord[ranking[i], 1])
            new_Coord[ranking[i], 2] -= np.floor(new_Coord[ranking[i], 2])

    return new_Coord


def move_all_mol_vectorized(ZMATRIX_centers, LATTICE, order, max_sigma):
    """向量化版本"""
    N = len(ZMATRIX_centers)
    new_Coord = ZMATRIX_centers.copy()
    new_Coord = new_Coord @ np.linalg.inv(LATTICE)

    if len(LATTICE.shape) == 1 or LATTICE.shape[0] == 6:
        lat = LATTICE if LATTICE.shape[0] == 3 else LATTICE[:3].reshape(3, 3)
    else:
        lat = LATTICE
    temp_potLat = np.array([np.linalg.norm(lat[0]),
                             np.linalg.norm(lat[1]),
                             np.linalg.norm(lat[2])])

    ranking = np.argsort(order)
    r1 = order[ranking[0]]
    rN = order[ranking[-1]]

    if rN > r1:
        rI_all = order[ranking]
        koef_all = (rN - rI_all) / (rN - r1)  # N×1

        deviat_dist = np.random.randn(N, 3) * (max_sigma * koef_all[:, np.newaxis])
        deviat_frac = deviat_dist / temp_potLat

        new_Coord[ranking] += deviat_frac
        new_Coord[ranking] -= np.floor(new_Coord[ranking])

    return new_Coord


# 测试
np.random.seed(42)
n_tests = 1000
agree = 0

for t in range(n_tests):
    N = 8
    ZMATRIX_centers = np.random.randn(N, 3) * 10
    LATTICE = np.array([[20, 2, 1], [0, 18, -3], [0, 0, 22]], dtype=float)
    order = np.random.randint(1, 100, N)
    max_sigma = 0.5

    # 用相同随机种子
    np.random.seed(t)
    r1 = move_all_mol_original(ZMATRIX_centers, LATTICE, order, max_sigma)

    np.random.seed(t)
    r2 = move_all_mol_vectorized(ZMATRIX_centers, LATTICE, order, max_sigma)

    if np.allclose(r1, r2):
        agree += 1
    else:
        print(f"✗ 测试 {t}: max diff = {np.max(np.abs(r1 - r2))}")
        break

print(f"一致: {agree}/{n_tests} ({agree/n_tests*100:.1f}%)")
