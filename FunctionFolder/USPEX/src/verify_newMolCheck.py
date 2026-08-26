#!/usr/bin/env python3
"""
验证 newMolCheck 矩阵化版本与原版数学等价性

算法对照:
  Part 1 (分子间): 原版四重循环 CalcDist → 新版距离矩阵 + 阈值矩阵一次性比较
  Part 2 (分子内): 原版 bsxfun + pdist2 → 新版广播 + min + 阈值向量

两者使用相同的 inv(LATTICE)/round/norm 路径，数学等价。
"""

import numpy as np
import time


def calc_dist_original(lattice, coord1, coord2):
    """原版 CalcDist 的 Python 等价实现"""
    inv = np.linalg.inv(lattice)
    c1 = coord1 @ inv
    c2 = coord2 @ inv
    diff = c1 - c2
    diff = diff - np.round(diff)
    return np.linalg.norm(diff @ lattice)


def newMolCheck_original(MOLECULES, LATTICE, MtypeLIST, STDMOL, miniMatrix):
    """原版 newMolCheck 的 Python 等价实现"""
    # Part 1: 分子间
    for i in range(len(MtypeLIST) - 1):
        n_i = len(STDMOL[MtypeLIST[i]]['types'])
        for j in range(n_i):
            for k in range(i + 1, len(MtypeLIST)):
                n_k = len(STDMOL[MtypeLIST[k]]['types'])
                for l in range(n_k):
                    dist = calc_dist_original(
                        LATTICE,
                        MOLECULES[i]['coords'][j],
                        MOLECULES[k]['coords'][l],
                    )
                    t1 = STDMOL[MtypeLIST[i]]['types'][j]
                    t2 = STDMOL[MtypeLIST[k]]['types'][l]
                    if dist <= miniMatrix[t1, t2]:
                        return 0
    # Part 2: 分子内
    X = np.arange(-2, 3)
    Y = np.arange(-2, 3)
    Z = np.arange(-2, 3)
    X2, Y2, Z2 = np.meshgrid(X, Y, Z, indexing='ij')
    Matrix = np.column_stack([X2.ravel(), Y2.ravel(), Z2.ravel()])
    Matrix = Matrix[~np.all(Matrix == 0, axis=1)]  # 124 × 3

    for ind in range(len(MtypeLIST)):
        n_atom = len(STDMOL[MtypeLIST[ind]]['types'])
        for m in range(n_atom - 1):
            for n in range(m + 1, n_atom):
                coor1 = MOLECULES[ind]['coords'][m]
                coor2 = MOLECULES[ind]['coords'][n]
                tmp = coor2 + Matrix @ LATTICE  # 124 × 3
                # pdist2: 124 个点到 coor1 的距离
                dists = np.sqrt(np.sum((tmp - coor1) ** 2, axis=1))
                min_dist = np.min(dists)
                t1 = STDMOL[MtypeLIST[ind]]['types'][m]
                t2 = STDMOL[MtypeLIST[ind]]['types'][n]
                if min_dist <= miniMatrix[t1, t2]:
                    return 0
    return 1


def newMolCheck_matrix(MOLECULES, LATTICE, MtypeLIST, STDMOL, miniMatrix):
    """矩阵化优化版本 — 与原版数学等价"""
    nMols = len(MtypeLIST)
    invLat = np.linalg.inv(LATTICE)

    # 预收集
    molCoords = [MOLECULES[i]['coords'] for i in range(nMols)]
    molTypes = [np.array(STDMOL[MtypeLIST[i]]['types']) for i in range(nMols)]

    # ========= Part 1: 分子间 (矩阵化) =========
    for i in range(nMols - 1):
        coords_i = molCoords[i]
        types_i = molTypes[i]
        n_i = len(coords_i)
        if n_i == 0:
            continue

        frac_i = coords_i @ invLat  # n_i × 3

        for k in range(i + 1, nMols):
            coords_k = molCoords[k]
            types_k = molTypes[k]
            n_k = len(coords_k)
            if n_k == 0:
                continue

            frac_k = coords_k @ invLat  # n_k × 3

            # 所有原子对的分数坐标差: n_i × n_k × 3
            diff_frac = frac_i[:, np.newaxis, :] - frac_k[np.newaxis, :, :]

            # 最小镜像约定
            diff_frac = diff_frac - np.round(diff_frac)

            # 转为笛卡尔距离: n_i × n_k
            diff_cart = diff_frac @ LATTICE  # n_i × n_k × 3
            dists = np.sqrt(np.sum(diff_cart ** 2, axis=2))  # n_i × n_k

            # 阈值矩阵
            thresh = miniMatrix[types_i[:, np.newaxis], types_k[np.newaxis, :]]

            if np.any(dists <= thresh):
                return 0

    # ========= Part 2: 分子内 (矩阵化) =========
    # 预计算 124 个周期偏移
    X = np.arange(-2, 3)
    Y = np.arange(-2, 3)
    Z = np.arange(-2, 3)
    X2, Y2, Z2 = np.meshgrid(X, Y, Z, indexing='ij')
    offsetMat = np.column_stack([X2.ravel(), Y2.ravel(), Z2.ravel()])
    offsetMat = offsetMat[~np.all(offsetMat == 0, axis=1)]  # 124 × 3
    offsetCart = offsetMat @ LATTICE  # 124 × 3

    for ind in range(nMols):
        coords = molCoords[ind]
        types = molTypes[ind]
        nAtoms = len(coords)

        if nAtoms < 2:
            continue

        for n in range(1, nAtoms):
            # 原子 n 的所有 124 个周期镜像
            images_n = coords[n] + offsetCart  # 124 × 3

            # 原子 0..n-1
            atoms_m = coords[:n]  # n × 3

            # 广播: (n) × 124 × 3
            diff = atoms_m[:, np.newaxis, :] - images_n[np.newaxis, :, :]
            dists = np.sqrt(np.sum(diff ** 2, axis=2))  # n × 124

            # 每个 m 取最小距离
            minDists = np.min(dists, axis=1)  # n × 1

            # 阈值
            thresh = miniMatrix[types[:n], types[n]]

            if np.any(minDists <= thresh):
                return 0

    return 1


def run_test():
    """运行正确性验证和性能测试"""
    np.random.seed(42)

    # 模拟 TNT/CL20 体系
    STDMOL = {
        1: {'types': np.array([0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0])},   # TNT: 21 atoms
        2: {'types': np.array([0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
                               0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3])},  # CL20: 36 atoms
    }

    miniMatrix = np.array([
        [1.74, 1.69, 1.71, 1.29],
        [1.69, 1.64, 1.66, 1.24],
        [1.71, 1.66, 1.69, 1.26],
        [1.29, 1.24, 1.26, 1.02],
    ])

    LATTICE = np.array([[25, 0, 0], [0, 25, 0], [0, 0, 25]], dtype=float)
    MtypeLIST = [1, 1, 1, 1, 2, 2, 2, 2]

    print("=== newMolCheck 矩阵化版本验证 ===")
    print(f"体系: TNT4·CL20₄, 8 分子, 228 原子")
    print(f"miniMatrix:\n{miniMatrix}")

    # 正确性验证
    nTests = 500
    agree = 0
    disagree = 0
    t_old = 0.0
    t_new = 0.0

    print(f"\n运行 {nTests} 组随机分子配置...")

    for t in range(nTests):
        MOLECULES = []
        for mt in MtypeLIST:
            nAtoms = len(STDMOL[mt]['types'])
            MOLECULES.append({'coords': np.random.rand(nAtoms, 3) * 20})

        t0 = time.perf_counter()
        r1 = newMolCheck_original(MOLECULES, LATTICE, MtypeLIST, STDMOL, miniMatrix)
        t_old += time.perf_counter() - t0

        t0 = time.perf_counter()
        r2 = newMolCheck_matrix(MOLECULES, LATTICE, MtypeLIST, STDMOL, miniMatrix)
        t_new += time.perf_counter() - t0

        if r1 == r2:
            agree += 1
        else:
            disagree += 1
            print(f"  ✗ 测试 {t}: 原版={r1}, 新版={r2}")

    print(f"\n=== 正确性 ===")
    print(f"一致: {agree}/{nTests} ({agree/nTests*100:.1f}%)")
    if disagree > 0:
        print(f"不一致: {disagree} (需检查浮点精度)")

    print(f"\n=== 性能 ===")
    print(f"原版总耗时: {t_old:.4f} s (平均 {t_old/nTests*1000:.3f} ms)")
    print(f"新版总耗时: {t_new:.4f} s (平均 {t_new/nTests*1000:.3f} ms)")
    if t_new > 0:
        print(f"加速比: {t_old/t_new:.1f}x")

    # 大规模测试
    print(f"\n=== 大规模压力测试 (2000 次) ===")
    t0 = time.perf_counter()
    for t in range(2000):
        MOLECULES = []
        for mt in MtypeLIST:
            nAtoms = len(STDMOL[mt]['types'])
            MOLECULES.append({'coords': np.random.rand(nAtoms, 3) * 20})
        newMolCheck_matrix(MOLECULES, LATTICE, MtypeLIST, STDMOL, miniMatrix)
    t_large = time.perf_counter() - t0
    print(f"2000 次新版耗时: {t_large:.3f} s (平均 {t_large/2000*1000:.3f} ms)")

    # 验证与原版输出完全一致的特殊测试
    print(f"\n=== 逐位验证 (固定种子) ===")
    np.random.seed(12345)
    MOLECULES = []
    for mt in MtypeLIST:
        nAtoms = len(STDMOL[mt]['types'])
        MOLECULES.append({'coords': np.random.rand(nAtoms, 3) * 20})

    r1 = newMolCheck_original(MOLECULES, LATTICE, MtypeLIST, STDMOL, miniMatrix)
    r2 = newMolCheck_matrix(MOLECULES, LATTICE, MtypeLIST, STDMOL, miniMatrix)
    print(f"固定种子: 原版={r1}, 新版={r2}, 一致={r1==r2}")

    # 边界测试: 距离恰好等于阈值
    print(f"\n=== 边界测试: 距离恰好等于阈值 ===")
    MOLECULES = []
    for mt in MtypeLIST:
        nAtoms = len(STDMOL[mt]['types'])
        MOLECULES.append({'coords': np.random.rand(nAtoms, 3) * 20})

    # 人为制造一个距离恰好等于阈值的原子对
    MOLECULES[0]['coords'][0] = np.array([0.0, 0.0, 0.0])
    MOLECULES[0]['coords'][1] = np.array([miniMatrix[0, 0], 0.0, 0.0])  # C-C = 1.74
    r1 = newMolCheck_original(MOLECULES, LATTICE, MtypeLIST, STDMOL, miniMatrix)
    r2 = newMolCheck_matrix(MOLECULES, LATTICE, MtypeLIST, STDMOL, miniMatrix)
    print(f"距离=阈值: 原版={r1}, 新版={r2}, 一致={r1==r2} (期望=0)")

    # 边界测试: 距离略大于阈值
    MOLECULES[0]['coords'][1] = np.array([miniMatrix[0, 0] + 0.01, 0.0, 0.0])
    r1 = newMolCheck_original(MOLECULES, LATTICE, MtypeLIST, STDMOL, miniMatrix)
    r2 = newMolCheck_matrix(MOLECULES, LATTICE, MtypeLIST, STDMOL, miniMatrix)
    print(f"距离>阈值: 原版={r1}, 新版={r2}, 一致={r1==r2} (期望=1)")


if __name__ == '__main__':
    run_test()
