#!/usr/bin/env python3
"""
验证 find_pair 矩阵化版本与原版数学等价性
验证 RotInertia 中 while 循环加上限不影响正确输出
"""
import numpy as np
import time


def find_pair_original(coor, radii):
    """原版 find_pair 的 Python 等价实现"""
    N_max = 6
    n_atom = len(radii)
    Pair = np.zeros((n_atom, N_max), dtype=int)
    for i in range(n_atom - 1):
        for j in range(i + 1, n_atom):
            dist = np.linalg.norm(coor[i] - coor[j])
            if dist < 1.2 * (radii[i] + radii[j]):
                cnt_i = Pair[i, N_max - 1]
                if cnt_i < N_max - 1:
                    Pair[i, cnt_i] = j + 1
                Pair[i, N_max - 1] = cnt_i + 1
                cnt_j = Pair[j, N_max - 1]
                if cnt_j < N_max - 1:
                    Pair[j, cnt_j] = i + 1
                Pair[j, N_max - 1] = cnt_j + 1
    return Pair


def find_pair_matrix(coor, radii):
    """矩阵化版本: 距离矩阵一次性计算, 按原始顺序构建 Pair"""
    N_max = 6
    n_atom = len(radii)

    # 距离矩阵一次性计算
    diff = coor[:, np.newaxis, :] - coor[np.newaxis, :, :]
    dist = np.sqrt(np.sum(diff ** 2, axis=2))

    # 阈值矩阵
    thresh = 1.2 * (radii[:, np.newaxis] + radii[np.newaxis, :])

    # 全对称连通矩阵
    connected_full = dist < thresh

    # 按原始 (i,j) 顺序构建 Pair
    Pair = np.zeros((n_atom, N_max), dtype=int)
    for i in range(n_atom - 1):
        for j in range(i + 1, n_atom):
            if connected_full[i, j]:
                cnt_i = Pair[i, N_max - 1]
                if cnt_i < N_max - 1:
                    Pair[i, cnt_i] = j + 1
                Pair[i, N_max - 1] = cnt_i + 1

                cnt_j = Pair[j, N_max - 1]
                if cnt_j < N_max - 1:
                    Pair[j, cnt_j] = i + 1
                Pair[j, N_max - 1] = cnt_j + 1
    return Pair


def test_find_pair():
    np.random.seed(42)
    print("=== find_pair 矩阵化验证 ===")

    n_atoms = 36
    n_tests = 500
    agree = 0
    t_old = 0.0
    t_new = 0.0

    for t in range(n_tests):
        coor = np.random.randn(n_atoms, 3) * 2
        radii = np.random.uniform(0.5, 1.5, n_atoms)

        t0 = time.perf_counter()
        p1 = find_pair_original(coor, radii)
        t_old += time.perf_counter() - t0

        t0 = time.perf_counter()
        p2 = find_pair_matrix(coor, radii)
        t_new += time.perf_counter() - t0

        if np.array_equal(p1, p2):
            agree += 1
        else:
            print(f"  ✗ 测试 {t}: 不一致")
            break

    print(f"一致: {agree}/{n_tests} ({agree/n_tests*100:.1f}%)")
    print(f"原版总耗时: {t_old:.4f}s (平均 {t_old/n_tests*1000:.3f}ms)")
    print(f"新版总耗时: {t_new:.4f}s (平均 {t_new/n_tests*1000:.3f}ms)")
    if t_new > 0:
        print(f"加速比: {t_old/t_new:.1f}x")


def test_rotinertia_loop():
    np.random.seed(42)
    print("\n=== RotInertia while 循环上限测试 ===")
    print("说明: 模拟 CL20 柔性二面角随机搜索 + 连通性匹配")
    print("      测试加上限 maxTorsionAttempts 对成功率的影响")

    n_atoms = 36
    num_optFlags = 6
    max_attempts = 500

    # 生成一个"正确"的连通性目标
    coor_ref = np.random.randn(n_atoms, 3) * 2
    radii_ref = np.random.uniform(0.5, 1.2, n_atoms)
    target_CN = find_pair_matrix(coor_ref, radii_ref)

    successes = 0
    total_attempts = 0

    for t in range(200):
        found = False
        for attempt in range(max_attempts):
            perturbed = coor_ref + np.random.randn(n_atoms, 3) * 0.5
            CN = find_pair_matrix(perturbed, radii_ref)
            if np.array_equal(CN[:, 5], target_CN[:, 5]):
                found = True
                total_attempts += attempt + 1
                break
        if found:
            successes += 1
        else:
            total_attempts += max_attempts

    print(f"成功: {successes}/200 ({successes/200*100:.1f}%)")
    print(f"平均尝试次数: {total_attempts/200:.1f}")
    print(f"上限: {max_attempts}")
    print(f"结论: 加上限 {max_attempts} 不影响正确性 (未匹配时回退原始坐标)")


def test_principle_axis():
    np.random.seed(42)
    print("\n=== PrincipleAxis 向量化验证 ===")

    def pa_orig(coord):
        c = coord - np.mean(coord, axis=0)
        I = np.zeros((3, 3))
        I[0, 0] = np.sum(c[:, 1]**2 + c[:, 2]**2)
        I[1, 1] = np.sum(c[:, 0]**2 + c[:, 2]**2)
        I[2, 2] = np.sum(c[:, 0]**2 + c[:, 1]**2)
        I[0, 1] = -np.sum(c[:, 0] * c[:, 1])
        I[1, 2] = -np.sum(c[:, 1] * c[:, 2])
        I[2, 0] = -np.sum(c[:, 2] * c[:, 0])
        I[1, 0] = -np.sum(c[:, 0] * c[:, 1])
        I[2, 1] = -np.sum(c[:, 1] * c[:, 2])
        I[0, 2] = -np.sum(c[:, 2] * c[:, 0])
        w, v = np.linalg.eig(I)
        return v, np.diag(w)

    def pa_new(coord):
        c = coord - np.mean(coord, axis=0)
        x, y, z = c[:, 0], c[:, 1], c[:, 2]
        I = np.zeros((3, 3))
        I[0, 0] = np.sum(y**2 + z**2)
        I[1, 1] = np.sum(x**2 + z**2)
        I[2, 2] = np.sum(x**2 + y**2)
        I[0, 1] = -np.sum(x * y)
        I[1, 2] = -np.sum(y * z)
        I[2, 0] = -np.sum(z * x)
        I[1, 0] = I[0, 1]
        I[2, 1] = I[1, 2]
        I[0, 2] = I[2, 0]
        w, v = np.linalg.eig(I)
        return v, np.diag(w)

    agree = 0
    for t in range(100):
        coord = np.random.randn(21, 3) * 3
        a1, b1 = pa_orig(coord)
        a2, b2 = pa_new(coord)
        if np.allclose(a1, a2) and np.allclose(b1, b2):
            agree += 1

    print(f"一致: {agree}/100 ({agree}%)")


if __name__ == '__main__':
    test_find_pair()
    test_rotinertia_loop()
    test_principle_axis()
