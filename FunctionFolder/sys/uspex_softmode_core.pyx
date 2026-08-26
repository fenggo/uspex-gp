# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True
# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

"""
uspex_softmode_core.pyx — USPEX calcSoftModes_molecules.m 的忠实 Cython 翻译

逐行对照 calcSoftModes_molecules.m (325行) 实现，不跳过任何计算步骤。
  1. 键搜索 + 插入排序 (行35-89)
  2. 键类型分类 + colors 连通性检查 (行90-222)
  3. nu_factor + 力常数 + 动力学矩阵 D (行223-322)
  4. eig(D) 对角化 (行323-325)
"""

import numpy as np
cimport numpy as cnp
from libc.math cimport sqrt, exp, abs as cabs, log
cimport cython

cnp.import_array()
ctypedef double f64_t


def calc_soft_modes_raw(f64_t[:, ::1] coords, f64_t[:, ::1] lat,
                         cnp.int32_t[:] at_types, f64_t[:] R_val,
                         f64_t[:] N_val_arr, f64_t[:] val_arr,
                         f64_t[:, ::1] goodBonds):
    """
    完整替代 calcSoftModes_molecules.m

    参数 (全部与 Octave 1-based 一致，除 at_types):
      coords:    (N, 3) 分数坐标
      lat:       (3, 3) 晶格矩阵
      at_types:  (N,) 0-based 原子类型索引 (Python侧)
      R_val:     (n_types,) 共价半径
      N_val_arr: (n_types,) 价电子数
      val_arr:   (n_types,) 价数
      goodBonds: (n_types, n_types) 矩阵 (与 ORG_STRUC.goodBonds 一致)

    返回 (与 Octave eig(D) 格式一致):
      freq_mat:   (3N, 3N) 对角矩阵, diag = abs(real(eigenvalues))
      eigvector:  (3N, 3N) 实部特征向量矩阵, 每列一个特征向量
    """
    cdef int N_i = coords.shape[0]
    cdef int n_types = R_val.shape[0]
    cdef f64_t same_bond = 0.05
    cdef f64_t max_bond = 5.0
    cdef int nL = 2
    cdef int nL1 = 2 * nL + 1  # = 5

    # 所有循环变量提前声明 (Cython 要求)
    cdef int i, j, k1, k2, k3, bb, bbb, it_idx, r1, r2, r3, m
    cdef int n1, n2, n3_idx, pos, bb2
    cdef int bond_types = 0, N_bonds = 0, N_bonds1 = 0
    cdef int max_bonds = N_i * N_i * 27
    cdef int n3 = 3 * N_i
    cdef int a, b, a3, b3, it2, bt
    cdef int a_type, b_type, bond_used
    cdef int ci_val, cj_val
    cdef bint connected
    cdef long c1, c2
    cdef int k, k1_idx, k2_idx, k3_idx
    cdef f64_t vx, vy, vz, dist, delta, delta_cutoff
    cdef f64_t R_a, R_b, nu, EN_a, EN_b, CN_a, CN_b, f_ab, X_ab, H
    cdef f64_t dx, dy, dz, dnorm, cos_x, cos_y, cos_z
    cdef f64_t nu_full
    cdef f64_t[:, ::1] sb, bonds_v, bs, D_view

    # ── 行14: small_bond = -0.37 * log(goodBonds)  (矩阵) ──
    small_bond = -0.37 * np.log(np.asarray(goodBonds))
    sb = np.ascontiguousarray(small_bond)

    # ═══════════════════════════════════════════════════════════
    # 行35-89: 键搜索 + 插入排序
    # ═══════════════════════════════════════════════════════════

    # 预分配 (替代 vertcat)
    bonds_buf = np.zeros((max_bonds, 7), dtype=np.float64)
    cdef f64_t[:, ::1] bonds = bonds_buf

    for i in range(N_i):
        for j in range(i, N_i):
            for k1 in range(-1, 2):
                for k2 in range(-1, 2):
                    for k3 in range(-1, 2):
                        if i == j and k1 == 0 and k2 == 0 and k3 == 0:
                            continue
                        vx = coords[i, 0] + k1 - coords[j, 0]
                        vy = coords[i, 1] + k2 - coords[j, 1]
                        vz = coords[i, 2] + k3 - coords[j, 2]
                        # vect * lat  (行60: sqrt(sum((vect*lat).^2)))
                        dist = sqrt((vx*lat[0,0]+vy*lat[1,0]+vz*lat[2,0])**2 +
                                    (vx*lat[0,1]+vy*lat[1,1]+vz*lat[2,1])**2 +
                                    (vx*lat[0,2]+vy*lat[1,2]+vz*lat[2,2])**2)
                        delta = dist - R_val[at_types[i]] - R_val[at_types[j]]
                        if delta < max_bond:
                            # 插入排序: 找到插入位置, 使 bonds[:N_bonds1] 按 delta 升序
                            # (行70-82 的插入排序逻辑)
                            pos = N_bonds1
                            for bb in range(N_bonds1):
                                if bonds[bb, 2] > delta:
                                    pos = bb
                                    break
                            # 后移
                            for bb2 in range(N_bonds1, pos, -1):
                                bonds[bb2, 0] = bonds[bb2-1, 0]
                                bonds[bb2, 1] = bonds[bb2-1, 1]
                                bonds[bb2, 2] = bonds[bb2-1, 2]
                                bonds[bb2, 3] = bonds[bb2-1, 3]
                                bonds[bb2, 4] = bonds[bb2-1, 4]
                                bonds[bb2, 5] = bonds[bb2-1, 5]
                                bonds[bb2, 6] = bonds[bb2-1, 6]
                            # 插入 (1-based, 与 Octave 一致)
                            bonds[pos, 0] = i + 1
                            bonds[pos, 1] = j + 1
                            bonds[pos, 2] = delta
                            bonds[pos, 3] = 0
                            bonds[pos, 4] = k1
                            bonds[pos, 5] = k2
                            bonds[pos, 6] = k3
                            N_bonds1 += 1

    if N_bonds1 == 0:
        return np.zeros((n3, n3)), np.eye(n3)

    # 截取有效部分
    bonds = bonds_buf[:N_bonds1, :].copy()
    bonds_v = bonds

    # ═══════════════════════════════════════════════════════════
    # 行36-46: colors 数组初始化 (5×5×5 × N_i)
    # ═══════════════════════════════════════════════════════════
    # colors[k1,k2,k3,i] = i + nL1^2*N_i*(k1-1) + nL1*N_i*(k2-1) + N_i*(k3-1)
    # (1-based indices in Octave; we use 0-based internally)
    colors = np.zeros((nL1, nL1, nL1, N_i), dtype=np.int64)
    cdef long[:, :, :, ::1] colors_v = colors
    cdef long[:, :, :, ::1] colors1

    for k1 in range(nL1):
        for k2 in range(nL1):
            for k3 in range(nL1):
                for i in range(N_i):
                    colors_v[k1, k2, k3, i] = (i + 1 +
                        nL1*nL1*N_i*k1 + nL1*N_i*k2 + N_i*k3)

    colors1 = colors.copy()

    # ═══════════════════════════════════════════════════════════
    # 行90-222: 键类型分类 + colors 连通性检查
    # ═══════════════════════════════════════════════════════════

    bonds_sorted = np.zeros((N_bonds1, 7), dtype=np.float64)
    bs = bonds_sorted

    for bb in range(N_bonds1):
        if bonds_v[bb, 3] == 0:
            bond_types += 1
            bonds_v[bb, 3] = bond_types
            N_bonds += 1
            # bonds_sorted(N_bonds,:) = bonds(bb,:)  (行98-100)
            bs[N_bonds-1, :] = bonds_v[bb, :]

            for bbb in range(bb + 1, N_bonds1):
                if bonds_v[bbb, 2] > bonds_v[bb, 2] + same_bond:
                    break
                if (bonds_v[bbb, 3] == 0 and
                    at_types[<int>bonds_v[bbb, 0]-1] == at_types[<int>bonds_v[bb, 0]-1] and
                    at_types[<int>bonds_v[bbb, 1]-1] == at_types[<int>bonds_v[bb, 1]-1]):
                    bonds_v[bbb, 3] = bond_types
                    N_bonds += 1
                    bs[N_bonds-1, :] = bonds_v[bbb, :]

    # 行114-222: colors 连通性检查循环
    it_idx = -1  # Octave it=0, 用-1因为下面+1后是0-based

    for bt in range(1, bond_types + 1):
        while it_idx + 1 < N_bonds and <int>bs[it_idx + 1, 3] == bt:
            it_idx += 1

            a_type = at_types[<int>bs[it_idx, 0]-1] + 1  # 1-based type
            b_type = at_types[<int>bs[it_idx, 1]-1] + 1
            bond_used = 0

            # 行123-151: colors union-find
            for m1 in range(nL1):
                for m2 in range(nL1):
                    for m3 in range(nL1):
                        n1 = m1 + <int>bs[it_idx, 4]  # bonds_sorted(it,5)=k1
                        n2 = m2 + <int>bs[it_idx, 5]  # k2
                        n3_idx = m3 + <int>bs[it_idx, 6]  # k3
                        if n1 < 0 or n1 >= nL1 or n2 < 0 or n2 >= nL1 or n3_idx < 0 or n3_idx >= nL1:
                            continue
                        # 行132: 检查是否已连通 + small_bond 矩阵检查
                        # colors(nL+1,nL+1,nL+1,bonds_sorted(it,2)) == colors(nL+1+k1,nL+1+k2,nL+1+k3,bonds_sorted(it,1))
                        # 注意: Octave nL+1=3 (1-based), Python nL=2 (0-based)
                        if (colors_v[nL, nL, nL, <int>bs[it_idx, 1]-1] ==
                            colors_v[nL+<int>bs[it_idx, 4],
                                      nL+<int>bs[it_idx, 5],
                                      nL+<int>bs[it_idx, 6],
                                      <int>bs[it_idx, 0]-1] and
                            bs[it_idx, 2] > sb[a_type-1, b_type-1]):
                            continue

                        bond_used = 1
                        # union: colors1(m1,m2,m3,atom2) -> colors1(n1,n2,n3,atom1)
                        cj_val = colors1[m1, m2, m3, <int>bs[it_idx, 1]-1]
                        ci_val = colors1[n1, n2, n3_idx, <int>bs[it_idx, 0]-1]
                        for r1 in range(nL1):
                            for r2 in range(nL1):
                                for r3 in range(nL1):
                                    for m in range(N_i):
                                        if colors1[r1, r2, r3, m] == cj_val:
                                            colors1[r1, r2, r3, m] = ci_val

            # 行153-156: bond_used=0 -> 标记无效
            if bond_used == 0:
                bs[it_idx, 0] = 0
                bs[it_idx, 1] = 0

            # 行161-187: 额外的 colors1 同步
            for i in range(N_i):
                for j in range(i + 1, N_i):
                    if colors1[nL, nL, nL, i] == colors1[nL, nL, nL, j]:
                        for m1 in range(nL1):
                            for m2 in range(nL1):
                                for m3 in range(nL1):
                                    if colors1[m1, m2, m3, i] != colors1[m1, m2, m3, j]:
                                        cj_val = colors1[m1, m2, m3, i]
                                        ci_val = colors1[m1, m2, m3, j]
                                        for r1 in range(nL1):
                                            for r2 in range(nL1):
                                                for r3 in range(nL1):
                                                    for m in range(N_i):
                                                        if colors1[r1, r2, r3, m] == cj_val:
                                                            colors1[r1, r2, r3, m] = ci_val

            # 行188-222: 连通性检查 + 跳过剩余键
            if it_idx + 1 < N_bonds:
                if bs[it_idx + 1, 2] < np.max(sb):
                    colors1 = colors1  # continue (行190 continue 跳到下一个 bt)
                    if it_idx + 1 < N_bonds:
                        continue  # 对应行190 continue

            colors_v = colors1.copy()
            connected = True
            c1 = colors_v[0, 0, 0, 0]
            c2 = 0
            for k1_idx in range(nL, nL + 3):
                for k2_idx in range(nL, nL + 3):
                    for k3_idx in range(nL, nL + 3):
                        for i in range(N_i):
                            if nL == 0 and colors_v[0, 0, 0, i] != c1:
                                connected = False
                                c2 = 1
                            else:
                                if colors_v[k1_idx, k2_idx, k3_idx, i] != c1 and c2 == 0:
                                    c2 = colors_v[k1_idx, k2_idx, k3_idx, i]
                                elif (colors_v[k1_idx, k2_idx, k3_idx, i] != c1 and
                                      colors_v[k1_idx, k2_idx, k3_idx, i] != c2):
                                    connected = False

            if connected or it_idx >= N_bonds - 1:
                for dumm in range(it_idx + 1, N_bonds):
                    bs[dumm, 0] = 0
                    bs[dumm, 1] = 0
                break  # 行220 break -> 跳出 bt 循环

    # ═══════════════════════════════════════════════════════════
    # 行223-232: nu_factor
    # ═══════════════════════════════════════════════════════════
    nu_factor = np.zeros(N_i, dtype=np.float64)
    cdef f64_t[:] nu_f = nu_factor

    for k in range(N_i):
        nu_full = 0.0
        for m in range(N_bonds):
            if <int>bs[m, 0] == k + 1 or <int>bs[m, 1] == k + 1:
                nu_full += exp(-1.0 * bs[m, 2] / 0.37)
        # 行231: nu_factor(k) = val(at_types(k)) / nu_full  (无除零保护!)
        nu_f[k] = val_arr[at_types[k]] / nu_full

    # ═══════════════════════════════════════════════════════════
    # 行233-322: 动力学矩阵 D
    # ═══════════════════════════════════════════════════════════
    D = np.zeros((n3, n3), dtype=np.float64)
    D_view = D

    it2 = -1
    for bt in range(1, bond_types + 1):
        while it2 + 1 < N_bonds and <int>bs[it2 + 1, 3] == bt:
            it2 += 1
            if <int>bs[it2, 0] == 0:
                if it2 == N_bonds - 1:
                    break
                continue

            # 行243-244: a = at_types(bonds_sorted(it,1)); b = at_types(bonds_sorted(it,2))
            # 注意: Octave at_types 是 1-based
            a = at_types[<int>bs[it2, 0]-1] + 1  # 1-based type
            b = at_types[<int>bs[it2, 1]-1] + 1

            delta = bs[it2, 2]
            R_a = R_val[a-1] + delta / 2
            R_b = R_val[b-1] + delta / 2

            # 行248-255
            if R_a < 0.05:
                R_b = R_b - 0.05 + R_a
                R_a = 0.05
            if R_b < 0.05:
                R_a = R_a - 0.05 + R_b
                R_b = 0.05

            nu = exp(-delta / 0.37)
            EN_a = 0.481 * N_val_arr[a-1] / R_a
            EN_b = 0.481 * N_val_arr[b-1] / R_b
            # 行259-260: 无除零保护 (与原版一致)
            CN_a = val_arr[a-1] / (nu * nu_f[<int>bs[it2, 0]-1])
            CN_b = val_arr[b-1] / (nu * nu_f[<int>bs[it2, 1]-1])

            f_ab = 0.25 * abs(EN_a - EN_b) / sqrt(EN_a * EN_b)
            X_ab = sqrt((EN_a * EN_b) / (CN_a * CN_b))

            # 行263-264: 交换 a 和 b !!!
            # b = bonds_sorted(it,1); a = bonds_sorted(it,2);
            b = <int>bs[it2, 0]  # 1-based atom index
            a = <int>bs[it2, 1]  # 1-based atom index

            k1 = <int>bs[it2, 4]
            k2 = <int>bs[it2, 5]
            k3 = <int>bs[it2, 6]

            # 行268-270: vect = coords(b,:) + [k1 k2 k3] - coords(a,:)
            vx = coords[b-1, 0] + k1 - coords[a-1, 0]
            vy = coords[b-1, 1] + k2 - coords[a-1, 1]
            vz = coords[b-1, 2] + k3 - coords[a-1, 2]
            # vect * lat
            dx = vx*lat[0,0] + vy*lat[1,0] + vz*lat[2,0]
            dy = vx*lat[0,1] + vy*lat[1,1] + vz*lat[2,1]
            dz = vx*lat[0,2] + vy*lat[1,2] + vz*lat[2,2]
            dnorm = sqrt(dx*dx + dy*dy + dz*dz)

            cos_x = dx / dnorm
            cos_y = dy / dnorm
            cos_z = dz / dnorm
            cos_x = round(cos_x * 1000000) / 1000000
            cos_y = round(cos_y * 1000000) / 1000000
            cos_z = round(cos_z * 1000000) / 1000000

            H = X_ab * exp(-2.7 * f_ab)

            # 行279-317: if a == b, skip; else 填充 D
            if a != b:
                a3 = (a - 1) * 3  # 0-based 矩阵索引
                b3 = (b - 1) * 3

                D_view[a3,   a3  ] += H*cos_x*cos_x
                D_view[a3,   a3+1] += H*cos_x*cos_y
                D_view[a3,   a3+2] += H*cos_x*cos_z
                D_view[a3,   b3  ] -= H*cos_x*cos_x
                D_view[a3,   b3+1] -= H*cos_x*cos_y
                D_view[a3,   b3+2] -= H*cos_x*cos_z

                D_view[a3+1, a3  ] += H*cos_y*cos_x
                D_view[a3+1, a3+1] += H*cos_y*cos_y
                D_view[a3+1, a3+2] += H*cos_y*cos_z
                D_view[a3+1, b3  ] -= H*cos_y*cos_x
                D_view[a3+1, b3+1] -= H*cos_y*cos_y
                D_view[a3+1, b3+2] -= H*cos_y*cos_z

                D_view[a3+2, a3  ] += H*cos_z*cos_x
                D_view[a3+2, a3+1] += H*cos_z*cos_y
                D_view[a3+2, a3+2] += H*cos_z*cos_z
                D_view[a3+2, b3  ] -= H*cos_z*cos_x
                D_view[a3+2, b3+1] -= H*cos_z*cos_y
                D_view[a3+2, b3+2] -= H*cos_z*cos_z

                D_view[b3,   b3  ] += H*cos_x*cos_x
                D_view[b3,   b3+1] += H*cos_x*cos_y
                D_view[b3,   b3+2] += H*cos_x*cos_z
                D_view[b3,   a3  ] -= H*cos_x*cos_x
                D_view[b3,   a3+1] -= H*cos_x*cos_y
                D_view[b3,   a3+2] -= H*cos_x*cos_z

                D_view[b3+1, b3  ] += H*cos_y*cos_x
                D_view[b3+1, b3+1] += H*cos_y*cos_y
                D_view[b3+1, b3+2] += H*cos_y*cos_z
                D_view[b3+1, a3  ] -= H*cos_y*cos_x
                D_view[b3+1, a3+1] -= H*cos_y*cos_y
                D_view[b3+1, a3+2] -= H*cos_y*cos_z

                D_view[b3+2, b3  ] += H*cos_z*cos_x
                D_view[b3+2, b3+1] += H*cos_z*cos_y
                D_view[b3+2, b3+2] += H*cos_z*cos_z
                D_view[b3+2, a3  ] -= H*cos_z*cos_x
                D_view[b3+2, a3+1] -= H*cos_z*cos_y
                D_view[b3+2, a3+2] -= H*cos_z*cos_z

            if it2 == N_bonds - 1:
                break

    # ═══════════════════════════════════════════════════════════
    # 行323-325: [eigvector, freq] = eig(D); eigvector=real(eigvector); freq=abs(freq)
    # ═══════════════════════════════════════════════════════════
    # numpy: w, V = eig(D)  (w=1D eigenvalues, V=eigenvectors column-wise)
    # Octave: [V, D] = eig(D)  (V=eigenvectors, D=diagonal matrix of eigenvalues)
    freq_raw, eigvector = np.linalg.eig(D)
    eigvector = np.real(eigvector)
    freq_raw = np.abs(np.real(freq_raw))

    # 返回对角矩阵 (与 Octave eig(D) 输出一致)
    freq_mat = np.diag(freq_raw)

    return freq_mat, eigvector
