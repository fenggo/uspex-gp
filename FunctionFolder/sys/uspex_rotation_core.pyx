# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True
# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

"""
uspex_rotation_core.pyx — USPEX Rotation_310 Cython 加速核心

加速 SmartRotInertia 中的热路径:
  1. zmatrix_to_coords: Z矩阵 -> 笛卡尔坐标 (替代 NEW_ZMATRIXCOORD)
  2. coords_to_zmatrix: 笛卡尔 -> Z矩阵 (替代 NEW_coord2Zmatrix)
  3. find_pair_fast: 找键连关系 (替代 find_pair)
  4. smart_rot_inertia: 完整二面角迭代循环内联
"""

import numpy as np
cimport numpy as cnp
from libc.math cimport sqrt, sin, cos, acos, atan2, abs as cabs, M_PI
cimport cython

cnp.import_array()

ctypedef double f64_t


# ═══════════════════════════════════════════════════════════════
# 1. 内坐标 -> 笛卡尔 (conStructCOO 内联)
# ═══════════════════════════════════════════════════════════════

cdef inline void _conStructCOO(const f64_t* RC, f64_t r, f64_t theta, f64_t phi,
                                f64_t* out) noexcept nogil:
    cdef f64_t xi = RC[0], yi = RC[1], zi = RC[2]
    cdef f64_t xj = RC[3], yj = RC[4], zj = RC[5]
    cdef f64_t xk = RC[6], yk = RC[7], zk = RC[8]
    cdef f64_t xji = xi - xj, yji = yi - yj, zji = zi - zj
    cdef f64_t rji = sqrt(xji*xji + yji*yji + zji*zji)
    xji /= rji; yji /= rji; zji /= rji
    cdef f64_t xki = xk - xj, yki = yk - yj, zki = zk - zj
    cdef f64_t xn = yji*zki - yki*zji
    cdef f64_t yn = zji*xki - zki*xji
    cdef f64_t zn = xji*yki - xki*yji
    cdef f64_t rn = sqrt(xn*xn + yn*yn + zn*zn)
    xn /= rn; yn /= rn; zn /= rn
    cdef f64_t xp = yn*zji - yji*zn
    cdef f64_t yp = zn*xji - zji*xn
    cdef f64_t zp = xn*yji - xji*yn
    rji = r * cos(theta)
    rn  = r * sin(theta) * sin(phi)
    cdef f64_t rp = r * sin(theta) * cos(phi)
    out[0] = xi - rji*xji + rn*xn + rp*xp
    out[1] = yi - rji*yji + rn*yn + rp*yp
    out[2] = zi - rji*zji + rn*zn + rp*zp


# ═══════════════════════════════════════════════════════════════
# 2. Z矩阵 -> 笛卡尔坐标 (替代 NEW_ZMATRIXCOORD)
# ═══════════════════════════════════════════════════════════════

def zmatrix_to_coords(f64_t[:, ::1] ZMATRIX, cnp.int32_t[:, ::1] fmt):
    """Z矩阵 -> 笛卡尔坐标 (替代 NEW_ZMATRIXCOORD.m)"""
    cdef int N = ZMATRIX.shape[0]
    cdef cnp.ndarray[f64_t, ndim=2] Coords_arr = np.zeros((N, 3), dtype=np.float64)
    cdef f64_t[:, ::1] C = Coords_arr
    cdef int ind, form0, form1, form2, inder
    cdef f64_t r, theta, phi, torsang, c_val, s_val, omc
    cdef f64_t rotVec0, rotVec1, rotVec2
    cdef f64_t RelatedC[9]
    cdef f64_t tx, ty, tz

    if N == 1:
        C[0, 0] = ZMATRIX[0, 0]; C[0, 1] = ZMATRIX[0, 1]; C[0, 2] = ZMATRIX[0, 2]
        return Coords_arr

    tx = ZMATRIX[0, 0]; ty = ZMATRIX[0, 1]; tz = ZMATRIX[0, 2]
    C[0, 0] = 0.0; C[0, 1] = 0.0; C[0, 2] = 0.0
    C[1, 0] = ZMATRIX[1, 0]; C[1, 1] = 0.0; C[1, 2] = 0.0

    if N == 2:
        # 旋转 (与 NEW_ZMATRIXCOORD.m 的 N==2 分支一致)
        _apply_rot2(C, ZMATRIX, 0)  # rotVec=[0,1,0], torsang=Z(2,2)-pi/2
        _apply_rot2(C, ZMATRIX, 1)  # rotVec=[0,0,1], torsang=Z(2,3)
    else:
        for ind in range(2, N):
            r = ZMATRIX[ind, 0]; theta = ZMATRIX[ind, 1]; phi = ZMATRIX[ind, 2]
            if ind == 2:
                RelatedC[0] = ZMATRIX[fmt[ind, 0]-1, 0]
                RelatedC[1] = 0.0; RelatedC[2] = 0.0
                RelatedC[3] = ZMATRIX[fmt[ind, 1]-1, 0]
                RelatedC[4] = 0.0; RelatedC[5] = 0.0
                RelatedC[6] = 0.0; RelatedC[7] = 0.0; RelatedC[8] = 1.0
            else:
                form0 = fmt[ind, 0]-1; form1 = fmt[ind, 1]-1; form2 = fmt[ind, 2]-1
                RelatedC[0] = C[form0, 0]; RelatedC[1] = C[form0, 1]; RelatedC[2] = C[form0, 2]
                RelatedC[3] = C[form1, 0]; RelatedC[4] = C[form1, 1]; RelatedC[5] = C[form1, 2]
                RelatedC[6] = C[form2, 0]; RelatedC[7] = C[form2, 1]; RelatedC[8] = C[form2, 2]
            _conStructCOO(RelatedC, r, theta, phi, &C[ind, 0])

            if ind == 3:
                # NEW_ZMATRIXCOORD.m 对 ind==3 做了特殊旋转
                # outer=1: rotVec=[0,1,0], torsang=Z(2,2)-pi/2
                _apply_rot2_all(C, ZMATRIX, 0, ind)
                # outer=2: rotVec=[0,0,1], torsang=Z(2,3)
                _apply_rot2_all(C, ZMATRIX, 1, ind)

    for ind in range(N):
        C[ind, 0] += tx; C[ind, 1] += ty; C[ind, 2] += tz
    return Coords_arr


cdef inline void _apply_rot2(f64_t[:, ::1] C, f64_t[:, ::1] ZMATRIX, int mode) noexcept nogil:
    """对前2个原子应用旋转"""
    cdef f64_t torsang, rv0, rv1, rv2, c_val, s_val, omc
    cdef int inder
    cdef f64_t x, y, z

    if mode == 0:
        rv0 = 0.0; rv1 = 1.0; rv2 = 0.0
        torsang = ZMATRIX[1, 1] - M_PI / 2
    else:
        rv0 = 0.0; rv1 = 0.0; rv2 = 1.0
        torsang = ZMATRIX[1, 2]

    c_val = cos(torsang); s_val = sin(torsang); omc = 1.0 - c_val
    for inder in range(2):
        x = C[inder, 0]; y = C[inder, 1]; z = C[inder, 2]
        C[inder, 0] = (c_val + rv0*rv0*omc)*x + (rv0*rv1*omc - rv2*s_val)*y + (rv0*rv2*omc + rv1*s_val)*z
        C[inder, 1] = (rv1*rv0*omc + rv2*s_val)*x + (c_val + rv1*rv1*omc)*y + (rv1*rv2*omc - rv0*s_val)*z
        C[inder, 2] = (rv2*rv0*omc - rv1*s_val)*x + (rv2*rv1*omc + rv0*s_val)*y + (c_val + rv2*rv2*omc)*z


cdef inline void _apply_rot2_all(f64_t[:, ::1] C, f64_t[:, ::1] ZMATRIX, int mode, int max_ind) noexcept nogil:
    """对前 max_ind 个原子应用旋转 (ind==3 的情况)"""
    cdef f64_t torsang, rv0, rv1, rv2, c_val, s_val, omc
    cdef int inder
    cdef f64_t x, y, z

    if mode == 0:
        rv0 = 0.0; rv1 = 1.0; rv2 = 0.0
        torsang = ZMATRIX[1, 1] - M_PI / 2
    else:
        rv0 = 0.0; rv1 = 0.0; rv2 = 1.0
        torsang = ZMATRIX[1, 2]

    c_val = cos(torsang); s_val = sin(torsang); omc = 1.0 - c_val
    for inder in range(1, max_ind + 1):
        x = C[inder, 0]; y = C[inder, 1]; z = C[inder, 2]
        C[inder, 0] = (c_val + rv0*rv0*omc)*x + (rv0*rv1*omc - rv2*s_val)*y + (rv0*rv2*omc + rv1*s_val)*z
        C[inder, 1] = (rv1*rv0*omc + rv2*s_val)*x + (c_val + rv1*rv1*omc)*y + (rv1*rv2*omc - rv0*s_val)*z
        C[inder, 2] = (rv2*rv0*omc - rv1*s_val)*x + (rv2*rv1*omc + rv0*s_val)*y + (c_val + rv2*rv2*omc)*z


# ═══════════════════════════════════════════════════════════════
# 3. 笛卡尔 -> Z矩阵 (替代 NEW_coord2Zmatrix)
# ═══════════════════════════════════════════════════════════════

def coords_to_zmatrix(f64_t[:, ::1] coords, cnp.int32_t[:, ::1] fmt):
    """笛卡尔坐标 -> Z矩阵 (替代 NEW_coord2Zmatrix.m)"""
    cdef int N = coords.shape[0]
    cdef cnp.ndarray[f64_t, ndim=2] Zmatrix = np.zeros((N, 3), dtype=np.float64)
    cdef f64_t[:, ::1] Z = Zmatrix
    cdef int ind, i, j, k
    cdef f64_t dx, dy, dz, d2, d, nl
    cdef f64_t v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z
    cdef f64_t n1x, n1y, n1z, n2x, n2y, n2z
    cdef f64_t dot_val, torsang, tmp

    if N == 1:
        Z[0, 0] = coords[0, 0]; Z[0, 1] = coords[0, 1]; Z[0, 2] = coords[0, 2]
        return Zmatrix

    # 拷贝并平移使第一个原子在原点
    cdef cnp.ndarray[f64_t, ndim=2] c_arr = np.array(coords, dtype=np.float64, copy=True)
    cdef f64_t[:, ::1] c = c_arr
    for i in range(N):
        c[i, 0] -= coords[0, 0]; c[i, 1] -= coords[0, 1]; c[i, 2] -= coords[0, 2]

    Z[0, 0] = coords[0, 0]; Z[0, 1] = coords[0, 1]; Z[0, 2] = coords[0, 2]

    # 第二个原子
    dx = c[1, 0]; dy = c[1, 1]; dz = c[1, 2]
    Z[1, 0] = sqrt(dx*dx + dy*dy + dz*dz)
    if cabs(c[1, 2]) < 1e-15:
        Z[1, 1] = M_PI / 2
    else:
        Z[1, 1] = acos(c[1, 2] / Z[1, 0])
    if cabs(c[1, 1]) < 1e-15:
        Z[1, 2] = 0.0
    else:
        Z[1, 2] = atan2(c[1, 1], c[1, 0])

    # 第3~N个原子
    for ind in range(2, N):
        i = fmt[ind, 0] - 1
        j = fmt[ind, 1] - 1

        # 键长
        dx = c[ind, 0] - c[i, 0]; dy = c[ind, 1] - c[i, 1]; dz = c[ind, 2] - c[i, 2]
        Z[ind, 0] = sqrt(dx*dx + dy*dy + dz*dz)

        # 键角: vectorsCrutch1 = coords(ind) - coords(i) = v2
        #        vectorsCrutch2 = coords(i) - coords(j) = -v1 (注意 v1 = coords(j)-coords(i))
        v1x = c[i, 0] - c[j, 0]; v1y = c[i, 1] - c[j, 1]; v1z = c[i, 2] - c[j, 2]
        v2x = c[ind, 0] - c[i, 0]; v2y = c[ind, 1] - c[i, 1]; v2z = c[ind, 2] - c[i, 2]

        nl = sqrt(v1x*v1x + v1y*v1y + v1z*v1z)
        d = sqrt(v2x*v2x + v2y*v2y + v2z*v2z)
        dot_val = v2x*(-v1x) + v2y*(-v1y) + v2z*(-v1z)
        if nl * d < 1e-15:
            Z[ind, 1] = 0.0
        else:
            tmp = dot_val / (d * nl)
            if tmp > 1.0: tmp = 1.0
            elif tmp < -1.0: tmp = -1.0
            Z[ind, 1] = acos(tmp)

        # 二面角
        # vectorsCrutch2 = coords(i) - coords(j) = v1' (注意符号: Octave代码用的是 coords(format(ind,1)) - coords(format(ind,2)))
        # 实际: vectorsCrutch2 = c(i) - c(j) = [v1x, v1y, v1z] (已定义)
        # vectorsCrutch3 = c(j) - c(k) where k = format(ind,3)
        if ind == 2:
            v3x = 0.0; v3y = 0.0; v3z = -1.0
        else:
            k = fmt[ind, 2] - 1
            v3x = c[j, 0] - c[k, 0]; v3y = c[j, 1] - c[k, 1]; v3z = c[j, 2] - c[k, 2]

        # normalVEC1 = cross(vectorsCrutch2, vectorsCrutch3) = cross(v1, v3)
        n1x = v1y*v3z - v3y*v1z
        n1y = v1z*v3x - v3z*v1x
        n1z = v1x*v3y - v3x*v1y

        # normalVEC2 = cross(vectorsCrutch1, vectorsCrutch2) = cross(v2, v1)
        n2x = v2y*v1z - v1y*v2z
        n2y = v2z*v1x - v1z*v2x
        n2z = v2x*v1y - v1x*v2y

        # torsang = atan2(norm(vectorsCrutch2) * dot(vectorsCrutch1, normalVEC1), dot(normalVEC2, normalVEC1))
        nl = sqrt(v1x*v1x + v1y*v1y + v1z*v1z)  # |vectorsCrutch2|
        dot_val = v2x*n1x + v2y*n1y + v2z*n1z    # dot(vectorsCrutch1, normalVEC1)
        tmp = n2x*n1x + n2y*n1y + n2z*n1z         # dot(normalVEC2, normalVEC1)

        torsang = atan2(nl * dot_val, tmp)
        Z[ind, 2] = -torsang
        if cabs(torsang + M_PI) < 0.0001:
            Z[ind, 2] = M_PI

    return Zmatrix


# ═══════════════════════════════════════════════════════════════
# 4. find_pair_fast — 键连关系查找
# ═══════════════════════════════════════════════════════════════

def find_pair_fast(f64_t[:, ::1] coor, f64_t[:] radii):
    """查找键连关系 (替代 find_pair.m)"""
    cdef int n = coor.shape[0]
    cdef int N_max = 6
    cdef cnp.ndarray[cnp.int32_t, ndim=2] Pair_arr = np.zeros((n, N_max), dtype=np.int32)
    cdef int[:, ::1] P = Pair_arr
    cdef int i, j, cnt_i, cnt_j
    cdef f64_t dx, dy, dz, dist, thresh

    for i in range(n - 1):
        for j in range(i + 1, n):
            dx = coor[i, 0] - coor[j, 0]
            dy = coor[i, 1] - coor[j, 1]
            dz = coor[i, 2] - coor[j, 2]
            dist = sqrt(dx*dx + dy*dy + dz*dz)
            thresh = 1.2 * (radii[i] + radii[j])
            if dist < thresh:
                cnt_i = P[i, N_max - 1]
                if cnt_i < N_max - 2:
                    P[i, cnt_i] = j + 1
                P[i, N_max - 1] = cnt_i + 1
                cnt_j = P[j, N_max - 1]
                if cnt_j < N_max - 2:
                    P[j, cnt_j] = i + 1
                P[j, N_max - 1] = cnt_j + 1

    return Pair_arr


# ═══════════════════════════════════════════════════════════════
# 5. smart_rot_inertia — 完整二面角迭代循环
# ═══════════════════════════════════════════════════════════════

def smart_rot_inertia(f64_t[:, ::1] coords,
                       cnp.int32_t[:, ::1] fmt,
                       int num_optFlags,
                       cnp.int32_t[:, ::1] flex_dihedral,
                       cnp.int32_t[:, ::1] CN_target,
                       f64_t[:] radii,
                       f64_t[:] inertia_eigvals,
                       f64_t[:, ::1] inertia_eigvecs,
                       f64_t angleMax,
                       f64_t transMax,
                       int maxIter,
                       seed=None):
    """
    SmartRotInertia 完整加速版本

    返回: (Zmatrix, MOLCOORS, goodRot, n_iter)
    """
    import numpy as rng_np
    if seed is not None:
        rng_np.random.seed(seed)

    cdef int N = coords.shape[0]
    cdef int i, j, loop, ind, it
    cdef f64_t angle, Ref, angleStep, tx, ty, tz
    cdef bint goodRot = 0
    cdef int n_iter = 0
    cdef int i1, i2
    cdef int N_max = 6
    cdef cnp.ndarray[f64_t, ndim=2] MOLCOORS, Zmatrix
    cdef cnp.ndarray[cnp.int32_t, ndim=2] CN_arr
    cdef f64_t[:, ::1] MC, Z
    cdef int[:, ::1] CN_view
    cdef cnp.int32_t[:, ::1] CN_t

    # 1. 惯量主轴旋转
    MOLCOORS = np.array(coords, dtype=np.float64, copy=True)
    MC = MOLCOORS

    if cabs(inertia_eigvals[0]) < 0.0001:
        Ref = inertia_eigvals[1]
    else:
        Ref = inertia_eigvals[0]

    for loop in range(3):
        if inertia_eigvals[loop] > 0.0001:
            angle = (rng_np.random.random() - 0.5) * 2 * angleMax * Ref / inertia_eigvals[loop]
            _rotate_rigid_body_inplace(MC, inertia_eigvecs, loop, angle)

    # 2. 随机平移
    tx = (rng_np.random.random() - 0.5) * 2 * transMax
    ty = (rng_np.random.random() - 0.5) * 2 * transMax
    tz = (rng_np.random.random() - 0.5) * 2 * transMax
    for i in range(N):
        MC[i, 0] += tx; MC[i, 1] += ty; MC[i, 2] += tz

    # 3. 转回 Z-matrix
    Zmatrix = coords_to_zmatrix(MOLCOORS, fmt)
    Z = Zmatrix

    # 4. 二面角迭代
    if num_optFlags > 0:
        CN_t = CN_target

        for it in range(maxIter):
            n_iter = it + 1
            for ind in range(num_optFlags):
                i1 = flex_dihedral[ind, 0] - 1
                i2 = flex_dihedral[ind, 1] - 1
                if it < maxIter // 3:
                    angleStep = (rng_np.random.random() - 0.5) * M_PI
                elif it < 2 * maxIter // 3:
                    angleStep = (rng_np.random.random() - 0.5) * M_PI / 3
                else:
                    angleStep = (rng_np.random.random() - 0.5) * M_PI / 6
                Z[i1, i2] += angleStep

            MOLCOORS = zmatrix_to_coords(Zmatrix, fmt)
            MC = MOLCOORS
            Zmatrix = coords_to_zmatrix(MC, fmt)
            Z = Zmatrix

            CN_arr = find_pair_fast(MC, radii)
            CN_view = CN_arr

            goodRot = 1
            for i in range(N):
                for j in range(N_max):
                    if CN_view[i, j] != CN_t[i, j]:
                        goodRot = 0
                        break
                if not goodRot:
                    break

            if goodRot:
                break

        if not goodRot:
            MC = np.array(coords, dtype=np.float64, copy=True)
            for loop in range(3):
                if inertia_eigvals[loop] > 0.0001:
                    angle = (rng_np.random.random() - 0.5) * 2 * angleMax * Ref / inertia_eigvals[loop]
                    _rotate_rigid_body_inplace(MC, inertia_eigvecs, loop, angle)
            tx = (rng_np.random.random() - 0.5) * 2 * transMax
            ty = (rng_np.random.random() - 0.5) * 2 * transMax
            tz = (rng_np.random.random() - 0.5) * 2 * transMax
            for i in range(N):
                MC[i, 0] += tx; MC[i, 1] += ty; MC[i, 2] += tz
            Zmatrix = coords_to_zmatrix(MC, fmt)
            Z = Zmatrix

    return Zmatrix, MOLCOORS, 1 if goodRot else 0, n_iter


cdef inline void _rotate_rigid_body_inplace(f64_t[:, ::1] coords, f64_t[:, ::1] eigvecs,
                                             int axis_idx, f64_t angle) noexcept nogil:
    """绕质心绕指定本征向量旋转"""
    cdef int N = coords.shape[0]
    cdef int i
    cdef f64_t cx = 0.0, cy = 0.0, cz = 0.0
    cdef f64_t ax, ay, az, ax2, ay2, az2
    cdef f64_t c_val = cos(angle), s_val = sin(angle), omc = 1.0 - c_val
    cdef f64_t x, y, z

    for i in range(N):
        cx += coords[i, 0]; cy += coords[i, 1]; cz += coords[i, 2]
    cx /= N; cy /= N; cz /= N

    ax = eigvecs[axis_idx, 0]; ay = eigvecs[axis_idx, 1]; az = eigvecs[axis_idx, 2]
    ax2 = ax*ax; ay2 = ay*ay; az2 = az*az

    for i in range(N):
        x = coords[i, 0] - cx; y = coords[i, 1] - cy; z = coords[i, 2] - cz
        coords[i, 0] = (c_val + ax2*omc)*x + (ax*ay*omc - az*s_val)*y + (ax*az*omc + ay*s_val)*z + cx
        coords[i, 1] = (ay*ax*omc + az*s_val)*x + (c_val + ay2*omc)*y + (ay*az*omc - ax*s_val)*z + cy
        coords[i, 2] = (az*ax*omc - ay*s_val)*x + (az*ay*omc + ax*s_val)*y + (c_val + az2*omc)*z + cz
