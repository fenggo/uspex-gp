#!/usr/bin/env python
"""
uspex_fingerprint_wrapper.py — Octave ↔ Python 桥接接口 v3

支持两种模式:
  1. 普通指纹（对应 ReadJobs_310.m 行51-52：分子质心指纹）
  2. 带分子内距离过滤的指纹（对应行58-64：全原子指纹）

用 scipy.io.savemat 写 .mat 文件，Octave 可直接 load。

用法:
  python3 uspex_fingerprint_wrapper.py --poscar=POSCAR --output=result.mat
  python3 uspex_fingerprint_wrapper.py --poscar=POSCAR --output=result.mat --intra-map=intra_map.npy
"""

import sys
import os
import time
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uspex_fast_core as ufc


def parse_poscar(filename):
    """解析 POSCAR 文件"""
    with open(filename, 'r') as f:
        lines = f.readlines()

    scale = float(lines[1].strip()) if lines[1].strip() else 1.0
    lattice = np.array([
        [float(x) for x in lines[2].split()],
        [float(x) for x in lines[3].split()],
        [float(x) for x in lines[4].split()],
    ], dtype=np.float64) * scale

    symbol_to_z = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Na': 11, 'Mg': 12,
                   'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'K': 19, 'Ca': 20,
                   'Ti': 22, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29,
                   'Zn': 30, 'Br': 35, 'I': 53}
    atom_symbols = lines[5].split()
    atomType = np.array([symbol_to_z.get(s, 0) for s in atom_symbols], dtype=np.int32)
    numIons = np.array([int(x) for x in lines[6].split()], dtype=np.int32)
    N = sum(numIons)

    coord_start = 8
    coords = np.zeros((N, 3), dtype=np.float64)
    for i in range(N):
        parts = lines[coord_start + i].split()
        coords[i, 0] = float(parts[0])
        coords[i, 1] = float(parts[1])
        coords[i, 2] = float(parts[2])

    return lattice, coords, numIons, atomType


def main():
    ap = argparse.ArgumentParser(description='USPEX fast fingerprint (Octave bridge v3)')
    ap.add_argument('--poscar', type=str, required=True)
    ap.add_argument('--output', type=str, default='fingerprint_result.mat')
    ap.add_argument('--rmax', type=float, default=12.0)
    ap.add_argument('--sigma', type=float, default=0.05)
    ap.add_argument('--delta', type=float, default=0.08)
    ap.add_argument('--dimension', type=int, default=3)
    ap.add_argument('--intra-map', type=str, default=None,
                    help='Intra_map .mat/.npy 文件（分子内距离标记矩阵），'
                         '用于过滤分子内距离')
    ap.add_argument('--soap', action='store_true',
                    help='同时计算 SOAP 指纹 (需要 dscribe)')
    ap.add_argument('--soap-r-cut', type=float, default=6.0)
    ap.add_argument('--soap-n-max', type=int, default=8)
    ap.add_argument('--soap-l-max', type=int, default=6)
    args = ap.parse_args()

    t_start = time.time()

    lattice, coords, numIons, atomType = parse_poscar(args.poscar)

    # 计算距离矩阵
    dist_arr, cc_idx, bc_idx, ti_arr, tj_arr, shift_arr, N_out, V, N = ufc.build_distance_matrix(
        coords, lattice, args.rmax, numIons, atomType
    )

    # 如果提供了 intra_map，过滤分子内距离
    # 注意：Octave 只对基本晶胞内的原子对（零位移）应用 Intra_map 过滤，
    # 周期镜像原子对的距离即使属于同一分子也保留。
    # 这对应 Octave ReadJobs_310.m 中:
    #   tmp = dist_matrix(1:sum(numIons), :);
    #   tmp(find(Intra_map==0)) = 0;
    #   dist_matrix(1:sum(numIons),:) = tmp;
    if args.intra_map and os.path.exists(args.intra_map):
        # 支持 .mat (Octave save) 和 .npy 两种格式
        if args.intra_map.endswith('.mat'):
            from scipy.io import loadmat
            intra_map = loadmat(args.intra_map)['Intra_map']
        else:
            intra_map = np.load(args.intra_map)
        intra_map = np.asarray(intra_map, dtype=np.int8)
        # 只对零位移的 pair 应用 Intra_map 过滤（匹配 Octave 行为）
        # 零位移 pair: shift_arr == 0, 且 intra_map[cc, bc] == 0 -> 过滤掉
        # 非零位移 pair: 无论 intra_map 如何都保留
        zero_shift_mask = (shift_arr == 0)
        intra_mask = (intra_map[cc_idx, bc_idx] == 1)
        keep_mask = np.where(zero_shift_mask, intra_mask, True)
        dist_arr = dist_arr[keep_mask]
        cc_idx = cc_idx[keep_mask]
        bc_idx = bc_idx[keep_mask]
        ti_arr = ti_arr[keep_mask]
        tj_arr = tj_arr[keep_mask]

    # 计算指纹
    order, fing, atom_fing, _ = ufc.fingerprint_calc(
        dist_arr, cc_idx, bc_idx, ti_arr, tj_arr, N_out, V, numIons,
        args.rmax, args.sigma, args.delta, args.dimension
    )

    t_end = time.time()

    # 写 .mat 文件
    from scipy.io import savemat
    out = {
        'order': order,
        'FINGERPRINT': fing,
        'atom_fing': atom_fing,
        'V': V,
        'n_pairs': len(dist_arr),
        'time_total': t_end - t_start,
    }

    # 可选：计算 SOAP 指纹
    if args.soap:
        try:
            from ase.io import read as ase_read
            atoms = ase_read(args.poscar, format='vasp')
            # 优先使用 uspexkit 包，否则用内联实现
            try:
                sys.path.insert(0, '/home/feng/uspexkit')
                from uspexkit.soap import soap_fingerprint
                soap_fp = soap_fingerprint(atoms, r_cut=args.soap_r_cut,
                                           n_max=args.soap_n_max,
                                           l_max=args.soap_l_max)
            except ImportError:
                from dscribe.descriptors import SOAP as SOAPDesc
                species = sorted(set(atoms.get_chemical_symbols()))
                soap_obj = SOAPDesc(species=species, r_cut=args.soap_r_cut,
                                    n_max=args.soap_n_max, l_max=args.soap_l_max,
                                    periodic=True, sparse=False, average='inner')
                soap_fp = np.asarray(soap_obj.create(atoms))
            out['soap_fp'] = soap_fp
        except Exception as e:
            import traceback
            print(f"# SOAP computation failed: {e}", file=sys.stderr)
            traceback.print_exc(file=sys.stderr)

    savemat(args.output, out)

    # print(f"# Fingerprint: {t_end - t_start:.4f}s, pairs={len(dist_arr)} -> {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()