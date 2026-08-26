#!/usr/bin/env python
"""
uspex_rotation_wrapper.py — Octave ↔ Python 桥接：并行旋转变异

被 Octave EA_310.m 通过 system() 调用，
并行执行多个 Rotation_310 变异操作。

用法:
  python3 uspex_rotation_wrapper.py --input=batch.json --output=result.mat

输入 JSON:
  {
    "structures": [
      {
        "coords": [[x,y,z], ...],      # 分子坐标 (8 molecules, each N×3)
        "format": [[i,j,k], ...],      # Z矩阵格式 (per molecule type)
        "flex_dihedral": [[i1,i2], ...],
        "CN_target": [[...]],           # 目标键连关系
        "radii": [...],                 # 共价半径
        "num_optFlags": int,
        "angleMax": float,
        "transMax": float,
        "maxIter": int
      },
      ...
    ]
  }

输出 .mat:
  {
    "results": [  # list of struct arrays
      {"MOLCOORS": [...], "ZMATRIX": [...], "goodRot": int, "n_iter": int},
      ...
    ]
  }
"""

import sys
import os
import json
import time
import argparse
import numpy as np
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uspex_rotation_core as urc


def process_one_structure(args):
    """处理一个结构的旋转变异（工作进程）"""
    idx, struct_data, seed = args

    coords = np.array(struct_data['coords'], dtype=np.float64)
    fmt = np.array(struct_data['format'], dtype=np.int32)
    flex_dihedral = np.array(struct_data['flex_dihedral'], dtype=np.int32) if struct_data['flex_dihedral'] else np.zeros((0, 2), dtype=np.int32)
    CN_target = np.array(struct_data['CN_target'], dtype=np.int32) if struct_data['CN_target'] else np.zeros(0, dtype=np.int32)
    radii = np.array(struct_data['radii'], dtype=np.float64)
    inertia_eigvals = np.array(struct_data['inertia_eigvals'], dtype=np.float64)
    inertia_eigvecs = np.array(struct_data['inertia_eigvecs'], dtype=np.float64)

    num_optFlags = struct_data['num_optFlags']
    angleMax = struct_data.get('angleMax', np.pi / 4)
    transMax = struct_data.get('transMax', 0.3)
    maxIter = struct_data.get('maxIter', 200)

    Zmatrix, MOLCOORS, goodRot, n_iter = urc.smart_rot_inertia(
        coords, fmt, num_optFlags, flex_dihedral, CN_target,
        radii, inertia_eigvals, inertia_eigvecs,
        angleMax, transMax, maxIter, seed=seed
    )

    return {
        'idx': idx,
        'MOLCOORS': MOLCOORS,
        'ZMATRIX': Zmatrix,
        'goodRot': int(goodRot),
        'n_iter': n_iter,
    }


def main():
    ap = argparse.ArgumentParser(description='USPEX parallel rotation mutation')
    ap.add_argument('--input', type=str, required=True, help='Input JSON file')
    ap.add_argument('--output', type=str, required=True, help='Output .mat file')
    ap.add_argument('--nproc', type=int, default=None, help='Number of processes')
    args = ap.parse_args()

    t_start = time.time()

    with open(args.input, 'r') as f:
        batch = json.load(f)

    structures = batch['structures']
    n = len(structures)
    nproc = args.nproc or min(n, os.cpu_count() or 4)

    # 准备参数
    tasks = [(i, s, i * 1000 + 42) for i, s in enumerate(structures)]

    # 并行处理
    if nproc > 1 and n > 1:
        with Pool(nproc) as pool:
            results = pool.map(process_one_structure, tasks)
    else:
        results = [process_one_structure(t) for t in tasks]

    t_end = time.time()

    # 按 idx 排序
    results.sort(key=lambda r: r['idx'])

    # 写 .mat 文件
    from scipy.io import savemat

    # 将结果组织为 Octave 友好的格式
    n_mol = len(results)
    MOLCOORS_all = np.array([r['MOLCOORS'] for r in results])  # (n_structs, N, 3)
    ZMATRIX_all = np.array([r['ZMATRIX'] for r in results])
    goodRot_all = np.array([r['goodRot'] for r in results], dtype=np.int32)
    n_iter_all = np.array([r['n_iter'] for r in results], dtype=np.int32)

    savemat(args.output, {
        'MOLCOORS_all': MOLCOORS_all,
        'ZMATRIX_all': ZMATRIX_all,
        'goodRot_all': goodRot_all,
        'n_iter_all': n_iter_all,
        'n_structs': n,
        'time_total': t_end - t_start,
    })

    n_good = sum(goodRot_all)
    # print(f"# Rotation: {n} structs, {n_good} good, {t_end-t_start:.2f}s with {nproc} procs -> {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
