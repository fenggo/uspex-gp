#!/usr/bin/env python
"""
uspex_softmode_wrapper.py — Octave ↔ Python 桥接：calcSoftModes 加速

被 Octave fast_calcSoftModes.m 通过 system() 调用。
只做一件事: 替代 calcSoftModes_molecules.m 的完整计算。

输入: .mat 文件 (coords, lat, at_types, R_val, N_val, val_arr, goodBonds)
输出: .mat 文件 (freq (3N×3N 对角矩阵), eigvector (3N×3N))

与 Octave calcSoftModes_molecules 返回格式完全一致:
  freq       = abs(real(diag(eigenvalues)))   (3N×3N 对角矩阵)
  eigvector  = real(eigenvectors)              (3N×3N)
"""

import sys
import os
import time
import numpy as np
from scipy.io import loadmat, savemat

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import uspex_softmode_core as usc


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', type=str, required=True)
    ap.add_argument('--output', type=str, required=True)
    args = ap.parse_args()

    t0 = time.time()

    data = loadmat(args.input)

    coords  = np.ascontiguousarray(data['coords'])              # (N, 3)
    lat     = np.ascontiguousarray(data['lat'])                 # (3, 3)
    at_types = np.ascontiguousarray(data['at_types'].flatten().astype(np.int32))
    R_val   = np.ascontiguousarray(data['R_val'].flatten())
    N_val   = np.ascontiguousarray(data['N_val_arr'].flatten())
    val_arr = np.ascontiguousarray(data['val_arr'].flatten())
    goodBonds = np.ascontiguousarray(data['goodBonds'])         # (n_types, n_types) 矩阵

    # 调用 Cython 核心 — 返回对角矩阵 (与 Octave eig(D) 一致)
    freq_mat, eigvector = usc.calc_soft_modes_raw(
        coords, lat, at_types, R_val, N_val, val_arr, goodBonds)

    savemat(args.output, {
        'freq': freq_mat,
        'eigvector': eigvector,
        'time_used': time.time() - t0,
    })

    # print(f"# calcSoftModes: {time.time()-t0:.3f}s -> {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
