#!/usr/bin/env python
"""
pyxtal_random.py — PyXtal random_crystal 生成分子中心位置

替代 random_cell_mol (Fortran 二进制):
  - PyXtal random_crystal 将分子看作"原子"生成晶体结构
  - 输出所有分子中心位置（分数坐标）和晶格
  - Octave 端从 MOL 文件放置分子（保证原子顺序正确）

用法:
  python3 pyxtal_random.py --mols MOL_1,MOL_2 --numMols 4,4 --volume 6000 --spg 14

输出格式 (纯文本):
  n_mols n_types
  lattice (3×3)
  每个分子: type type_idx frac_x frac_y frac_z
"""

import argparse
import os
import sys
import numpy as np


def generate_molecule_centers(num_mols, spg, volume, factor=1.5,
                                max_attempts=50, seed=None):
    """
    用 PyXtal random_crystal 生成分子中心。

    每种分子类型用不同"原子"元素区分:
      TNT(4) + CL20(4) → species=['C','N'], numIons=[4,4]
    """
    if seed is not None:
        np.random.seed(seed)

    from pyxtal.crystal import random_crystal
    from pyxtal.symmetry import Group

    n_types = len(num_mols)
    species_map = ['C', 'N', 'O', 'F', 'S', 'P', 'Cl', 'Br', 'I', 'H']
    species = [species_map[i % len(species_map)] for i in range(n_types)]

    spg_fallback = [1, 2, 14, 19, 33, 61]
    trial_spgs = [spg] + [s for s in spg_fallback if s != spg]

    for sg in trial_spgs:
        try:
            group = Group(sg)
        except Exception:
            continue

        for attempt in range(max_attempts):
            try:
                s = seed + attempt * 1009 + sg * 10007 if seed is not None else None
                if s is not None:
                    np.random.seed(s)

                crystal = random_crystal(
                    dim=3, group=group,
                    species=species, numIons=list(num_mols),
                    factor=factor,
                )

                if not getattr(crystal, 'valid', False):
                    continue

                # 提取所有原子位置（按物种分组）
                return _extract_centers(crystal, species, num_mols)

            except Exception:
                continue

    return None


def _extract_centers(crystal, species, num_mols):
    """提取所有分子中心位置，按物种分组"""
    species_to_type = {s: i for i, s in enumerate(species)}
    lattice = np.array(crystal.lattice.matrix, dtype=float)

    # 按物种分组收集所有分数坐标
    centers_by_type = {i: [] for i in range(len(species))}

    for site in crystal.atom_sites:
        wp = site.wp
        type_idx = species_to_type.get(site.specie, 0)
        all_pos = wp.get_all_positions(site.position)
        for pos in all_pos:
            centers_by_type[type_idx].append(pos)

    return centers_by_type, lattice


def main():
    parser = argparse.ArgumentParser(
        description='PyXtal random_crystal → USPEX molecule centers')
    parser.add_argument('--mols', required=True,
                        help='MOL file paths, comma separated')
    parser.add_argument('--numMols', required=True,
                        help='molecule counts, comma separated')
    parser.add_argument('--volume', type=float, default=6000.0)
    parser.add_argument('--spg', type=int, default=1)
    parser.add_argument('--max-attempts', type=int, default=50)
    parser.add_argument('--seed', type=int, default=None)
    parser.add_argument('--outdir', default='.')
    args = parser.parse_args()

    mol_files = [x.strip() for x in args.mols.split(',')]
    num_each = [int(x.strip()) for x in args.numMols.split(',')]

    # 计算 factor: PyXtal 用 factor 缩放原子体积
    # 对于 8 个伪原子，factor=80 对应 ~6000-9000 Å³
    # 公式: factor = V / 75 (经验校准)
    total_mols = sum(num_each)
    if total_mols > 0:
        factor = args.volume / 75.0
    else:
        factor = 80.0
    factor = max(10.0, min(factor, 300.0))

    os.makedirs(args.outdir, exist_ok=True)

    result = generate_molecule_centers(
        num_mols=num_each, spg=args.spg, volume=args.volume,
        factor=factor, max_attempts=args.max_attempts, seed=args.seed,
    )

    if result is None:
        outpath = os.path.join(args.outdir, 'output')
        with open(outpath, 'w') as f:
            f.write('  Cycle:      0 Energy:         100000.00000000\n')
        sys.exit(1)

    centers_by_type, lattice = result

    # 写入数据文件
    data_path = os.path.join(args.outdir, 'pyxtal_data.txt')
    total_mols = sum(len(v) for v in centers_by_type.values())
    n_types = len(num_each)

    with open(data_path, 'w') as f:
        f.write(f"{total_mols} {n_types}\n")
        # Lattice
        for i in range(3):
            f.write(f"{lattice[i,0]:.12f} {lattice[i,1]:.12f} {lattice[i,2]:.12f}\n")
        # 每个分子: type_idx frac_x frac_y frac_z
        # 按类型顺序输出（TNT在前，CL20在后）
        for type_idx in range(n_types):
            for pos in centers_by_type[type_idx]:
                f.write(f"{type_idx} {pos[0]:.12f} {pos[1]:.12f} {pos[2]:.12f}\n")

    # 成功标记
    outpath = os.path.join(args.outdir, 'output')
    with open(outpath, 'w') as f:
        f.write('  Cycle:      0 Energy:              0.00000000\n')

    print(f"PyXtal: {total_mols} molecules, "
          f"SPG #{args.spg}, vol={np.linalg.det(lattice):.1f} Å³")


if __name__ == '__main__':
    main()