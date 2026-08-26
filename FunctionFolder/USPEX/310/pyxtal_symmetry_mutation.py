#!/usr/bin/env python3
"""
pyxtal_symmetry_mutation.py — PyXtal 对称性引导变异算子

策略: 高对称空间群生成 → 降为 P1 → 小扰动破缺对称性
产生保持分子 packing 但对称性不同的新结构

算法:
  1. 在高对称空间群 (如 P2₁/c, P-1, P2₁2₁2₁) 生成父代
  2. 将分子坐标提取为 P1 表示 (去对称化)
  3. 对每个分子施加小随机位移 + 旋转 (破缺残留对称性)
  4. 距离检查 → 输出变体

用法:
  python3 pyxtal_symmetry_mutation.py --mols MOL_1,MOL_2 --numMols 4,4 --volume 6000 --outdir .
"""

import argparse, json, os, sys, time
import numpy as np


def read_uspex_mol(mol_file):
    with open(mol_file, 'r') as f:
        lines = [l.rstrip() for l in f if l.strip()]
    name = lines[0].strip()
    n_atoms = int(lines[1].split(':')[1].strip())
    symbols, positions = [], []
    for i in range(2, 2 + n_atoms):
        parts = lines[i].split()
        symbols.append(parts[0])
        positions.append([float(parts[1]), float(parts[2]), float(parts[3])])
    arr = np.array(positions, dtype=float)
    arr -= arr.mean(axis=0)
    return name, symbols, arr


def generate_in_spg(xtal_mols, num_each, spg, factor, seed):
    """在指定空间群生成分子晶体"""
    from pyxtal.molecular_crystal import molecular_crystal
    from pyxtal.symmetry import Group

    for attempt in range(30):
        try:
            s = seed + attempt * 1009 if seed is not None else None
            crystal = molecular_crystal(
                dim=3, group=Group(spg),
                molecules=xtal_mols, numMols=list(num_each),
                factor=factor, seed=s,
            )
            if getattr(crystal, 'valid', False):
                return crystal
        except Exception:
            continue
    return None


def extract_molecules(crystal, xtal_mols, num_each):
    """从 PyXtal crystal 提取所有独立分子坐标 (P1 表示)

    PyXtal 的 mol_sites 中每个 site 包含一个 Wyckoff 位置的所有对称等价分子。
    用 get_coords_and_species 获取全部原子后按分子原子数拆分。
    """
    lattice = np.array(crystal.lattice.matrix)

    # 构建分子原子数映射
    mol_natoms = {}
    for i, (n, mol) in enumerate(zip(num_each, xtal_mols)):
        mol_natoms[i] = len(mol.mol)

    mols = []
    for site in crystal.mol_sites:
        coords, species = site.get_coords_and_species(absolute=True, unitcell=True)
        mol_type = site.type
        n_per_mol = mol_natoms[mol_type]
        n_total = len(species)
        n_mols = n_total // n_per_mol

        for i in range(n_mols):
            s = i * n_per_mol
            e = s + n_per_mol
            mols.append({
                'coords': coords[s:e],
                'species': list(species[s:e]),
            })

    return mols, lattice


def min_distance(coords_a, coords_b, lattice, inv_lattice):
    """两分子间最小原子距离"""
    min_d = float('inf')
    for ca in coords_a:
        for cb in coords_b:
            diff = (ca - cb) @ inv_lattice
            diff = diff - np.round(diff)
            d = np.linalg.norm(diff @ lattice)
            if d < min_d:
                min_d = d
    return min_d


def perturb_molecules(mols, lattice, max_trans=1.0, max_rot_deg=10.0, min_dist=1.5):
    """
    对分子施加小扰动破缺对称性
    - 随机平移 (max_trans Å)
    - 随机旋转 (max_rot_deg 度)
    - 距离检查
    """
    from scipy.spatial.transform import Rotation as R
    inv = np.linalg.inv(lattice)

    new_mols = []
    for i, mol in enumerate(mols):
        best_coords = None
        best_score = -float('inf')

        for _ in range(20):
            # 随机旋转
            axis = np.random.randn(3)
            axis /= np.linalg.norm(axis)
            angle = np.random.uniform(-max_rot_deg, max_rot_deg) * np.pi / 180
            rot = R.from_rotvec(axis * angle).as_matrix()

            # 随机平移
            trans = np.random.uniform(-max_trans, max_trans, 3)

            center = mol['coords'].mean(axis=0)
            rotated = (mol['coords'] - center) @ rot.T + center
            test_coords = rotated + trans

            # 距离检查 (与已放置分子)
            ok = True
            min_d = float('inf')
            for pm in new_mols:
                d = min_distance(test_coords, pm['coords'], lattice, inv)
                if d < min_dist:
                    ok = False
                    break
                if d < min_d:
                    min_d = d

            if ok:
                score = min_d
                if score > best_score:
                    best_score = score
                    best_coords = test_coords.copy()

        if best_coords is None:
            best_coords = mol['coords'].copy()  # 保持原位

        new_mols.append({'coords': best_coords, 'species': mol['species']})

    return new_mols


def write_output(mols, lattice, outdir, variant_name):
    """输出 POSCAR + optimized.structure + output"""
    os.makedirs(outdir, exist_ok=True)

    all_syms, all_pos = [], []
    for mol in mols:
        all_syms.extend(mol['species'])
        all_pos.extend(mol['coords'])
    all_pos = np.array(all_pos)
    inv = np.linalg.inv(lattice)

    # POSCAR
    order = ['C', 'O', 'N', 'H']
    idx = [i for s in order for i, sym in enumerate(all_syms) if sym == s]
    counts = [sum(1 for s in all_syms if s == o) for o in order]

    with open(os.path.join(outdir, 'POSCAR'), 'w') as f:
        f.write(' '.join(order) + '\n1.0\n')
        for row in lattice:
            f.write(f'  {row[0]:16.10f}  {row[1]:16.10f}  {row[2]:16.10f}\n')
        f.write(' '.join(order) + '\n')
        f.write(' '.join(str(c) for c in counts) + '\nDirect\n')
        for i in idx:
            frac = all_pos[i] @ inv - np.floor(all_pos[i] @ inv)
            f.write(f'  {frac[0]:16.12f}  {frac[1]:16.12f}  {frac[2]:16.12f}\n')

    # GULP
    a = np.linalg.norm(lattice[0]); b = np.linalg.norm(lattice[1]); c = np.linalg.norm(lattice[2])
    alpha = np.arccos(np.dot(lattice[1], lattice[2]) / (b*c)) * 180/np.pi
    beta  = np.arccos(np.dot(lattice[0], lattice[2]) / (a*c)) * 180/np.pi
    gamma = np.arccos(np.dot(lattice[0], lattice[1]) / (a*b)) * 180/np.pi
    scaled = all_pos @ inv

    with open(os.path.join(outdir, 'optimized.structure'), 'w') as f:
        f.write('opti nosymmetry conp qiterative conjugate\n\ncell\n')
        f.write(f' {a:12.8f} {b:12.8f} {c:12.8f} {alpha:12.8f} {beta:12.8f} {gamma:12.8f}\n')
        f.write('fractional  1\n')
        for sym, pos in zip(all_syms, scaled):
            f.write(f'{sym:<2s}    core {pos[0]:14.9f} {pos[1]:14.9f} {pos[2]:14.9f}    0.0 1.0 0.0\n')
        f.write('\ndump every      1 optimized.structure\n')

    with open(os.path.join(outdir, 'output'), 'w') as f:
        f.write('  Cycle:      0 Energy:             0.00000000\n')

    # molecules.txt (MATLAB 读取)
    with open(os.path.join(outdir, 'molecules.txt'), 'w') as f:
        f.write(f"{len(mols)} 0\n")
        for row in lattice:
            f.write(f"{row[0]:.12f} {row[1]:.12f} {row[2]:.12f}\n")
        for mol in mols:
            f.write(f"{len(mol['coords'])} 0\n")
            for c in mol['coords']:
                f.write(f"{c[0]:.12f} {c[1]:.12f} {c[2]:.12f}\n")

    meta = {
        'method': 'pyxtal_symmetry_mutation',
        'variant': variant_name,
        'volume': float(abs(np.linalg.det(lattice))),
        'atoms': len(all_syms),
        'num_molecules': len(mols),
    }
    with open(os.path.join(outdir, 'crystal_meta.json'), 'w') as f:
        json.dump(meta, f, indent=2)

    return meta


def main():
    parser = argparse.ArgumentParser(description='PyXtal Symmetry Mutation')
    parser.add_argument('--mols', required=True)
    parser.add_argument('--numMols', required=True)
    parser.add_argument('--volume', type=float, default=6000.0)
    parser.add_argument('--spg', type=int, default=14, help='父代空间群')
    parser.add_argument('--n-variants', type=int, default=3, help='变体数量')
    parser.add_argument('--max-trans', type=float, default=1.5, help='最大平移 (Å)')
    parser.add_argument('--max-rot', type=float, default=15.0, help='最大旋转 (度)')
    parser.add_argument('--min-dist', type=float, default=1.5, help='最小原子间距')
    parser.add_argument('--seed', type=int, default=None)
    parser.add_argument('--outdir', default='symmetry_variants')
    args = parser.parse_args()

    mol_files = [x.strip() for x in args.mols.split(',')]
    num_each = [int(x.strip()) for x in args.numMols.split(',')]

    from pymatgen.core import Molecule
    from pyxtal.molecule import pyxtal_molecule

    xtal_mols = []
    for mf in mol_files:
        if not os.path.exists(mf):
            print(f"Error: {mf} not found", file=sys.stderr)
            sys.exit(1)
        name, syms, pos = read_uspex_mol(mf)
        xtal_mols.append(pyxtal_molecule(Molecule(syms, pos), symmetrize=False))

    mol_vol = sum(n * m.volume for n, m in zip(num_each, xtal_mols))
    factor = args.volume / mol_vol if mol_vol > 0 else 1.5

    t0 = time.time()

    # Step 1: 高对称空间群生成父代
    print(f"Generating parent in SPG #{args.spg}...", flush=True)
    parent = generate_in_spg(xtal_mols, num_each, args.spg, factor, args.seed)
    if parent is None:
        print("Failed to generate parent", file=sys.stderr)
        sys.exit(1)
    print(f"  Parent: SPG {parent.group.symbol} (#{parent.group.number}), "
          f"Vol={parent.lattice.volume:.1f} Å³, Sites={len(parent.mol_sites)}")

    # Step 2: 提取分子 (P1 表示)
    mols, lattice = extract_molecules(parent, xtal_mols, num_each)
    print(f"  Extracted {len(mols)} molecules in P1")

    # Step 3: 生成多个扰动变体
    os.makedirs(args.outdir, exist_ok=True)
    for vi in range(args.n_variants):
        vseed = (args.seed or 0) + vi * 10007
        np.random.seed(vseed)

        perturbed = perturb_molecules(mols, lattice,
                                      max_trans=args.max_trans,
                                      max_rot_deg=args.max_rot,
                                      min_dist=args.min_dist)

        vdir = os.path.join(args.outdir, f'variant_{vi:03d}')
        meta = write_output(perturbed, lattice, vdir, f'variant_{vi}')
        print(f"  Variant {vi}: Atoms={meta['atoms']}, Vol={meta['volume']:.1f} Å³")

    print(f"\nDone: {args.n_variants} variants, time={time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
