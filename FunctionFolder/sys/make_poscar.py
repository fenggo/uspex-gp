#!/usr/bin/env python3
"""生成测试用 POSCAR 文件"""
import numpy as np

data = np.load('/tmp/test_structure.npz')
lattice = data['lattice']
coords = data['coords']
numIons = data['numIons']

with open('/tmp/test_POSCAR', 'w') as f:
    f.write("EA1    20.350 20.988 21.074 90.32 90.02 89.93 Sym.group:    1\n")
    f.write("1.0000\n")
    for i in range(3):
        f.write(f"  {lattice[i,0]:12.6f}  {lattice[i,1]:12.6f}  {lattice[i,2]:12.6f}\n")
    f.write("   C   O   N   H\n")
    f.write(f"  {numIons[0]}   {numIons[1]}   {numIons[2]}   {numIons[3]} \n")
    f.write("Direct\n")
    for i in range(coords.shape[0]):
        f.write(f"  {coords[i,0]:12.6f}  {coords[i,1]:12.6f}  {coords[i,2]:12.6f}\n")

print("Saved /tmp/test_POSCAR")
