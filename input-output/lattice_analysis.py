import numpy as np
from ase import Atoms
from ase.visualize import view

column_datax = []
column_datay = []
column_dataz = []
r0 = []
timesteps = int(0.04/0.00001)
from ase.io import write

with open('lattice_output.txt', 'r') as f:
    for line in f:
        # Assuming space-separated values
        values = line.strip().split()
        # Assuming you want the 2nd column (index 1)
        column_datax.append(float(values[1]))
        column_datay.append(float(values[2]))
        column_dataz.append(0)

array_2d = np.array([column_datax, column_datay, column_dataz]).T # Transpose to get column-wise combination
print(array_2d) # Output: [[1 4] [2 5] [3 6]]

def extract_species_flag(filepath, species_index=0):
    values = []
    inside_config = False

    with open(filepath, 'r') as f:
        for line in f:
            stripped = line.strip()

            if stripped.startswith("configuration"):
                inside_config = True
                continue

            if inside_config:
                parts = stripped.split()
                if len(parts) >= 5 and parts[0].isdigit():
                    try:
                        # Surface species flags are usually at the end of line
                        value = float(parts[-3 + species_index])  # e.g., -3, -2, -1
                        values.append(value)
                    except (ValueError, IndexError):
                        continue

    return values

def extract_sites_flag(filepath, species_index=0):
    values = []
    inside_config = False

    with open(filepath, 'r') as f:
        for line in f:
            stripped = line.strip()

            if stripped.startswith("configuration"):
                inside_config = True
                continue

            if inside_config:
                parts = stripped.split()
                if len(parts) >= 5 and parts[0].isdigit():
                    try:
                        # Surface species flags are usually at the end of line
                        value = float(parts[1])  # e.g., -3, -2, -1
                        values.append(value)
                    except (ValueError, IndexError):
                        continue

    return values

vals = extract_species_flag("history_output.txt", species_index=0)  # GeH3*
sites = extract_sites_flag("history_output.txt", species_index=0)  # GeH3*
print(vals)

for i in range(timesteps):
    ats = len(column_datax)*['Ge']
    for j in range(int(10000*i),int(10000*(i+1))):
        if vals[j] == 1: #GeH3
            ats[int(sites[j])] = 'Au'
        elif vals[j] == 2: #H
            ats[int(sites[j])] = 'H'
        elif vals[j] == 3: #GeH2
            ats[int(sites[j])] = 'Co'


    atoms = Atoms(symbols=ats, positions=array_2d, cell = [4.082*50,4.269*50,10], pbc=True)
    write("POSCAR{}".format(i), atoms, format="vasp")

