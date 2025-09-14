import os
from rdkit import Chem
from rdkit.Chem import Draw
from pathlib import Path
import re

def visualize_scps(frequent_subgraph_file, cmine_results_file, output_folder='scp_images'):
    os.makedirs(output_folder, exist_ok=True)

    # Step 1: Parse cmine_results.txt
    subgraph_coverage = {}
    with open(cmine_results_file, 'r') as f:
        for line in f:
            if line.startswith("Subgraph IDs:"):
                parts = line.strip().split("|")
                id_part = parts[0].split(":")[1].strip().strip("{}")
                # print(id_part)
                coverage_part = float(parts[1].split(":")[1].strip())
                # print(coverage_part)
                id_part = id_part.strip()
                # print(id_part)
                id_part = re.sub(r'\D', '', id_part)  # Removes anything that isn't a digit
                if id_part.isdigit():
                    # print(id_part)
                    subgraph_coverage[int(id_part)] = coverage_part

    # Limit to top 10 subgraphs
    top_subgraphs = sorted(subgraph_coverage.items(), key=lambda x: -x[1])[:10]
    # print(top_subgraphs)

    # Step 2: Parse frequentSubgraphs.txt into a dictionary of {id: (nodes, edges)}
    subgraph_dict = {}
    with open(frequent_subgraph_file, 'r') as f:
        lines = f.readlines()

    current_id = None
    nodes = []
    edges = []
    for line in lines:
        line = line.strip()
        if line.startswith("t #"):
            if current_id is not None:
                subgraph_dict[current_id] = (nodes, edges)
            parts = line.split()
            current_id = int(parts[2])
            nodes = []
            edges = []
        elif line.startswith("v"):
            _, idx, atom = line.split()
            nodes.append((int(idx), int(atom)))
        elif line.startswith("e"):
            _, src, dst, bond = line.split()
            edges.append((int(src), int(dst), int(bond)))
    if current_id is not None:
        subgraph_dict[current_id] = (nodes, edges)  # Add last subgraph

    # Step 3: Process each subgraph and visualize
    for sub_id, coverage in top_subgraphs:
        if sub_id not in subgraph_dict:
            print(f"Subgraph ID {sub_id} not found in frequent subgraphs.")
            continue

        nodes, edges = subgraph_dict[sub_id]
        mol = Chem.RWMol()
        idx_map = {}

        for node_id, atomic_num in nodes:
            idx_map[node_id] = mol.AddAtom(Chem.Atom(atomic_num))

        bond_type_dict = {
            1: Chem.BondType.SINGLE,
            2: Chem.BondType.DOUBLE,
            3: Chem.BondType.TRIPLE,
            4: Chem.BondType.AROMATIC
        }

        for src, dst, bond in edges:
            bond_type = bond_type_dict.get(bond, Chem.BondType.SINGLE)
            try:
                mol.AddBond(idx_map[src], idx_map[dst], bond_type)
            except:
                continue  # Skip invalid bonds

        final_mol = mol.GetMol()
        smiles = Chem.MolToSmiles(final_mol)

        # Step 4: Draw with coverage annotation
        img = Draw.MolToImage(final_mol, size=(300, 300), legend=f"ID: {sub_id}, Cov: {coverage:.2f}")
        img_path = os.path.join(output_folder, f"scp_{sub_id}.png")
        img.save(img_path)

    print(f"Saved top {len(top_subgraphs)} SCP visualizations to: {output_folder}")


