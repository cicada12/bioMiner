import ast
import os
from rdkit import Chem
from rdkit.Chem import Draw

# RDKit-style bond type mapping
BOND_TYPE_DICT = {
    1: Chem.BondType.SINGLE,
    2: Chem.BondType.DOUBLE,
    3: Chem.BondType.TRIPLE,
    4: Chem.BondType.AROMATIC
}

def plot_subgraphs_from_txt(txt_file, image_dir, show_labels=True):
    os.makedirs(image_dir, exist_ok=True)

    # Clear existing images
    for file in os.listdir(image_dir):
        file_path = os.path.join(image_dir, file)
        if os.path.isfile(file_path) and file_path.endswith(('.png', '.jpg')):
            os.remove(file_path)

    # Read subgraph lines
    with open(txt_file, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if not line.strip():
            continue

        # Parse line
        _, subgraph_str = line.strip().split(":", 1)
        edges = ast.literal_eval(subgraph_str.strip())

        # Create an editable RDKit molecule
        mol = Chem.RWMol()
        idx_map = {}

        # Add nodes as atoms
        for edge in edges:
            v1, v2, v1_label, edge_label, v2_label = edge
            if v1 not in idx_map:
                idx_map[v1] = mol.AddAtom(Chem.Atom(int(v1_label)))
            if v2 not in idx_map:
                idx_map[v2] = mol.AddAtom(Chem.Atom(int(v2_label)))

        # Add bonds
        for edge in edges:
            v1, v2, _, bond_label, _ = edge
            bond_type = BOND_TYPE_DICT.get(bond_label, Chem.BondType.SINGLE)
            try:
                mol.AddBond(idx_map[v1], idx_map[v2], bond_type)
            except:
                continue  # skip invalid bonds

        final_mol = mol.GetMol()

        # Draw molecule as image
        legend_text = line.split(":")[0].strip() if show_labels else ""
        img = Draw.MolToImage(final_mol, size=(300, 300), legend=legend_text)
        img_path = os.path.join(image_dir, f"{legend_text.replace(' ', '_')}.png")
        img.save(img_path)
