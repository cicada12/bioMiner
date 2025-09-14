from rdkit import Chem

def smiles_to_gtd(smiles_file, output_path):
    """
    Convert a file of SMILES strings to graph transaction format used in graph mining.

    Each graph is formatted like:
        t # graph_id
        v node_id atom_number
        e source_id target_id bond_type

    Atom labels are converted to atomic numbers (integers), and bond types are encoded as:
        single=1, double=2, triple=3, aromatic=4

    Parameters:
        smiles_file (str): Path to the text file containing SMILES strings, one per line.
        output_path (str): Path to the output file to save the graph transaction format.

    Returns:
        None
    """
    output_lines = []
    graph_id = 0

    with open(smiles_file, 'r') as f:
        for line in f:
            smiles = line.strip()
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue  # Skip invalid SMILES

            lines = [f"t # {graph_id}"]
            graph_id += 1

            # Add atoms as vertices with atomic numbers
            for atom in mol.GetAtoms():
                atom_id = atom.GetIdx()
                atom_number = atom.GetAtomicNum()
                lines.append(f"v {atom_id} {atom_number}")

            # Add bonds as edges with bond type encoding
            for bond in mol.GetBonds():
                start_atom = bond.GetBeginAtomIdx()
                end_atom = bond.GetEndAtomIdx()
                bond_type = bond.GetBondType()

                if bond_type == Chem.rdchem.BondType.SINGLE:
                    bond_type_str = "1"
                elif bond_type == Chem.rdchem.BondType.DOUBLE:
                    bond_type_str = "2"
                elif bond_type == Chem.rdchem.BondType.TRIPLE:
                    bond_type_str = "3"
                elif bond_type == Chem.rdchem.BondType.AROMATIC:
                    bond_type_str = "4"
                else:
                    bond_type_str = "1"  # Default

                lines.append(f"e {start_atom} {end_atom} {bond_type_str}")

            output_lines.append("\n".join(lines))

    with open(output_path, 'w') as out_f:
        out_f.write("\n".join(output_lines))
