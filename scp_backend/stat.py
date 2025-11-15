from rdkit import Chem
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

file_path = "/home/nida/Desktop/BTP/btp2/bioMiner/scp_backend/uploads/graph_transactional_dataset2.txt"  # change to your file

def graph_density(V, E):
    if V <= 1:
        return 0.0
    return (2 * E) / (V * (V - 1))

def mol_to_nx(mol):
    G = nx.Graph()
    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx())
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    return G

smiles = []
with open(file_path, "r") as f:
    smiles = [line.strip() for line in f if line.strip()]

stats = []
invalid = 0

for smi in smiles:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        invalid += 1
        continue

    V = mol.GetNumAtoms()
    E = mol.GetNumBonds()

    G = mol_to_nx(mol)

    # clustering coefficient (average over nodes)
    clustering = nx.average_clustering(G)

    density = graph_density(V, E)

    stats.append({
        "smiles": smi,
        "vertices": V,
        "edges": E,
        "density": density,
        "clustering_coeff": clustering
    })

df = pd.DataFrame(stats)

# Save the data
df.to_csv("graph_stats.csv", index=False)
print("Saved detailed stats → graph_stats.csv")

# Print summary
print("\n===== GRAPH STATISTICS =====")
print(f"Total SMILES: {len(smiles)}")
print(f"Valid molecules: {len(df)}")
print(f"Invalid molecules: {invalid}")

print("\nAverage vertices:", df["vertices"].mean())
print("Average edges:", df["edges"].mean())
print("Average density:", df["density"].mean())
print("Average clustering coefficient:", df["clustering_coeff"].mean())

# # ------------ HISTOGRAMS ------------ #

# def plot_hist(column, title, filename, bins=30):
#     plt.figure(figsize=(7,5))
#     plt.hist(df[column], bins=bins)
#     plt.title(title)
#     plt.xlabel(column)
#     plt.ylabel("Frequency")
#     plt.grid(True, alpha=0.3)
#     plt.tight_layout()
#     plt.savefig(filename)
#     plt.close()

# plot_hist("vertices", "Histogram of Number of Vertices (Atoms)", "hist_vertices.png")
# plot_hist("edges", "Histogram of Number of Edges (Bonds)", "hist_edges.png")
# plot_hist("density", "Histogram of Graph Density", "hist_density.png")
# plot_hist("clustering_coeff", "Histogram of Clustering Coefficient", "hist_clustering.png")

# print("\nSaved histograms:")
# print(" - hist_vertices.png")
# print(" - hist_edges.png")
# print(" - hist_density.png")
# print(" - hist_clustering.png")
