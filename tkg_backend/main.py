import uvicorn
from fastapi import FastAPI, UploadFile, File, Form
from fastapi.middleware.cors import CORSMiddleware
from pathlib import Path
from fastapi.responses import JSONResponse
from tusm import TUSM
import base64
from plot import plot_subgraphs_from_txt
import networkx as nx
import numpy as np
from collections import Counter

# ============================================================
#  GRAPH PARSER + STATISTICS
# ============================================================

def load_graphs_from_txt(path):
    graphs = []
    current = None
    
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("t"):  # new graph starts
                if current:
                    graphs.append(current)
                current = {"nodes": {}, "edges": []}

            elif line.startswith("v"):
                _, vid, lab = line.split()
                current["nodes"][int(vid)] = int(lab)

            elif line.startswith("e"):
                _, a, b, lab, w = line.split()
                current["edges"].append((int(a), int(b), int(lab), float(w)))

    if current:
        graphs.append(current)

    return graphs


def compute_graph_stats(graphs):
    num_graphs = len(graphs)

    node_counts = []
    edge_counts = []
    label_distribution = Counter()
    densities = []
    clusterings = []

    for g in graphs:
        node_counts.append(len(g["nodes"]))
        edge_counts.append(len(g["edges"]))
        label_distribution.update(g["nodes"].values())

        G = nx.Graph()
        for n, lab in g["nodes"].items():
            G.add_node(n, label=lab)
        for a, b, lab, w in g["edges"]:
            G.add_edge(a, b, label=lab, weight=w)

        densities.append(nx.density(G))
        clusterings.append(nx.average_clustering(G) if len(G.nodes()) > 1 else 0)

    stats = {
        "total_graphs": num_graphs,
        "avg_nodes": float(np.mean(node_counts)),
        "avg_edges": float(np.mean(edge_counts)),
        "min_graph_size": min(node_counts),
        "max_graph_size": max(node_counts),
        "label_distribution": dict(label_distribution),
        "avg_density": float(np.mean(densities)),
        "avg_clustering": float(np.mean(clusterings)),
    }

    return stats


# ============================================================
#  FASTAPI SETUP
# ============================================================

UPLOAD_DIR = Path("uploads1")
UPLOAD_DIR.mkdir(parents=True, exist_ok=True)

IMAGE_DIR = Path("images")
IMAGE_DIR.mkdir(parents=True, exist_ok=True)

app = FastAPI()

# ------------------------------------------------------------
# CORS MUST COME RIGHT AFTER APP CREATION
# ------------------------------------------------------------
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ============================================================
#  ⭐ File Upload Endpoint (Required for Frontend)
# ============================================================

@app.post("/uploadfile/")
async def upload_file(file_upload: UploadFile = File(...)):
    file_path = UPLOAD_DIR / file_upload.filename
    data = await file_upload.read()

    with open(file_path, "wb") as f:
        f.write(data)

    return {"filename": file_upload.filename}


# ============================================================
#  ⭐ Compute Graph Statistics Only (No TUSM)
# ============================================================

@app.post("/compute_graph_stats/")
async def compute_graph_stats_endpoint(
    filename: str = Form(None),
    use_dummy: bool = Form(False),
    dummy_name: str = Form(None)
):
    # choose file source
    if use_dummy and dummy_name:
        file_path = UPLOAD_DIR / dummy_name
    elif filename:
        file_path = UPLOAD_DIR / filename
    else:
        return {"error": "No dataset selected."}

    if not file_path.exists():
        return {"error": "File not found."}

    graphs = load_graphs_from_txt(str(file_path))
    stats = compute_graph_stats(graphs)

    return JSONResponse(content=stats)


# ============================================================
#  ⭐ TUSM MINING + STATS + IMAGE GENERATION
# ============================================================

@app.post("/run-tusm/")
async def run_tusm_endpoint(
    file_upload: UploadFile = File(...),
    k: int = Form(...),
    epsilon: float = Form(...),
    delta: float = Form(...)
):
    # save uploaded graph dataset
    file_path = UPLOAD_DIR / file_upload.filename
    data = await file_upload.read()
    with open(file_path, "wb") as f:
        f.write(data)

    print("=== CHECKPOINT: File Saved ===")

    # compute stats
    graphs = load_graphs_from_txt(str(file_path))
    graph_stats = compute_graph_stats(graphs)

    print("=== Graph Stats Computed ===")
    print(graph_stats)

    # run TUSM
    tusm_instance = TUSM(str(file_path))
    min_queue, _ = tusm_instance.mine(k=k, eps=epsilon, delta=delta)

    output_file = Path(__file__).resolve().parent / "mined_subgraphs.txt"

    top_k_subgraphs = []
    idx = 1

    with open(output_file, "w") as f:
        while not min_queue.is_empty():
            _, _, _, dfs_code = min_queue.pop()
            f.write(f"Subgraph {idx}: {str(dfs_code)}\n\n")
            top_k_subgraphs.append(str(dfs_code))
            idx += 1

    # delete old images
    for img_file in IMAGE_DIR.glob("*.png"):
        img_file.unlink()

    # create new images
    plot_subgraphs_from_txt(str(output_file), str(IMAGE_DIR))

    # convert images to base64
    image_files = sorted(IMAGE_DIR.glob("*.png"))
    image_base64_list = []

    for img_path in image_files:
        with open(img_path, "rb") as img_f:
            encoded = base64.b64encode(img_f.read()).decode("utf-8")
            image_base64_list.append({
                "filename": img_path.name,
                "data": encoded
            })

    return JSONResponse(content={
        "filename": file_upload.filename,
        "k": k,
        "epsilon": epsilon,
        "delta": delta,
        "graph_stats": graph_stats,
        "top_k_subgraphs": top_k_subgraphs[::-1],
        "images": image_base64_list
    })


# ============================================================

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8001)
