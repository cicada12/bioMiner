import uvicorn
from fastapi import FastAPI, UploadFile, File, Form
from fastapi.middleware.cors import CORSMiddleware
from pathlib import Path
from pydantic import BaseModel
from logic import run_algorithm
from fastapi.staticfiles import StaticFiles
from fastapi.responses import JSONResponse
import os
from rdkit import Chem
import networkx as nx
import pandas as pd

# ---------------------- DIRECTORIES --------------------------
UPLOAD_DIR = Path("uploads")                # contains dummy datasets + user files
USER_UPLOADS_DIR = Path("uploaded_files")   # user-uploaded datasets
UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
USER_UPLOADS_DIR.mkdir(parents=True, exist_ok=True)

# Store temporary session data for SCP algorithm
session_data = {
    "filename": None,
    "algorithm": None,
    "parameters": {}
}

# ---------------------- FASTAPI SETUP ------------------------
app = FastAPI()

# serve generated scp images
app.mount("/images", StaticFiles(directory="scp_images"), name="images")

# CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ------------------------- 1. UPLOAD FILE -------------------------
@app.post("/uploadfile/")
async def create_upload_file(file_upload: UploadFile = File(...)):
    data = await file_upload.read()
    save_to = UPLOAD_DIR / file_upload.filename
    with open(save_to, "wb") as f:
        f.write(data)

    # also store in uploaded_files/ for consistency
    save_to_user = USER_UPLOADS_DIR / file_upload.filename
    with open(save_to_user, "wb") as f:
        f.write(data)

    session_data["filename"] = file_upload.filename
    return {"filenames": file_upload.filename}


# --------------------- 2. SELECT ALGORITHM ------------------------
class AlgorithmSelection(BaseModel):
    algorithm: str

@app.post("/set-algorithm")
async def set_algorithm(selection: AlgorithmSelection):
    session_data["algorithm"] = selection.algorithm
    return {"message": "Algorithm received", "algorithm": selection.algorithm}


# -------------------- 3. SET PARAMETERS ---------------------------
class ParameterConfig(BaseModel):
    min_support: float
    max_overlap: float
    min_coverage: float

@app.post("/set-parameters")
async def set_parameters(params: ParameterConfig):
    session_data["parameters"] = params.dict()
    return {"message": "Parameters received", "params": params}


# -------------------- 4. RUN SCP ALGORITHM ------------------------
@app.post("/run")
async def run_algorithm_endpoint():
    filename = session_data["filename"]
    algorithm = session_data["algorithm"]
    params = session_data["parameters"]

    if not filename or not algorithm or not params:
        return {"error": "Missing data. Ensure file, algorithm, and parameters are set."}

    file_path = UPLOAD_DIR / filename
    if not file_path.exists():
        return {"error": "File not found"}

    # Run SCP algorithm
    result = run_algorithm(
        str(file_path),
        algorithm,
        params["min_support"],
        params["max_overlap"],
        params["min_coverage"]
    )

    return {"message": "Pipeline executed", "result": result}


# ---------------------- 5. GET OUTPUT IMAGES ------------------------
@app.get("/get-images")
async def get_generated_images():
    image_folder = "scp_images"
    if not os.path.exists(image_folder):
        return {"images": []}

    images = sorted([
        f for f in os.listdir(image_folder)
        if f.endswith(".png") or f.endswith(".jpg")
    ])

    image_urls = [f"/images/{img}" for img in images[:10]]
    return JSONResponse(content={"images": image_urls})


# =====================================================================
# =====================   SMILES  STATS  API   ========================
# =====================================================================

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


@app.post("/compute_stats/")
async def compute_stats(
    filename: str = Form(None),
    use_dummy: bool = Form(False),
    dummy_name: str = Form(None)
):
    """
    Computes descriptive SMILES graph statistics from either:
    - dummy dataset (uploads/)
    - uploaded dataset (uploaded_files/)
    """

    # ------------------ SELECT WHICH DATASET ------------------
    if use_dummy and dummy_name:
        file_path = UPLOAD_DIR / dummy_name
        if not file_path.exists():
            return {"error": f"Dummy dataset {dummy_name} not found."}

    elif filename:
        file_path = USER_UPLOADS_DIR / filename
        if not file_path.exists():
            return {"error": f"Uploaded file {filename} not found."}

    else:
        return {"error": "No dataset selected."}

    # Read SMILES from the file
    try:
        with open(file_path, "r") as f:
            smiles_list = [line.strip() for line in f if line.strip()]
    except:
        return {"error": "Could not read file."}

    total = len(smiles_list)
    invalid = 0
    stats_rows = []

    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)

        if mol is None:
            invalid += 1
            continue

        V = mol.GetNumAtoms()
        E = mol.GetNumBonds()

        G = mol_to_nx(mol)
        clustering = nx.average_clustering(G) if V > 0 else 0.0
        density = graph_density(V, E)

        stats_rows.append({
            "vertices": V,
            "edges": E,
            "density": density,
            "clustering": clustering,
        })

    if not stats_rows:
        return {"error": "No valid SMILES detected in dataset."}

    df = pd.DataFrame(stats_rows)

    # Summary statistics
    summary = {
        # " total_smiles ": total,
        # " valid_smiles ": int(len(df)),
        # " invalid_smiles ": int(invalid),

        "avg_vertices": float(df["vertices"].mean()),
        "avg_edges": float(df["edges"].mean()),
        "avg_density": float(df["density"].mean()),
        "avg_clustering": float(df["clustering"].mean()),

        "min_vertices": int(df["vertices"].min()),
        "max_vertices": int(df["vertices"].max()),

        "min_edges": int(df["edges"].min()),
        "max_edges": int(df["edges"].max()),
    }

    return summary
