import uvicorn
from fastapi import FastAPI, UploadFile, File, Form
from fastapi.middleware.cors import CORSMiddleware
from pathlib import Path
from fastapi.responses import JSONResponse
from tusm import TUSM  # Import your TUSM class
from vis import draw_subgraph
import os
import base64
import io
import matplotlib.pyplot as plt

# Define upload directory
UPLOAD_DIR = Path("uploads1")
UPLOAD_DIR.mkdir(parents=True, exist_ok=True)

# Initialize FastAPI app
app = FastAPI()

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
@app.post("/run-tusm/")
async def run_tusm_endpoint(
    file_upload: UploadFile = File(...),
    k: int = Form(...),
    epsilon: float = Form(...),
    delta: float = Form(...)
):
    # Save uploaded file
    file_path = UPLOAD_DIR / file_upload.filename
    data = await file_upload.read()
    with open(file_path, "wb") as f:
        f.write(data)

    # --- CHECKPOINTS ---
    print("=== CHECKPOINT: File Saved ===")
    print(f"Absolute path where file was written: {file_path.resolve()}")
    print(f"Uploaded filename: {file_upload.filename}")
    print(f"Parameters received -> k: {k}, epsilon: {epsilon}, delta: {delta}")
    print("==============================")

    # Run TUSM algorithm
    tusm_instance = TUSM(str(file_path))

    # --- CHECKPOINT ---
    print(f"Running TUSM on file: {file_path.resolve()}")

    min_queue, _ = tusm_instance.mine(k=k, eps=epsilon, delta=delta)

    # Extract DFS codes from min_queue
    top_k_subgraphs = []
    images_base64 = []
    idx = 1
    # while not min_queue.is_empty():
    #     _, _, _, dfs_code = min_queue.pop()

    #     # Use the edges stored in DFSCode
    #     edges = dfs_code.edges  
    #     top_k_subgraphs.append(str(dfs_code))

    #     # Generate image
    #     fig, ax = plt.subplots()
    #     draw_subgraph(ax, edges, title=f"Subgraph {idx}")

    #     # Save to buffer
    #     buf = io.BytesIO()
    #     fig.savefig(buf, format='png', bbox_inches='tight')
    #     buf.seek(0)
    #     plt.close(fig)

    #     # Convert to base64
    #     img_b64 = base64.b64encode(buf.read()).decode("utf-8")
    #     images_base64.append(img_b64)
    #     idx += 1


    # --- CHECKPOINT ---
    print("=== CHECKPOINT: TUSM run completed ===")
    print(f"Number of top-k subgraphs: {len(top_k_subgraphs)}")

    return JSONResponse(content={
        "filename": file_upload.filename,
        "k": k,
        "epsilon": epsilon,
        "delta": delta,
        "top_k_subgraphs": top_k_subgraphs[::-1],
        "images": images_base64[::-1]
    })

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)
