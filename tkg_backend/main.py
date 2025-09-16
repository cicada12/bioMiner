import uvicorn
from fastapi import FastAPI, UploadFile, File, Form
from fastapi.middleware.cors import CORSMiddleware
from pathlib import Path
from fastapi.responses import JSONResponse
from tusm import TUSM
import os
import base64
from plot import plot_subgraphs_from_txt  # your plotting function
import glob

# Define upload directory
UPLOAD_DIR = Path("uploads1")
UPLOAD_DIR.mkdir(parents=True, exist_ok=True)

# Directory to save images
IMAGE_DIR = Path("images")
IMAGE_DIR.mkdir(parents=True, exist_ok=True)

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

    print("=== CHECKPOINT: File Saved ===")
    print(f"Absolute path where file was written: {file_path.resolve()}")
    print(f"Parameters received -> k: {k}, epsilon: {epsilon}, delta: {delta}")
    print("==============================")

    # Run TUSM algorithm
    tusm_instance = TUSM(str(file_path))
    print(f"Running TUSM on file: {file_path.resolve()}")
    min_queue, _ = tusm_instance.mine(k=k, eps=epsilon, delta=delta)

    # Save results into a text file
    BASE_DIR = Path(__file__).resolve().parent
    output_file = BASE_DIR / "mined_subgraphs.txt"

    top_k_subgraphs = []
    idx = 1

    with open(output_file, "w") as f:
        while not min_queue.is_empty():
            _, _, _, dfs_code = min_queue.pop()
            subgraph_str = f"Subgraph {idx}: {str(dfs_code)}\n"
            f.write(subgraph_str + "\n")
            top_k_subgraphs.append(str(dfs_code))
            idx += 1

    print(f"Saved {len(top_k_subgraphs)} subgraphs to {output_file}")

    # --- NEW: Generate images ---
    # Clear old images first
    for img_file in IMAGE_DIR.glob("*.png"):
        img_file.unlink()

    # Call plotting function (it will save images in IMAGE_DIR)
    plot_subgraphs_from_txt(str(output_file), str(IMAGE_DIR))

    # Encode images as base64 to send to frontend
    image_files = sorted(IMAGE_DIR.glob("*.png"))  # sort by name
    image_base64_list = []

    for img_path in image_files:
        with open(img_path, "rb") as img_f:
            encoded = base64.b64encode(img_f.read()).decode("utf-8")
            image_base64_list.append({
                "filename": img_path.name,
                "data": encoded
            })

    print("=== CHECKPOINT: Images generated and encoded ===")

    return JSONResponse(content={
        "filename": file_upload.filename,
        "k": k,
        "epsilon": epsilon,
        "delta": delta,
        "top_k_subgraphs": top_k_subgraphs[::-1],
        "images": image_base64_list
    })

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8001)
