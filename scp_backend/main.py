import uvicorn
from fastapi import FastAPI, UploadFile, File
from fastapi.middleware.cors import CORSMiddleware
from pathlib import Path
from pydantic import BaseModel
from logic import run_algorithm
from fastapi.staticfiles import StaticFiles
from fastapi.responses import JSONResponse
import os

# Define upload directory
UPLOAD_DIR = Path("uploads")
UPLOAD_DIR.mkdir(parents=True, exist_ok=True)  # Create if doesn't exist

# Store temporary session data
session_data = {
    "filename": None,
    "algorithm": None,
    "parameters": {}
}

# Initialize FastAPI app
app = FastAPI()

# Mount the folder containing generated images (static serving)
app.mount("/images", StaticFiles(directory="scp_images"), name="images")

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Upload file endpoint
@app.post("/uploadfile/")
async def create_upload_file(file_upload: UploadFile = File(...)):
    data = await file_upload.read()
    save_to = UPLOAD_DIR / file_upload.filename
    with open(save_to, 'wb') as f:
        f.write(data)
    session_data["filename"] = file_upload.filename
    return {"filenames": file_upload.filename}

# Select algorithm endpoint
class AlgorithmSelection(BaseModel):
    algorithm: str

@app.post("/set-algorithm")
async def set_algorithm(selection: AlgorithmSelection):
    session_data["algorithm"] = selection.algorithm
    return {"message": "Algorithm received", "algorithm": selection.algorithm}

# Set parameters endpoint
class ParameterConfig(BaseModel):
    min_support: float
    max_overlap: float
    min_coverage: float

@app.post("/set-parameters")
async def set_parameters(params: ParameterConfig):
    session_data["parameters"] = params.dict()
    return {"message": "Parameters received", "params": params}

# Run algorithm endpoint
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

    # Call your SCP algorithm function
    result = run_algorithm(
        str(file_path),
        algorithm,
        params["min_support"],
        params["max_overlap"],
        params["min_coverage"]
    )

    return {"message": "Pipeline executed", "result": result}

@app.get("/get-images")
async def get_generated_images():
    image_folder = "scp_images"
    images = sorted([f for f in os.listdir(image_folder) if f.endswith(".png") or f.endswith(".jpg")])
    image_urls = [f"/images/{img}" for img in images[:10]]
    return JSONResponse(content={"images": image_urls})
