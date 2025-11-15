# 🧬 BioMiner  
*A Subgraph Coverage Pattern Generator for Chemical Compound Analysis & Drug Discovery*

BioMiner is a powerful tool designed to analyze chemical compound structures and generate **Subgraph Coverage Patterns (SCPs)**. These patterns help researchers identify key molecular features, understand structure–property relationships, and accelerate **drug discovery** workflows.

---

## 🚀 Features

- 🔍 **Extract meaningful subgraphs** from chemical structures  
- 🧠 **Generates Subgraph Coverage Patterns (SCPs)** for compound analysis  
- ⚗️ **Designed for cheminformatics & drug discovery research**  
- 🧬 Uses **gSpan / TKG** subgraph mining algorithms  
- 📊 Supports graph statistics (density, clustering, label distributions)  
- ⚡ Fast, scalable backend built with **FastAPI**  
- 💻 Modern, responsive **React** UI  
- 🧩 Integrates RDKit and NetworkX for molecular graph processing  

---

## 🧱 Tech Stack

### **Backend**
- Python  
- FastAPI  
- gSpan
- TKG algorithm implementation

### **Frontend**
- React JS
- Responsive UI components

---

## 📁 How to Run:
- frontend:
    - npm run dev
    - npm install
- Backend:
    - uvicorn main:app --reload --port 8000 for scp algorithm
    - uvicorn main:app --reload --port 8001 for tkg algorithm





