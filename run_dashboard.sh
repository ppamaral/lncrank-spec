#!/bin/bash
# Startup script for lncRank-Spec Streamlit Dashboard and FastAPI backend

echo "=========================================================="
echo "Setting up lncRank-Spec Multi-Agent Platform Environment"
echo "=========================================================="

# 1. Check Python version
python3 --version

# 2. Install dependencies
pip install -r requirements.txt

# 3. Launch FastAPI backend in background
echo "Starting FastAPI Backend on http://localhost:8000 ..."
python3 -m uvicorn fastapi_multiagent_api:app --host 0.0.0.0 --port 8000 &

# 4. Launch Streamlit dashboard
echo "Starting Streamlit Dashboard on http://localhost:8501 ..."
python3 -m streamlit run prescreening-dashboard.py --server.port 8501 --server.address 0.0.0.0
