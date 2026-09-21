[README.md](https://github.com/user-attachments/files/32475654/README.md)
# 🧬 lncRank-Spec: Multi-Agent lncRNA Target Discovery & Specificity Scoring Suite

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![FastAPI](https://img.shields.io/badge/FastAPI-0.100%2B-009688.svg)](https://fastapi.tiangolo.com/)
[![Streamlit](https://img.shields.io/badge/Streamlit-1.25%2B-FF4B4B.svg)](https://streamlit.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**lncRank-Spec** is an open-source computational biology platform, data engine, and multi-agent AI framework designed to mine, process, and rank **long non-coding RNAs (lncRNAs)** based on tissue-, cell-type-, and disease-specificity. 

While traditional gene therapies target protein-coding genes (mRNAs) that are frequently expressed ubiquitously across normal organs—often leading to severe off-target toxicity—lncRNAs exhibit ultra-restricted, context-dependent expression profiles. **lncRank-Spec** leverages multi-omic transcriptomics (GTEx, TCGA, GEO, DepMap, lncATLAS, FANTOM6, and RNAcentral) to systematically isolate "zero-toxicity" therapeutic targets and diagnostic liquid biopsy biomarkers.

---

## 🌟 Key Features

* **Dual-Filter Specificity & Actionability Engine:** Applies dynamic detectability thresholds ($\text{TPM} > 1.5$) to eliminate background transcriptional noise, paired with Yanai's **Tau ($\tau$) index** to penalize "leaky" transcripts.
* **The lncScore Mathematical Motor:** Computes target priority using the quadratic specificity formula:
  $$\text{lncScore} = \tau^2 \times \log_2(\text{Mean Expression} + 1)$$
* **Multi-Agent AI Architecture:** Inspired by the *Virtual Biotech* framework (*Science*, 2026), orchestrating specialized AI divisions for Data Ingestion, Math Scoring, Spatial Subcellular Localization, Contextual Vulnerability, and Cross-Species Translation.
* **Unified Cross-Species Validation:** Evaluates orthology and expression profiles across **Humans, Mice, and Rats** under a single mathematical model to de-risk pre-clinical animal studies.
* **Interactive Streamlit Web Dashboard:** No-code graphical interface for real-time target exploration, multi-organ safety profiling, and custom data ingestion.
* **High-Performance FastAPI REST Microservice:** Programmatic endpoints enabling external bioinformatic workflows and agentic tool-calling via the Model Context Protocol (MCP).
* **Live Public API Connectors:** Built-in integration layer for **cBioPortal (TCGA)**, **GTEx Portal REST API (v2)**, and **DepMap Portal (26Q1)**.

---

## 📁 Repository Structure

```text
lncRank-Spec/
├── prescreening-dashboard.py   # Streamlit interactive web application
├── fastapi_multiagent_api.py  # FastAPI RESTful microservice & MCP endpoints
├── multi_agent_prescreening.py # Core multi-agent pipeline engine & CLI
├── requirements.txt            # Python dependencies
├── run_dashboard.sh            # One-click startup script (macOS/Linux/WSL)
├── docker-compose.yml          # Container orchestration configuration
└── README.md                   # Repository documentation
```

---

## 🚀 Quick Start Guide

### Prerequisites
* Python 3.10 or higher
* `pip` package manager
* (Optional) Docker & Docker Compose

### 1. Installation

Clone the repository and install dependencies:
```bash
git clone https://github.com/your-username/lncRank-Spec.git
cd lncRank-Spec
pip install -r requirements.txt
```

### 2. Running Locally

#### Option A: One-Click Startup Script (macOS / Linux / WSL)
```bash
chmod +x run_dashboard.sh
./run_dashboard.sh
```

#### Option B: Manual Terminal Execution
Launch the FastAPI backend in one terminal window:
```bash
uvicorn fastapi_multiagent_api:app --reload --port 8000
```
*(Interactive API documentation available at `http://localhost:8000/docs`)*

In a second terminal window, launch the Streamlit dashboard:
```bash
streamlit run prescreening-dashboard.py
```
*(Dashboard opens automatically at `http://localhost:8501`)*

### 3. Containerized Setup (Docker Compose)
To run both the backend API and frontend Streamlit dashboard in isolated Docker containers:
```bash
docker-compose up --build
```

---

## 🧬 Multi-Agent AI Pipeline Workflow

```text
+-----------------------------------------------------------------------------------+
|                         EXECUTIVE ORCHESTRATOR (CEO AGENT)                         |
+-----------------------------------------------------------------------------------+
                                          |
     +------------------------------------+------------------------------------+
     |                                    |                                    |
+----+-------------------+      +---------+----------+      +------------------+----+
|  1. DATA INGESTION     |      |  2. lncSCORE       |      |  3. SPATIAL      |
|     DIVISION           | ---> |     MATH ENGINE    | ---> |     LOCALIZATION |
|  - GTEx, TCGA, GEO     |      |  - Tau Index (τ)   |      |  - lncATLAS RCI  |
|  - Noise (TPM > 1.5)   |      |  - lncScore Formula|      |  - ASO vs siRNA  |
+------------------------+      +--------------------+      +------------------+----+
                                                                               |
     +------------------------------------+------------------------------------+
     |                                    |
+----+-------------------+      +---------+----------+
|  4. DEPENDENCY & SAFETY| ---> |  5. CROSS-SPECIES  |
|  - DepMap Chronos      |      |     TRANSLATION    |
|  - Multi-Organ Safety  |      |  - Human/Mouse/Rat |
+------------------------+      +--------------------+
```

1. **Ingestion & Noise Filtering Agent:** Harvests expression matrices and enforces $TPM > 1.5$.
2. **lncScore Math Engine Agent:** Calculates Yanai's $\tau$ index and computes weighted lncScores.
3. **Spatial & Localization Agent:** Evaluates lncATLAS Relative Concentration Index (RCI) to route candidates to **Antisense Oligonucleotides (ASOs)** for nuclear targets ($RCI \ge 1.0$) or **siRNA/miRNA sponges** for cytoplasmic targets ($RCI \le -1.0$).
4. **Dependency & Safety Agent:** Cross-references DepMap CRISPR Chronos scores ($\le -0.6$) and verifies zero expression across healthy organ baselines.
5. **Cross-Species Translation Agent:** Aligns mouse/rat ortholog Tau scores to validate animal models for pre-clinical trials.

---

## 📡 API Endpoints Reference

The FastAPI service exposes the following RESTful endpoints:

* `GET /`: API health check.
* `GET /api/v1/integrations/status`: Real-time status monitor for live TCGA, GTEx, and DepMap connectors.
* `POST /api/v1/agents/ingestion`: Data ingestion & noise filtering.
* `POST /api/v1/agents/math-engine`: Specificity scoring ($\tau$) and lncScore calculation.
* `POST /api/v1/agents/spatial`: Subcellular localization & drug modality routing.
* `POST /api/v1/agents/vulnerability`: DepMap dependency & multi-organ safety profiling.
* `POST /api/v1/agents/orthology`: Cross-species orthology alignment.
* `POST /api/v1/orchestrator/prescreen`: Full end-to-end multi-agent target pre-screening.

---

## 🌐 Deploying to Streamlit Community Cloud

1. Push this repository to GitHub.
2. Visit [share.streamlit.io](https://share.streamlit.io/) and log in with GitHub.
3. Click **"New app"**.
4. Select repository `lncRank-Spec`, branch `main`, and set main file path to `prescreening-dashboard.py`.
5. Click **Deploy**!

---

## 📜 Scientific Foundations & References

1. **Cell-Type Specificity in Drug Success:** Zhang et al., *The Virtual Biotech: A multi-agent AI framework for therapeutic discovery and development*, **Science** (2026).
2. **Tissue Specificity Benchmark:** Kryuchkova-Mostacci N. & Robinson-Rechavi M., *A benchmark of gene expression tissue-specificity metrics*, **Briefings in Bioinformatics** (2017).
3. **Cancer Dependency Map:** Broad Institute DepMap Portal (26Q1 Release).
4. **Subcellular Localization:** Mas-Ponte D. et al., *LncATLAS database for subcellular localization of long noncoding RNAs*, **RNA** (2017).

---

## 📄 License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.
