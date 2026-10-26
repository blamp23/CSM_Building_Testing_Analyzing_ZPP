"""
Chapter pipeline path config — single source of truth for the 9 chapter scripts.

Repo-internal paths (functions/, data/, chapter_assets/) resolve from __file__.
External large-data paths (model .mat files, pooled hopsy samples, hpf-matched
exchange envelopes) live outside git. By default they resolve to the OneDrive
working location on the author's laptop; override with the CSM_DATA_ROOT
environment variable to point at a different mirror.

    export CSM_DATA_ROOT=/path/to/Discrete_Models

Required external artifacts (produced by pipeline/00_bake, 01_extract, 02_sampling,
03_loop_removal, 04_matched_floor and the sampling-side pooling step):
    $CSM_DATA_ROOT/reviewer_packet_v4_slim/v7_models/trans_rfastcormics_*.mat
    $CSM_DATA_ROOT/reviewer_packet_v4_slim/s130_v7_ledger.csv
    $CSM_DATA_ROOT/grace_v4/results/pooled_v7/{cell}_pooled.csv
    $CSM_DATA_ROOT/grace_v4/results/subsets_v7/{cell}_subset_membership.csv
"""
from pathlib import Path
import os

# ----- repo-internal (committed) -----
# This file lives at: pipeline/05_chapter/scripts/config.py
CSM_ROOT = Path(__file__).resolve().parents[3]
CHAPTER_HOME = Path(__file__).resolve().parents[1]
FUNCTIONS = CSM_ROOT / "functions"
DATA = CSM_ROOT / "data"
CHAPTER_ASSETS = CHAPTER_HOME / "chapter_assets"
COMPARISON_TABLES = CHAPTER_ASSETS / "comparison_tables"
TABLES = CHAPTER_ASSETS / "tables"
FIGURES = CHAPTER_ASSETS / "figures"

# ----- external (large binary artifacts; not in git) -----
_DEFAULT_EXTERNAL = Path(
    "/Users/lamp_b/Library/CloudStorage/OneDrive-TexasA&MUniversity/"
    "Hala, David's files - Benji_COBRA/Tanguay_Data/Discrete_Models"
)
DATA_ROOT = Path(os.environ.get("CSM_DATA_ROOT", _DEFAULT_EXTERNAL))

V7_MODELS = DATA_ROOT / "reviewer_packet_v4_slim" / "v7_models"
BASE_MODELS = DATA_ROOT / "ocr_anchored_extraction" / "models"
LEDGER = DATA_ROOT / "reviewer_packet_v4_slim" / "s130_v7_ledger.csv"
POOLED = DATA_ROOT / "grace_v4" / "results" / "pooled_v7"
SUBSETS = DATA_ROOT / "grace_v4" / "results" / "subsets_v7"
SUSCEPTIBILITY = DATA_ROOT / "grace_v4" / "results" / "v7_susceptibility"
ATP_SWEEP = DATA_ROOT / "grace_v4" / "results" / "atp_sweep_v7"
EXCHANGE_BOUNDS = DATA / "phase2_exchange_bounds_long.csv"

# ----- canonical cell ordering -----
CONDS = ("BL", "D", "LD")
HPFS = (24, 48, 72, 96, 120)
CELLS = tuple(f"{c}_{h}" for c in CONDS for h in HPFS)
