# -*- coding: utf-8 -*-
"""Input and working locations shared by all scripts.

LRC_DATASET  root of the input data store, laid out as <LRC_DATASET>/01_SOURCE/... (source products as downloaded) and
             <LRC_DATASET>/02_CURATED/... (the event catalogue and the 30 m rasters of the fourth national forest inventory);
             the expected sub-folders are listed in README.md, section "Input data".
LRC_WORK     working directory for intermediate files and results (q01/, s01/, v1/, v2/, cases/, figures/). The analysis
             scripts from s02 onwards read and write paths relative to this directory, so run them with it as the current
             directory (README.md, section "Running the pipeline").
"""
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATASET = os.environ.get('LRC_DATASET', os.path.join(ROOT, 'input'))
SRC = os.path.join(DATASET, '01_SOURCE')
CUR = os.path.join(DATASET, '02_CURATED')
WORK = os.environ.get('LRC_WORK', os.path.join(ROOT, 'work'))

# project 30 m grid: TWD97 / TM2 zone 121 (EPSG:3826), upper-left corner (X0, Y1), W columns x H rows; cell = row * W + col
X0, Y1, RES, W, H = 151050.0, 2798910.0, 30.0, 6680, 12558
