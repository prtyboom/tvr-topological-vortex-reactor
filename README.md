# Topological Vortex Reactor (TVR) v7.1

DOI (Zenodo): https://doi.org/10.5281/zenodo.17797496

This repository contains a purely quantum mechanical model of the Topological Vortex Reactor (TVR). The core is a nonlinear parametric oscillator coupled to an electromagnetic helicity operator, plus a tunneling degree of freedom for studying Floquet-assisted modulation of tunneling rates.

## Files

- `tvr_v7_1.tex` — LaTeX source of the article.  
- `tvr_v7_1.pdf` — compiled PDF (not stored in git, but can be built locally).  
- `requirements.txt` — minimal Python dependencies.  
- `python/double_well_static.py` — static double-well Hamiltonian example.  
- `python/double_well_floquet.py` — Floquet-driven double-well tunneling example.

## Requirements

Python 3.12 (or later) is recommended.

Minimal Python packages:

- numpy  
- scipy  
- matplotlib  

Install (in a terminal / cmd):

    python -m pip install --upgrade pip
    python -m pip install -r requirements.txt

MiKTeX or another LaTeX distribution is required to build the PDF.

## Build the PDF

From the repository directory:

    pdflatex tvr_v7_1.tex
    pdflatex tvr_v7_1.tex

This produces `tvr_v7_1.pdf`.

## Run the static double-well example

From the repository directory:

    python python\double_well_static.py

Expected output (dimensionless energies for default parameters):

    Lowest energy levels (dimensionless):
      E[0] = 0.373779
      E[1] = 0.379814
      E[2] = 0.941599
      E[3] = 1.084628

A matplotlib window should show the double-well potential and the two lowest eigenfunctions.

## Run the Floquet-driven double-well example

From the repository directory:

    python python\double_well_floquet.py

Default parameters inside the script (as of v7.1):

- epsilon = 0.2  
- mu = 0.3  
- N = 1024  
- n_periods = 10  
- steps_per_period = 200  

Example of console output:

    Simulation finished: epsilon=0.2, mu=0.3
    Total time in tau units: 62.832
    Final probabilities: P_L=0.7639, P_R=0.2259

A matplotlib window shows the probabilities in the left (y < 0) and right (y > 0) halves as functions of the dimensionless time.

## Links

GitHub repository:  
https://github.com/prtyboom/tvr-topological-vortex-reactor

Zenodo archive (DOI):  
https://doi.org/10.5281/zenodo.17797496