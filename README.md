# Reproducibility package for article: Spin structures and phase diagrams of the spin-5/2 triangular-lattice antiferromagnet Na₂BaMn(PO₄)₂ under magnetic field

[![License: CC BY-4.0](https://img.shields.io/badge/License-CC%20BY--4.0-lightgrey.svg)](LICENSE)
[![Launch Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/me2d09/nbmpo-diffraction-paper/HEAD?labpath=index.ipynb)

This repository contains all Python code and raw data needed to **reproduce every
figure** in the manuscript

> “Spin structures and phase diagrams of the spin-5/2 triangular-lattice
> antiferromagnet Na₂BaMn(PO₄)₂ under magnetic field”.

The core of the workflow is the Jupyter notebook  
`index.ipynb` — running it from top to bottom regenerates the publication-ready
PDF versions of Figures 2–8 in the main text and Figure S3 in the Supplement.

## Repository layout
```
.
├── index.ipynb          # one-click pipeline that orchestrates everything
├── 2dmodel/             # Monte-Carlo outputs for the 2D easy-axis triangular Heisenberg model (no interlayer J)
├── 2dmodel/             # Monte-Carlo outputs for the 3D model with frustrated interlayer couplings (Ja, J2)
├── rawdata/             # raw neutron diffraction files from D23, these data could e downloaded again from ILL data repository
├── rawhc/               # raw heat-capacity files, including puck calibration, measured at MGML facility
├── rawsim/              # calculated magnetization curves
├── requirements.txt     # python requirements environment freezed for python 3.10 for full reproducibility
└── LICENSE              # CC-BY-4.0
````

### Data provenance

* **`rawdata/`** — Original neutron diffraction files recorded on the ILL D23 diffractometer. These are unmodified exports from the instrument and can be re-downloaded from the ILL data repository. Experiment identifiers are referenced in the python script, together with the code to downlaod them. ILL username and password is needed. This step is skippeb by defualt since the data are mirrored in this repository. The notebooks read this files **read-only**.

* **`rawhc/`** — Original heat-capacity data measured on PPMS at the MGML facility, including puck calibration runs. Any addenda subtraction/processing happens in the notebooks; the files here remain untouched.

* **`rawsim/`** — Calculated magnetization curves used in the figures/analysis. These are model outputs saved to disk so the full pipeline can be reproduced without rerunning long simulations. Used for Fig. 6.

* **`2dmodel/`** — Monte-Carlo outputs for the **2D easy-axis triangular Heisenberg model** (nearest-neighbour intralayer exchange only; **no interlayer coupling**). Contains the simulation products consumed by the plotting notebook.

* **`3dmodel/`** — Monte-Carlo outputs for the **3D model with frustrated interlayer couplings** (competing $J_a$, $J_2$ on top of intralayer exchange). Also consumed directly by the plotting notebook.

* **`index.ipynb`** — Orchestrates the full, one-click rebuild of all figures from the above inputs (no manual edits to the raw folders). In addition, raw trusted neutron data([this](https://doi.ill.fr/10.5291/ILL-DATA.5-41-1252) and [this](https://doi.ill.fr/10.5291/ILL-DATA.CRG-3064)) can be downloaded from the [https://data.ill.eu/] portal when the 3 years embargo period will be over. Ehis step is not needed as the copy of the needed data are included in rawdata folder.

* **`requirements.txt`** — Frozen Python 3.10 environment to reproduce the exact software stack used to generate the figures.

* **`LICENSE`** — CC-BY-4.0; please cite the manuscript if you reuse data or figures.


---

## Quick start

### A) Local setup (Python 3.10 + requirements.txt)

```bash
# 1) Clone
git clone https://github.com/me2d09/nbmpo-diffraction-paper.git
cd nbmpo-diffraction-paper

# 2) Create a Python 3.10 environment

# Option A: venv (Linux/macOS)
python3.10 -m venv .venv
source .venv/bin/activate

# Option A (Windows PowerShell)
# py -3.10 -m venv .venv
# .\.venv\Scripts\Activate.ps1

# Option B: conda
# conda create -n nbmpo python=3.10 -y
# conda activate nbmpo

# 3) Install dependencies
pip install --upgrade pip
pip install -r requirements.txt

# 4) Launch Jupyter
jupyter lab     # or: jupyter notebook

# 5) Open index.ipynb and "Run All Cells"
```

### B) Run in the browser (Binder)

[![Launch Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/me2d09/nbmpo-diffraction-paper/HEAD?labpath=index.ipynb)

Click the badge to launch a JupyterLab session in your browser (configured from `requirements.txt` and `runtime.txt` for Python 3.10), then open **`index.ipynb`** and run all cells.


### Key Python dependencies

| Package                    | Purpose                                        |
| -------------------------- | ---------------------------------------------- |
| `numpy`, `scipy`, `pandas` | numerics and data wrangling                    |
| `matplotlib`               | publication graphics (vector + raster)         |
| `ufit`                     | neutron data loading and fitting               |
| `periodictable`            | molar masses                                   |
| `illdata`                  | downloading data from ILL portal (optional)    |
| `LongHCPulse`              | treating heat capacity with long pulse method  |

The complete, version-pinned list is in `requirements.txt`.

---

## How the notebook is organized

| Notebook section | What happens |
| --- | --- |
| **Neutron measurement on D23** | Defines the D23 NUMOR ranges for both field directions; (optionally) downloads data via `IllSftp` (guarded by `if False:`); loads the scans and composes a multi-panel plot of the neutron data with field-direction labels. Outputs `fig2.pdf` and `fig3.pdf`. |
| **Specific heat capacity measurement** | Loads PPMS raw files from `rawhc/` and processes them with `LongHCPulse`; interpolates onto a \(T\)–\(H\) grid and renders heat-map plots (imshow) with appropriate axis limits. Outputs `fig4_nolines.pdf`, lines to guide the eye in the paper are added manually. |
| **Susceptibility** | Loads calculated magnetization curves from `rawsim/` and plots \(M(H)\) for \(H\parallel c\) and \(H\perp c\); adds an “exp. data [14]” proxy to the legend; saves the plot as `fig6_nolines.pdf`, lines to guide the eye in the paper are added manually. |
| **Specific heat from the model** | Renders model-based \(T\)–\(H\) colormaps using `plot_phase_diagram` (property `mean-specific-heat`) from folders in `2dmodel/` and `3dmodel/` (pattern `"{dim}dmodel/output-Field8-0T-B{bdir}_high_res_{suffix}"`), arranged in a 2×2 subplot layout, `fig7.pdf`. |
| **Showing magnetic structures** | Composes side-by-side panels of model vs. experiment by loading pre-rendered PNGs (helper `_load_png`), cropping ROIs, and adding compass/labels and phase names (`phase_naming`), for both field orientations and phases. Outputs `fig8.pdf`, `fig5.pdf` and `fig-s3.pdf`. |

Every cell is commented so the logic is transparent; feel free to
`Shift-Enter` step-wise if you prefer.

---

## Authors & contributions

* **F. J. dos Santos** — *Theory & simulations*: implemented and ran all model calculations; produced the Monte-Carlo / spin-dynamics outputs used in `2dmodel/` and `3dmodel/`.
* **David Sviták** ([@DavidSvitak](https://github.com/DavidSvitak)) — *Thermodynamics*: performed all specific-heat measurements and provided the PPMS datasets in `rawhc/`.
* **Petr Čermák** — *Integration & visualization*: integrated datasets, authored the notebook workflow, generated figures, and carried out plotting and fit/refinement steps; repository maintenance.


---

## License and citation

All code and processed data are released under the
[Creative Commons BY 4.0](LICENSE) license.
If you find this repository useful, please cite the journal and data article once
published.

---

Enjoy the science!