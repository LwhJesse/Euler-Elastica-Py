### 🇺🇸 English(`README.md`)
# Euler Elastica & Nonlinear Beam Analyzer
[English](README.md) | [Français](README.fr.md) |  [Español](README.es.md) |[繁體中文](README.zh-TW.md) | [简体中文](README.zh-CN.md) | [Deutsch](README.de.md) | [日本語](README.ja.md) | [한국어](README.ko.md)
<div align="center">
  <img src="results/batch_runs/3d_renders/3D_Stress_FEM_Case8_F-20.0N.jpg" width="800" alt="3D Isometric Render of Nonlinear Beam">
  <p><em>Large deformation of a simply supported beam under extreme midpoint load (Corotational FEM).</em></p>
</div>

A highly robust, Python-based computational mechanics framework designed to solve, simulate, and validate **large-deformation geometrically nonlinear beams (Euler Elastica)**. 

By cross-verifying three distinct mathematical dimensions (Analytical, Runge-Kutta Shooting, and Corotational FEM), this tool explicitly maps the exact numerical boundaries between linear assumptions and nonlinear reality.

## 🧮 Theoretical Basis

Traditional linear mechanics (Euler-Bernoulli Beam Theory) simplifies the exact curvature equation by assuming infinitesimal rotations ($w' \approx 0$), leading to severe overestimations of deflection under extreme loads. This framework fundamentally solves the **Euler Elastica** problem by retaining the full geometric nonlinearity:

$$ \kappa = \frac{M(x)}{EI} = \frac{w''}{(1 + (w')^2)^{3/2}} $$

This highly coupled, nonlinear differential equation lacks closed-form analytical solutions for most complex load cases, necessitating the advanced numerical solvers implemented in this project.

## ⚙️ Technical Details & Algorithms

The framework achieves high-fidelity validation through three independent solver engines and advanced spatial mapping:

### 1. Tri-fold Mathematical Cross-Validation & Lagrangian Alignment
* **Linear Baseline:** Closed-form solutions based on classical Euler-Bernoulli small-deflection theory. This serves as the comparative baseline to quantify the threshold of nonlinear divergence.
* **Exact Nonlinear Engine:** Runge-Kutta integration of the exact Euler Elastica differential equations.
* **FEM Gold Standard:** OpenSeesPy utilizing `Corotational` geometric nonlinear formulations.
* **Atom-to-Atom Lagrangian Alignment:** Due to extreme bending, the beam's horizontal projection shrinks drastically. All error evaluations and visualizations are strictly performed in the **Lagrangian coordinate system** (based on the initial arc length $s$ via `scipy.interpolate.interp1d`). This effectively resolves the spatial misalignment issues inherent to extreme geometric deformations in Eulerian coordinates.

### 2. RK Shooting Method Engine
* **State-Space Formulation**: Transforms the higher-order differential equation into a system of 1st-order ODEs, establishing the state vector $\mathbf{y} =[w, \theta, M, V, N]^T$.
* **BVP to IVP Conversion**: Solves the Boundary Value Problem (BVP) via the Shooting Method. It uses the **Levenberg-Marquardt (`lm`) algorithm** (`scipy.optimize.least_squares`) to iteratively refine initial guesses, effectively preventing Jacobian matrix singularity in deep-bending regions.
* **Discontinuity Handling**: Implements piecewise continuous boundary matching conditions to smoothly resolve internal force/moment jumps.

### 3. Automated Boundary Search Engine
* **Root-Finding Algorithm**: Integrates **Brent's method** (`scipy.optimize.brentq`) to dynamically hunt for the exact applied load (or load position) where the relative error between the linear analytical model and the nonlinear FEM model hits a strict **5% threshold**.

### 4. Zero-Dependency High-Fidelity 3D Rendering
* Avoids heavy 3D scientific visualization libraries (like VTK, Mayavi, or ParaView). It mathematically reconstructs 1D beam elements into 3D physical solids and maps the equivalent stress tensors onto the surfaces, generating publication-ready 3D isometric stress contours using **pure Matplotlib**.

## 📦 Environment Setup

This project is a collection of computational Python scripts. You can run them directly once the dependencies are satisfied.

### Option A: Standard Setup (Windows / macOS / Debian-based / Red Hat-based)
For standard OS environments, using a virtual environment is optional but recommended to avoid dependency conflicts.
```bash
# 1. Clone the repository
git clone https://github.com/LwhJesse/Euler-Elastica-Py.git
cd Euler-Elastica-Py

# 2. (Optional) Create and activate a virtual environment
python3 -m venv venv
source venv/bin/activate  # On Windows use: venv\Scripts\activate

# 3. Install core dependencies
pip install --upgrade pip
pip install -r requirements.txt
```
> **⚠️ Note for macOS Apple Silicon (M1/M2/M3) users:** `openseespy` is a C++ wrapped framework. If `pip` fails to find a compatible ARM64 wheel, you may need to run your terminal using Rosetta 2 or refer to the [official OpenSeesPy installation documentation](https://openseespydoc.readthedocs.io/en/latest/src/installation.html).

### Option B: Arch-based Linux (AUR)
If you are on Arch-based Linux, global `pip` installations are externally managed (PEP 668). You can safely skip the virtual environment and install the dependencies directly via the system package manager. 
*(Note: The `python-openseespy` AUR package is officially maintained by the author of this repository [@LwhJesse](https://aur.archlinux.org/packages/python-openseespy)).*
```bash
# 1. Clone the repository
git clone https://github.com/LwhJesse/Euler-Elastica-Py.git
cd Euler-Elastica-Py

# 2. Install dependencies via pacman and your AUR helper (e.g., yay)
sudo pacman -S python-numpy python-pandas python-scipy python-matplotlib python-rich
yay -S python-openseespy
```

## 🚀 Usage Guide

The project relies on a single source of truth for parameter configuration located in `core/config.py`.

### 1. Single Case Analysis
Run a single case to generate 2D Lagrangian comparison plots:
```bash
python run_single.py
```
Generate a 3D isometric stress contour for the current configuration:
```bash
python run_3d_render.py
```

### 2. Batch Execution & Benchmarking
Execute all 10 predefined benchmark cases to generate a comprehensive validation suite:
```bash
python run_batch.py
python run_batch_3d.py
```
<details>
<summary><b>Click to View: Exact RK vs FEM Multi-physics Alignment (Case 8)</b></summary>
<br>
<div align="center">
  <img src="results/batch_runs/2d_plots/case8-1_comparison.jpg" width="800" alt="2D Comparison Plot">
  <p><em>Notice how the linear analytical solution (green) drastically diverges, while RK (red) and FEM (black) perfectly align.</em></p>
</div>
</details>

### 3. Critical Boundary Search (Live Cluster)
Utilize a `rich`-powered live terminal dashboard to map out the 5% nonlinear error phase boundaries across all CPU cores concurrently:
```bash
python run_multiprocess.py
```
<details>
<summary><b>Click to View: Bi-variable Safe/Danger Zone Map (Case 3)</b></summary>
<br>
<div align="center">
  <img src="results/boundary_analysis/pure_critical_boundary_fem_case3.jpg" width="600" alt="Boundary Envelope Diagram">
  <p><em>The red critical boundary explicitly dictates when the geometrically nonlinear model must be adopted.</em></p>
</div>
</details>

## 🛠️ Supported Benchmark Load Cases

| Case ID | Boundary Condition | Load Type | Active Variables |
| :---: | :--- | :--- | :--- |
| **1** | Cantilever | End Bending Moment | $M_e$ |
| **2** | Cantilever | End Concentrated Force | $F$ |
| **3** | Cantilever | Intermediate Concentrated Force | $F, a$ |
| **4** | Cantilever | Uniformly Distributed Load | $q$ |
| **5** | Simply Supported | Left End Moment | $M_e$ |
| **6** | Simply Supported | Right End Moment | $M_e$ |
| **7** | Simply Supported | Intermediate Moment | $M_e, a$ |
| **8** | Simply Supported | Midpoint Concentrated Force | $F$ |
| **9** | Simply Supported | Intermediate Concentrated Force | $F, a$ |
| **10** | Simply Supported | Uniformly Distributed Load | $q$ |

## 📁 Output Directory Structure
Results are perfectly isolated from the source code and organized automatically:
- `results/single_runs/` : Snapshots and renders of manual configurations.
- `results/batch_runs/` : Full validation suites of the 10 benchmark cases.
- `results/boundary_analysis/` : Dynamic error evolution and critical threshold curves.
- `results/mesh_convergence/` : Mesh independence verification plots.

## 📄 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
