### 🇩🇪 German (`README.de.md`)
# Euler Elastica & Nichtlinearer Balken-Analysator
[English](README.md) | [Français](README.fr.md) |  [Español](README.es.md) |[繁體中文](README.zh-TW.md) | [简体中文](README.zh-CN.md) | [Deutsch](README.de.md) | [日本語](README.ja.md) | [한국어](README.ko.md)
<div align="center">
  <img src="results/batch_runs/3d_renders/3D_Stress_FEM_Case8_F-20.0N.jpg" width="800" alt="3D-isometrische Darstellung eines nichtlinearen Balkens">
  <p><em>Große Verformung eines einfach gelagerten Balkens unter extremer Mittellast (Korotations-FEM).</em></p>
</div>

Ein hochrobustes, Python-basiertes Framework für computergestützte Mechanik, das für die Lösung, Simulation und Validierung von **geometrisch nichtlinearen Balken mit großen Verformungen (Euler Elastica)** entwickelt wurde.

Durch die Kreuzvalidierung von drei unterschiedlichen mathematischen Dimensionen (analytisch, Runge-Kutta-Shooting-Verfahren und Korotations-FEM) bildet dieses Werkzeug explizit die exakten numerischen Grenzen zwischen linearen Annahmen und der nichtlinearen Realität ab.

## 🧮 Theoretische Grundlagen

Die traditionelle lineare Mechanik (Euler-Bernoulli-Balkentheorie) vereinfacht die exakte Krümmungsgleichung durch die Annahme infinitesimaler Rotationen ($w' \approx 0$), was zu erheblichen Überschätzungen der Durchbiegung unter extremen Lasten führt. Dieses Framework löst das **Euler-Elastica**-Problem grundlegend, indem es die volle geometrische Nichtlinearität beibehält:

$$ \kappa = \frac{M(x)}{EI} = \frac{w''}{(1 + (w')^2)^{3/2}} $$

Diese stark gekoppelte, nichtlineare Differentialgleichung besitzt für die meisten komplexen Lastfälle keine geschlossenen analytischen Lösungen, was die in diesem Projekt implementierten fortgeschrittenen numerischen Löser erforderlich macht.

## ⚙️ Technische Details & Algorithmen

Das Framework erreicht eine hochpräzise Validierung durch drei unabhängige Lösungs-Engines und eine fortschrittliche räumliche Abbildung:

### 1. Dreifache mathematische Kreuzvalidierung & Lagrange-Ausrichtung
* **Lineare Basis:** Geschlossene Lösungen basierend auf der klassischen Euler-Bernoulli-Theorie für kleine Durchbiegungen. Dies dient als Vergleichsbasis, um die Schwelle der nichtlinearen Divergenz zu quantifizieren.
* **Exakte nichtlineare Engine:** Runge-Kutta-Integration der exakten Euler-Elastica-Differentialgleichungen.
* **FEM-Goldstandard:** OpenSeesPy unter Verwendung von `Corotational` geometrisch nichtlinearen Formulierungen.
* **Atom-zu-Atom-Lagrange-Ausrichtung:** Aufgrund der extremen Biegung schrumpft die horizontale Projektion des Balkens drastisch. Alle Fehlerbewertungen und Visualisierungen werden streng im **Lagrange-Koordinatensystem** durchgeführt (basierend auf der anfänglichen Bogenlänge $s$ mittels `scipy.interpolate.interp1d`). Dies löst effektiv die Probleme der räumlichen Fehlausrichtung, die bei extremen geometrischen Verformungen in Euler-Koordinaten auftreten.

### 2. RK-Shooting-Verfahren-Engine
* **Zustandsraumformulierung:** Wandelt die Differentialgleichung höherer Ordnung in ein System von gewöhnlichen Differentialgleichungen 1. Ordnung um und definiert den Zustandsvektor $\mathbf{y} =[w, \theta, M, V, N]^T$.
* **Konvertierung von Randwertproblem (BVP) zu Anfangswertproblem (IVP):** Löst das Randwertproblem (BVP) mittels des Shooting-Verfahrens. Es verwendet den **Levenberg-Marquardt (`lm`) Algorithmus** (`scipy.optimize.least_squares`), um Anfangsschätzungen iterativ zu verfeinern und so die Singularität der Jacobi-Matrix in Bereichen starker Biegung zu verhindern.
* **Behandlung von Diskontinuitäten:** Implementiert stückweise stetige Randanpassungsbedingungen, um interne Kraft-/Momentsprünge reibungslos aufzulösen.

### 3. Automatisierte Grenzsuch-Engine
* **Wurzelfindungsalgorithmus:** Integriert die **Brent-Methode** (`scipy.optimize.brentq`), um dynamisch die exakte aufgebrachte Last (oder Lastposition) zu finden, bei der der relative Fehler zwischen dem linearen analytischen Modell und dem nichtlinearen FEM-Modell eine strikte **5%-Schwelle** erreicht.

### 4. Abhängigkeitsfreies, hochpräzises 3D-Rendering
* Vermeidet schwere 3D-wissenschaftliche Visualisierungsbibliotheken (wie VTK, Mayavi oder ParaView). Es rekonstruiert mathematisch 1D-Balkenelemente in 3D-Festkörper und bildet die äquivalenten Spannungstensoren auf die Oberflächen ab, wodurch publikationsreife 3D-isometrische Spannungskonturen mit **reinem Matplotlib** erzeugt werden.

## 📦 Umgebungseinrichtung

Dieses Projekt ist eine Sammlung von computergestützten Python-Skripten. Sie können sie direkt ausführen, sobald die Abhängigkeiten erfüllt sind.

### Option A: Standard-Setup (Windows / macOS / Debian-basiert / Red Hat-basiert)
Für Standard-Betriebssystemumgebungen ist die Verwendung einer virtuellen Umgebung optional, wird aber empfohlen, um Abhängigkeitskonflikte zu vermeiden.
```bash
# 1. Das Repository klonen
git clone https://github.com/LwhJesse/Euler-Elastica-Py.git
cd Euler-Elastica-Py

# 2. (Optional) Eine virtuelle Umgebung erstellen und aktivieren
python3 -m venv venv
source venv/bin/activate  # Unter Windows: venv\Scripts\activate

# 3. Kernabhängigkeiten installieren
pip install --upgrade pip
pip install -r requirements.txt
```
> **⚠️ Hinweis für macOS Apple Silicon (M1/M2/M3) Benutzer:** `openseespy` ist ein in C++ gewrapptes Framework. Wenn `pip` kein kompatibles ARM64-Wheel findet, müssen Sie Ihr Terminal möglicherweise mit Rosetta 2 ausführen oder die [offizielle OpenSeesPy-Installationsdokumentation](https://openseespydoc.readthedocs.io/en/latest/src/installation.html) konsultieren.

### Option B: Arch-basiertes Linux (AUR)
Wenn Sie ein Arch-basiertes Linux verwenden, werden globale `pip`-Installationen extern verwaltet (PEP 668). Sie können die virtuelle Umgebung sicher überspringen und die Abhängigkeiten direkt über den System-Paketmanager installieren.
*(Hinweis: Das AUR-Paket `python-openseespy` wird offiziell vom Autor dieses Repositories [@LwhJesse](https://aur.archlinux.org/packages/python-openseespy) gepflegt).*
```bash
# 1. Das Repository klonen
git clone https://github.com/LwhJesse/Euler-Elastica-Py.git
cd Euler-Elastica-Py

# 2. Abhängigkeiten über pacman und Ihren AUR-Helfer (z.B. yay) installieren
sudo pacman -S python-numpy python-pandas python-scipy python-matplotlib python-rich
yay -S python-openseespy
```

## 🚀 Anwendungsleitfaden

Das Projekt stützt sich auf eine einzige Quelle der Wahrheit für die Parameterkonfiguration, die sich in `core/config.py` befindet.

### 1. Einzelfallanalyse
Führen Sie einen einzelnen Fall aus, um 2D-Lagrange-Vergleichsdiagramme zu erstellen:
```bash
python run_single.py
```
Erstellen Sie eine 3D-isometrische Spannungskontur für die aktuelle Konfiguration:
```bash
python run_3d_render.py
```

### 2. Stapelverarbeitung & Benchmarking
Führen Sie alle 10 vordefinierten Benchmark-Fälle aus, um eine umfassende Validierungssuite zu erstellen:
```bash
python run_batch.py
python run_batch_3d.py
```
<details>
<summary><b>Klicken zum Anzeigen: Exakte RK vs. FEM Multiphysik-Ausrichtung (Fall 8)</b></summary>
<br>
<div align="center">
  <img src="results/batch_runs/2d_plots/case8-1_comparison.jpg" width="800" alt="2D-Vergleichsdiagramm">
  <p><em>Beachten Sie, wie die lineare analytische Lösung (grün) drastisch abweicht, während RK (rot) und FEM (schwarz) perfekt übereinstimmen.</em></p>
</div>
</details>

### 3. Kritische Grenzsuche (Live-Cluster)
Nutzen Sie ein `rich`-betriebenes Live-Terminal-Dashboard, um die 5%-Phasengrenzen des nichtlinearen Fehlers auf allen CPU-Kernen gleichzeitig abzubilden:
```bash
python run_multiprocess.py
```
<details>
<summary><b>Klicken zum Anzeigen: Bivariates Sicherheits-/Gefahrenzonen-Diagramm (Fall 3)</b></summary>
<br>
<div align="center">
  <img src="results/boundary_analysis/pure_critical_boundary_fem_case3.jpg" width="600" alt="Grenz-Hüllkurven-Diagramm">
  <p><em>Die rote kritische Grenze schreibt explizit vor, wann das geometrisch nichtlineare Modell angewendet werden muss.</em></p>
</div>
</details>

## 🛠️ Unterstützte Benchmark-Lastfälle

| Fall-ID | Randbedingung | Lasttyp | Aktive Variablen |
| :---: | :--- | :--- | :--- |
| **1** | Kragarm | End-Biegemoment | $M_e$ |
| **2** | Kragarm | End-Einzellast | $F$ |
| **3** | Kragarm | Zwischen-Einzellast | $F, a$ |
| **4** | Kragarm | Gleichmäßig verteilte Last | $q$ |
| **5** | Einfach gelagert | Moment am linken Ende | $M_e$ |
| **6** | Einfach gelagert | Moment am rechten Ende | $M_e$ |
| **7** | Einfach gelagert | Zwischen-Moment | $M_e, a$ |
| **8** | Einfach gelagert | Mittelpunkt-Einzellast | $F$ |
| **9** | Einfach gelagert | Zwischen-Einzellast | $F, a$ |
| **10** | Einfach gelagert | Gleichmäßig verteilte Last | $q$ |

## 📁 Ausgabeverzeichnisstruktur
Die Ergebnisse sind perfekt vom Quellcode isoliert und werden automatisch organisiert:
- `results/single_runs/` : Snapshots und Renderings manueller Konfigurationen.
- `results/batch_runs/` : Vollständige Validierungssuites der 10 Benchmark-Fälle.
- `results/boundary_analysis/` : Dynamische Fehlerentwicklung und kritische Schwellenkurven.
- `results/mesh_convergence/` : Diagramme zur Überprüfung der Netzunabhängigkeit.

## 📄 Lizenz
Dieses Projekt ist unter der MIT-Lizenz lizenziert - siehe die Datei [LICENSE](LICENSE) für Details.