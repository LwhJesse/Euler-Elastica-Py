### 🇫🇷 French (`README_fr.md`)

# Analyseur de Poutres Non-linéaires & Euler Elastica
[English](README.md) | [Français](README.fr.md) |  [Español](README.es.md) |[繁體中文](README.zh-TW.md) | [简体中文](README.zh-CN.md) | [Deutsch](README.de.md) | [日本語](README.ja.md) | [한국어](README.ko.md)
<div align="center">
  <img src="results/batch_runs/3d_renders/3D_Stress_FEM_Case8_F-20.0N.jpg" width="800" alt="Rendu isométrique 3D d'une poutre non-linéaire">
  <p><em>Grande déformation d'une poutre simplement appuyée sous une charge ponctuelle extrême au milieu (FEM corotationnelle).</em></p>
</div>

Un framework de mécanique computationnelle très robuste, basé sur Python, conçu pour résoudre, simuler et valider des **poutres géométriquement non-linéaires à grande déformation (Euler Elastica)**.

En vérifiant de manière croisée trois dimensions mathématiques distinctes (Analytique, méthode de tir de Runge-Kutta, et FEM corotationnelle), cet outil cartographie explicitement les frontières numériques exactes entre les hypothèses linéaires et la réalité non-linéaire.

## 🧮 Bases Théoriques

La mécanique linéaire traditionnelle (théorie des poutres d'Euler-Bernoulli) simplifie l'équation de courbure exacte en supposant des rotations infinitésimales ($w' \approx 0$), ce qui conduit à de sévères surestimations de la flèche sous des charges extrêmes. Ce framework résout fondamentalement le problème de l'**Euler Elastica** en conservant la non-linéarité géométrique complète :

$$ \kappa = \frac{M(x)}{EI} = \frac{w''}{(1 + (w')^2)^{3/2}} $$

Cette équation différentielle non-linéaire et fortement couplée n'a pas de solution analytique de forme fermée pour la plupart des cas de charge complexes, ce qui nécessite les solveurs numériques avancés mis en œuvre dans ce projet.

## ⚙️ Détails Techniques & Algorithmes

Le framework atteint une validation de haute fidélité grâce à trois moteurs de résolution indépendants et une cartographie spatiale avancée :

### 1. Triple Validation Croisée Mathématique & Alignement Lagrangien
* **Référence Linéaire :** Solutions de forme fermée basées sur la théorie classique des petites déformations d'Euler-Bernoulli. Ceci sert de référence comparative pour quantifier le seuil de divergence non-linéaire.
* **Moteur Non-linéaire Exact :** Intégration par Runge-Kutta des équations différentielles exactes de l'Euler Elastica.
* **Référence Absolue (FEM) :** OpenSeesPy utilisant des formulations géométriques non-linéaires `Corotational`.
* **Alignement Lagrangien Point par Point :** En raison de la flexion extrême, la projection horizontale de la poutre se réduit considérablement. Toutes les évaluations d'erreur et visualisations sont strictement effectuées dans le **système de coordonnées Lagrangien** (basé sur la longueur d'arc initiale $s$ via `scipy.interpolate.interp1d`). Cela résout efficacement les problèmes de désalignement spatial inhérents aux déformations géométriques extrêmes en coordonnées Eulériennes.

### 2. Moteur de Méthode de Tir de Runge-Kutta (RK)
* **Formulation en Espace d'États :** Transforme l'équation différentielle d'ordre supérieur en un système d'EDO de 1er ordre, établissant le vecteur d'état $\mathbf{y} =[w, \theta, M, V, N]^T$.
* **Conversion de Problème aux Limites (BVP) en Problème à Valeur Initiale (IVP) :** Résout le problème aux limites (BVP) via la méthode de tir. Il utilise l'algorithme de **Levenberg-Marquardt (`lm`)** (`scipy.optimize.least_squares`) pour affiner itérativement les estimations initiales, empêchant ainsi la singularité de la matrice Jacobienne dans les zones de forte flexion.
* **Gestion des Discontinuités :** Met en œuvre des conditions de raccordement aux limites continues par morceaux pour gérer sans heurt les sauts de forces/moments internes.

### 3. Moteur de Recherche Automatisée de Frontières
* **Algorithme de Recherche de Racine :** Intègre la **méthode de Brent** (`scipy.optimize.brentq`) pour rechercher dynamiquement la charge appliquée exacte (ou la position de la charge) où l'erreur relative entre le modèle analytique linéaire et le modèle FEM non-linéaire atteint un seuil strict de **5%**.

### 4. Rendu 3D Haute Fidélité Sans Dépendances
* Évite les bibliothèques de visualisation scientifique 3D lourdes (comme VTK, Mayavi ou ParaView). Il reconstruit mathématiquement des éléments de poutre 1D en solides physiques 3D et cartographie les tenseurs de contraintes équivalentes sur les surfaces, générant des contours de contraintes isométriques 3D prêts pour la publication en utilisant **uniquement Matplotlib**.

## 📦 Configuration de l'Environnement

Ce projet est une collection de scripts de calcul en Python. Vous pouvez les exécuter directement une fois les dépendances satisfaites.

### Option A : Installation Standard (Windows / macOS / Basé sur Debian / Basé sur Red Hat)
Pour les environnements de système d'exploitation standard, l'utilisation d'un environnement virtuel est facultative mais recommandée pour éviter les conflits de dépendances.
```bash
# 1. Cloner le dépôt
git clone https://github.com/LwhJesse/Euler-Elastica-Py.git
cd Euler-Elastica-Py

# 2. (Optionnel) Créer et activer un environnement virtuel
python3 -m venv venv
source venv/bin/activate  # Sous Windows, utiliser : venv\Scripts\activate

# 3. Installer les dépendances principales
pip install --upgrade pip
pip install -r requirements.txt
```
> **⚠️ Note pour les utilisateurs de macOS avec Apple Silicon (M1/M2/M3) :** `openseespy` est un framework en C++. Si `pip` ne parvient pas à trouver une roue compatible ARM64, vous devrez peut-être exécuter votre terminal avec Rosetta 2 ou consulter la [documentation officielle d'installation d'OpenSeesPy](https://openseespydoc.readthedocs.io/en/latest/src/installation.html).

### Option B : Linux Basé sur Arch (AUR)
Si vous êtes sur un système Linux basé sur Arch, les installations globales de `pip` sont gérées de manière externe (PEP 668). Vous pouvez ignorer l'environnement virtuel et installer les dépendances directement via le gestionnaire de paquets du système.
*(Note : Le paquet AUR `python-openseespy` est officiellement maintenu par l'auteur de ce dépôt [@LwhJesse](https://aur.archlinux.org/packages/python-openseespy)).*
```bash
# 1. Cloner le dépôt
git clone https://github.com/LwhJesse/Euler-Elastica-Py.git
cd Euler-Elastica-Py

# 2. Installer les dépendances via pacman et votre assistant AUR (ex: yay)
sudo pacman -S python-numpy python-pandas python-scipy python-matplotlib python-rich
yay -S python-openseespy
```

## 🚀 Guide d'Utilisation

Le projet s'appuie sur une source unique de vérité pour la configuration des paramètres, située dans `core/config.py`.

### 1. Analyse d'un Cas Unique
Exécutez un cas unique pour générer des graphiques de comparaison Lagrangiens 2D :
```bash
python run_single.py
```
Générez un contour de contraintes isométrique 3D pour la configuration actuelle :
```bash
python run_3d_render.py
```

### 2. Exécution par Lots & Benchmarking
Exécutez les 10 cas de benchmark prédéfinis pour générer une suite de validation complète :
```bash
python run_batch.py
python run_batch_3d.py
```
<details>
<summary><b>Cliquez pour Voir : Alignement Multi-physique RK Exact vs FEM (Cas 8)</b></summary>
<br>
<div align="center">
  <img src="results/batch_runs/2d_plots/case8-1_comparison.jpg" width="800" alt="Graphique de Comparaison 2D">
  <p><em>Notez comment la solution analytique linéaire (en vert) diverge considérablement, tandis que RK (en rouge) et FEM (en noir) s'alignent parfaitement.</em></p>
</div>
</details>

### 3. Recherche de Frontière Critique (Cluster en Direct)
Utilisez un tableau de bord de terminal en direct alimenté par `rich` pour cartographier les frontières de la phase d'erreur non-linéaire de 5% sur tous les cœurs de CPU simultanément :
```bash
python run_multiprocess.py
```
<details>
<summary><b>Cliquez pour Voir : Carte de Zone Sûre/Dangereuse à Deux Variables (Cas 3)</b></summary>
<br>
<div align="center">
  <img src="results/boundary_analysis/pure_critical_boundary_fem_case3.jpg" width="600" alt="Diagramme d'Enveloppe de Frontière">
  <p><em>La frontière critique rouge dicte explicitement quand le modèle géométriquement non-linéaire doit être adopté.</em></p>
</div>
</details>

## 🛠️ Cas de Charge de Benchmark Supportés

| ID Cas | Condition aux Limites | Type de Charge | Variables Actives |
| :---: | :--- | :--- | :--- |
| **1** | Cantilever | Moment Fléchissant en Bout | $M_e$ |
| **2** | Cantilever | Force Concentrée en Bout | $F$ |
| **3** | Cantilever | Force Concentrée Intermédiaire | $F, a$ |
| **4** | Cantilever | Charge Uniformément Répartie | $q$ |
| **5** | Simplement Appuyée | Moment à l'Extrémité Gauche | $M_e$ |
| **6** | Simplement Appuyée | Moment à l'Extrémité Droite | $M_e$ |
| **7** | Simplement Appuyée | Moment Intermédiaire | $M_e, a$ |
| **8** | Simplement Appuyée | Force Concentrée au Milieu | $F$ |
| **9** | Simplement Appuyée | Force Concentrée Intermédiaire | $F, a$ |
| **10** | Simplement Appuyée | Charge Uniformément Répartie | $q$ |

## 📁 Structure du Répertoire de Sortie
Les résultats sont parfaitement isolés du code source et organisés automatiquement :
- `results/single_runs/` : Captures et rendus des configurations manuelles.
- `results/batch_runs/` : Suites de validation complètes des 10 cas de benchmark.
- `results/boundary_analysis/` : Évolution dynamique de l'erreur et courbes de seuil critiques.
- `results/mesh_convergence/` : Graphiques de vérification de l'indépendance du maillage.

## 📄 Licence
Ce projet est sous licence MIT - voir le fichier [LICENSE](LICENSE) pour plus de détails.

