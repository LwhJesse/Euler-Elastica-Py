### 🇪🇸 Spanish (`README.es.md`)
# Elastica de Euler & Analizador de Vigas No Lineales
[English](README.md) | [Français](README.fr.md) |  [Español](README.es.md) |[繁體中文](README.zh-TW.md) | [简体中文](README.zh-CN.md) | [Deutsch](README.de.md) | [日本語](README.ja.md) | [한국어](README.ko.md)
<div align="center">
  <img src="results/batch_runs/3d_renders/3D_Stress_FEM_Case8_F-20.0N.jpg" width="800" alt="Render isométrico 3D de una viga no lineal">
  <p><em>Gran deformación de una viga simplemente apoyada bajo una carga puntual extrema en el centro (FEM Corrotacional).</em></p>
</div>

Un framework de mecánica computacional altamente robusto, basado en Python, diseñado para resolver, simular y validar **vigas geométricamente no lineales con grandes deformaciones (Elastica de Euler)**.

Mediante la verificación cruzada de tres dimensiones matemáticas distintas (Analítica, método de disparo de Runge-Kutta y FEM Corrotacional), esta herramienta mapea explícitamente los límites numéricos exactos entre las suposiciones lineales y la realidad no lineal.

## 🧮 Base Teórica

La mecánica lineal tradicional (Teoría de vigas de Euler-Bernoulli) simplifica la ecuación de curvatura exacta asumiendo rotaciones infinitesimales ($w' \approx 0$), lo que lleva a sobreestimaciones severas de la deflexión bajo cargas extremas. Este framework resuelve fundamentalmente el problema de la **Elastica de Euler** al retener la no linealidad geométrica completa:

$$ \kappa = \frac{M(x)}{EI} = \frac{w''}{(1 + (w')^2)^{3/2}} $$

Esta ecuación diferencial no lineal y altamente acoplada carece de soluciones analíticas de forma cerrada para la mayoría de los casos de carga complejos, lo que requiere los solucionadores numéricos avanzados implementados en este proyecto.

## ⚙️ Detalles Técnicos y Algoritmos

El framework logra una validación de alta fidelidad a través de tres motores de solución independientes y un mapeo espacial avanzado:

### 1. Triple Verificación Cruzada Matemática y Alineación Lagrangiana
* **Línea de Base Lineal:** Soluciones de forma cerrada basadas en la teoría clásica de pequeñas deflexiones de Euler-Bernoulli. Esto sirve como la línea de base comparativa para cuantificar el umbral de divergencia no lineal.
* **Motor No Lineal Exacto:** Integración por Runge-Kutta de las ecuaciones diferenciales exactas de la Elastica de Euler.
* **Estándar de Oro (FEM):** OpenSeesPy utilizando formulaciones geométricas no lineales `Corotational`.
* **Alineación Lagrangiana Átomo a Átomo:** Debido a la flexión extrema, la proyección horizontal de la viga se reduce drásticamente. Todas las evaluaciones de error y visualizaciones se realizan estrictamente en el **sistema de coordenadas Lagrangiano** (basado en la longitud de arco inicial $s$ a través de `scipy.interpolate.interp1d`). Esto resuelve eficazmente los problemas de desalineación espacial inherentes a las deformaciones geométricas extremas en coordenadas Eulerianas.

### 2. Motor del Método de Disparo de Runge-Kutta (RK)
* **Formulación en Espacio de Estados:** Transforma la ecuación diferencial de orden superior en un sistema de EDOs de 1er orden, estableciendo el vector de estado $\mathbf{y} =[w, \theta, M, V, N]^T$.
* **Conversión de Problema de Valor de Frontera (BVP) a Problema de Valor Inicial (IVP):** Resuelve el Problema de Valor de Frontera (BVP) mediante el Método de Disparo. Utiliza el algoritmo de **Levenberg-Marquardt (`lm`)** (`scipy.optimize.least_squares`) para refinar iterativamente las conjeturas iniciales, previniendo eficazmente la singularidad de la matriz Jacobiana en regiones de flexión profunda.
* **Manejo de Discontinuidades:** Implementa condiciones de coincidencia de frontera continuas a trozos para resolver suavemente los saltos de fuerza/momento internos.

### 3. Motor de Búsqueda de Fronteras Automatizado
* **Algoritmo de Búsqueda de Raíces:** Integra el **método de Brent** (`scipy.optimize.brentq`) para buscar dinámicamente la carga aplicada exacta (o la posición de la carga) donde el error relativo entre el modelo analítico lineal y el modelo FEM no lineal alcanza un umbral estricto del **5%**.

### 4. Renderizado 3D de Alta Fidelidad sin Dependencias
* Evita bibliotecas pesadas de visualización científica 3D (como VTK, Mayavi o ParaView). Reconstruye matemáticamente elementos de viga 1D en sólidos físicos 3D y mapea los tensores de tensión equivalentes en las superficies, generando contornos de tensión isométricos 3D listos para publicación utilizando **puramente Matplotlib**.

## 📦 Configuración del Entorno

Este proyecto es una colección de scripts de computación en Python. Puedes ejecutarlos directamente una vez que se satisfagan las dependencias.

### Opción A: Configuración Estándar (Windows / macOS / Basado en Debian / Basado en Red Hat)
Para entornos de sistema operativo estándar, usar un entorno virtual es opcional pero recomendado para evitar conflictos de dependencias.
```bash
# 1. Clonar el repositorio
git clone https://github.com/LwhJesse/Euler-Elastica-Py.git
cd Euler-Elastica-Py

# 2. (Opcional) Crear y activar un entorno virtual
python3 -m venv venv
source venv/bin/activate  # En Windows usar: venv\Scripts\activate

# 3. Instalar dependencias principales
pip install --upgrade pip
pip install -r requirements.txt
```
> **⚠️ Nota para usuarios de macOS con Apple Silicon (M1/M2/M3):** `openseespy` es un framework envuelto en C++. Si `pip` no logra encontrar una "wheel" compatible con ARM64, es posible que necesites ejecutar tu terminal usando Rosetta 2 o consultar la [documentación oficial de instalación de OpenSeesPy](https://openseespydoc.readthedocs.io/en/latest/src/installation.html).

### Opción B: Linux Basado en Arch (AUR)
Si estás en un sistema Linux basado en Arch, las instalaciones globales de `pip` se gestionan externamente (PEP 668). Puedes omitir el entorno virtual de forma segura e instalar las dependencias directamente a través del gestor de paquetes del sistema.
*(Nota: El paquete de AUR `python-openseespy` es mantenido oficialmente por el autor de este repositorio [@LwhJesse](https://aur.archlinux.org/packages/python-openseespy)).*
```bash
# 1. Clonar el repositorio
git clone https://github.com/LwhJesse/Euler-Elastica-Py.git
cd Euler-Elastica-Py

# 2. Instalar dependencias vía pacman y tu ayudante de AUR (ej. yay)
sudo pacman -S python-numpy python-pandas python-scipy python-matplotlib python-rich
yay -S python-openseespy
```

## 🚀 Guía de Uso

El proyecto se basa en una única fuente de verdad para la configuración de parámetros ubicada en `core/config.py`.

### 1. Análisis de un Solo Caso
Ejecuta un solo caso para generar gráficos de comparación Lagrangianos 2D:
```bash
python run_single.py
```
Genera un contorno de tensión isométrico 3D para la configuración actual:
```bash
python run_3d_render.py
```

### 2. Ejecución por Lotes y Benchmarking
Ejecuta los 10 casos de benchmark predefinidos para generar un conjunto completo de validación:
```bash
python run_batch.py
python run_batch_3d.py
```
<details>
<summary><b>Haz clic para Ver: Alineación Multifísica Exacta RK vs FEM (Caso 8)</b></summary>
<br>
<div align="center">
  <img src="results/batch_runs/2d_plots/case8-1_comparison.jpg" width="800" alt="Gráfico de Comparación 2D">
  <p><em>Observa cómo la solución analítica lineal (verde) diverge drásticamente, mientras que RK (rojo) y FEM (negro) se alinean perfectamente.</em></p>
</div>
</details>

### 3. Búsqueda de Frontera Crítica (Clúster en Vivo)
Utiliza un panel de terminal en vivo impulsado por `rich` para mapear los límites de la fase de error no lineal del 5% en todos los núcleos de la CPU simultáneamente:
```bash
python run_multiprocess.py
```
<details>
<summary><b>Haz clic para Ver: Mapa de Zona Segura/Peligrosa Bivariable (Caso 3)</b></summary>
<br>
<div align="center">
  <img src="results/boundary_analysis/pure_critical_boundary_fem_case3.jpg" width="600" alt="Diagrama de Envolvente de Frontera">
  <p><em>La frontera crítica roja dicta explícitamente cuándo se debe adoptar el modelo geométricamente no lineal.</em></p>
</div>
</details>

## 🛠️ Casos de Carga de Benchmark Soportados

| ID de Caso | Condición de Frontera | Tipo de Carga | Variables Activas |
| :---: | :--- | :--- | :--- |
| **1** | Voladizo | Momento Flector en Extremo | $M_e$ |
| **2** | Voladizo | Fuerza Concentrada en Extremo | $F$ |
| **3** | Voladizo | Fuerza Concentrada Intermedia | $F, a$ |
| **4** | Voladizo | Carga Uniformemente Distribuida | $q$ |
| **5** | Simplemente Apoyada | Momento en Extremo Izquierdo | $M_e$ |
| **6** | Simplemente Apoyada | Momento en Extremo Derecho | $M_e$ |
| **7** | Simplemente Apoyada | Momento Intermedio | $M_e, a$ |
| **8** | Simplemente Apoyada | Fuerza Concentrada en Punto Medio | $F$ |
| **9** | Simplemente Apoyada | Fuerza Concentrada Intermedia | $F, a$ |
| **10** | Simplemente Apoyada | Carga Uniformemente Distribuida | $q$ |

## 📁 Estructura del Directorio de Salida
Los resultados están perfectamente aislados del código fuente y se organizan automáticamente:
- `results/single_runs/` : Instantáneas y renders de configuraciones manuales.
- `results/batch_runs/` : Conjuntos de validación completos de los 10 casos de benchmark.
- `results/boundary_analysis/` : Evolución dinámica del error y curvas de umbral crítico.
- `results/mesh_convergence/` : Gráficos de verificación de independencia de la malla.

## 📄 Licencia
Este proyecto está licenciado bajo la Licencia MIT - vea el archivo [LICENSE](LICENSE) para más detalles.
