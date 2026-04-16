### 🇨🇳 Simplified Chinese (`README.zh-CN.md`)

# 欧拉弹性力学 & 非线性梁分析器
[English](README.md) | [Français](README.fr.md) |  [Español](README.es.md) |[繁體中文](README.zh-TW.md) | [简体中文](README.zh-CN.md) | [Deutsch](README.de.md) | [日本語](README.ja.md) | [한국어](README.ko.md)
<div align="center">
  <img src="results/batch_runs/3d_renders/3D_Stress_FEM_Case8_F-20.0N.jpg" width="800" alt="非线性梁的3D等轴测渲染">
  <p><em>简支梁在极端中点载荷下的大变形（共旋转FEM）。</em></p>
</div>

一个基于 Python 的高鲁棒性计算力学框架，旨在求解、模拟和验证**大变形几何非线性梁（欧拉弹性体）**。

通过对三个不同的数学维度（解析解、龙格-库塔打靶法和共旋转FEM）进行交叉验证，该工具明确地描绘了线性假设与非线性现实之间的精确数值边界。

## 🧮 理论基础

传统的线性力学（欧拉-伯努利梁理论）通过假设无穷小转角（$w' \approx 0$）来简化精确的曲率方程，这导致在极端载荷下严重高估挠度。本框架通过保留完整的几何非线性，从根本上解决了**欧拉弹性体**问题：

$$ \kappa = \frac{M(x)}{EI} = \frac{w''}{(1 + (w')^2)^{3/2}} $$

这个高度耦合的非线性微分方程对于大多数复杂的载荷情况缺乏闭合形式的解析解，因此需要本项目中实现的高级数值求解器。

## ⚙️ 技术细节与算法

该框架通过三个独立的求解器引擎和先进的空间映射实现了高保真度验证：

### 1. 三重数学交叉验证与拉格朗日对齐
* **线性基准：** 基于经典欧拉-伯努利小挠度理论的闭合形式解。这作为量化非线性发散阈值的比较基准。
* **精确非线性引擎：** 对精确的欧拉弹性体微分方程进行龙格-库塔积分。
* **FEM黄金标准：** 使用 `Corotational` 几何非线性公式的 OpenSeesPy。
* **原子级拉格朗日对齐：** 由于极端弯曲，梁的水平投影会急剧缩小。所有的误差评估和可视化都严格在**拉格朗日坐标系**中执行（基于初始弧长 $s$，通过 `scipy.interpolate.interp1d` 实现）。这有效地解决了欧拉坐标系中极端几何变形固有的空间错位问题。

### 2. RK打靶法引擎
* **状态空间公式化：** 将高阶微分方程转换为一阶常微分方程组，建立状态向量 $\mathbf{y} =[w, \theta, M, V, N]^T$。
* **边值问题(BVP)到初值问题(IVP)的转换：** 通过打靶法求解边值问题(BVP)。它使用**Levenberg-Marquardt (`lm`) 算法** (`scipy.optimize.least_squares`) 迭代地优化初始猜测，有效防止在深度弯曲区域雅可比矩阵出现奇异。
* **不连续性处理：** 实现分段连续的边界匹配条件，以平滑地解决内力/力矩的跳跃。

### 3. 自动化边界搜索引擎
* **求根算法：** 集成**Brent方法** (`scipy.optimize.brentq`)，动态寻找线性解析模型与非线性FEM模型之间相对误差达到严格**5%阈值**时的精确施加载荷（或载荷位置）。

### 4. 零依赖高保真3D渲染
* 避免了重型的3D科学可视化库（如VTK、Mayavi或ParaView）。它通过数学方式将一维梁单元重建为三维实体，并将等效应力张量映射到表面，使用**纯Matplotlib**生成可用于出版的3D等轴测应力云图。

## 📦 环境设置

本项目是一系列计算性Python脚本的集合。一旦满足依赖关系，您就可以直接运行它们。

### 选项A：标准设置 (Windows / macOS / 基于Debian / 基于Red Hat)
对于标准操作系统环境，使用虚拟环境是可选的，但建议使用以避免依赖冲突。
```bash
# 1. 克隆仓库
git clone https://github.com/LwhJesse/Euler-Elastica-Py.git
cd Euler-Elastica-Py

# 2. (可选) 创建并激活虚拟环境
python3 -m venv venv
source venv/bin/activate  # 在Windows上使用: venv\Scripts\activate

# 3. 安装核心依赖
pip install --upgrade pip
pip install -r requirements.txt
```
> **⚠️ macOS Apple Silicon (M1/M2/M3) 用户注意：** `openseespy` 是一个C++封装的框架。如果 `pip` 无法找到兼容的ARM64 wheel文件，您可能需要使用Rosetta 2运行您的终端，或参考 [OpenSeesPy官方安装文档](https://openseespydoc.readthedocs.io/en/latest/src/installation.html)。

### 选项B：基于Arch的Linux (AUR)
如果您使用的是基于Arch的Linux，全局 `pip` 安装是外部管理的 (PEP 668)。您可以安全地跳过虚拟环境，直接通过系统包管理器安装依赖。
*(注：`python-openseespy` AUR包由本仓库作者 [@LwhJesse](https://aur.archlinux.org/packages/python-openseespy) 官方维护)。*
```bash
# 1. 克隆仓库
git clone https://github.com/LwhJesse/Euler-Elastica-Py.git
cd Euler-Elastica-Py

# 2. 通过pacman和您的AUR助手 (例如yay) 安装依赖
sudo pacman -S python-numpy python-pandas python-scipy python-matplotlib python-rich
yay -S python-openseespy
```

## 🚀 使用指南

项目依赖位于 `core/config.py` 中的单一参数配置文件。

### 1. 单一案例分析
运行单个案例以生成2D拉格朗日比较图：
```bash
python run_single.py
```
为当前配置生成3D等轴测应力云图：
```bash
python run_3d_render.py
```

### 2. 批量执行与基准测试
执行所有10个预定义的基准案例，以生成一套全面的验证套件：
```bash
python run_batch.py
python run_batch_3d.py
```
<details>
<summary><b>点击查看：精确RK vs FEM多物理场对齐 (案例8)</b></summary>
<br>
<div align="center">
  <img src="results/batch_runs/2d_plots/case8-1_comparison.jpg" width="800" alt="2D比较图">
  <p><em>请注意，线性解析解（绿色）出现了剧烈偏差，而RK（红色）和FEM（黑色）则完美对齐。</em></p>
</div>
</details>

### 3. 临界边界搜索 (实时集群)
利用一个由 `rich` 驱动的实时终端仪表板，在所有CPU核心上并发地描绘出5%非线性误差的相边界：
```bash
python run_multiprocess.py
```
<details>
<summary><b>点击查看：双变量安全/危险区域图 (案例3)</b></summary>
<br>
<div align="center">
  <img src="results/boundary_analysis/pure_critical_boundary_fem_case3.jpg" width="600" alt="边界包络图">
  <p><em>红色的临界边界明确地指示了何时必须采用几何非线性模型。</em></p>
</div>
</details>

## 🛠️ 支持的基准载荷案例

| 案例ID | 边界条件 | 载荷类型 | 活动变量 |
| :---: | :--- | :--- | :--- |
| **1** | 悬臂梁 | 端部弯矩 | $M_e$ |
| **2** | 悬臂梁 | 端部集中力 | $F$ |
| **3** | 悬臂梁 | 中间集中力 | $F, a$ |
| **4** | 悬臂梁 | 均布载荷 | $q$ |
| **5** | 简支梁 | 左端力矩 | $M_e$ |
| **6** | 简支梁 | 右端力矩 | $M_e$ |
| **7** | 简支梁 | 中间力矩 | $M_e, a$ |
| **8** | 简支梁 | 中点集中力 | $F$ |
| **9** | 简支梁 | 中间集中力 | $F, a$ |
| **10** | 简支梁 | 均布载荷 | $q$ |

## 📁 输出目录结构
结果与源代码完全隔离，并自动组织：
- `results/single_runs/` : 手动配置的快照和渲染图。
- `results/batch_runs/` : 10个基准案例的完整验证套件。
- `results/boundary_analysis/` : 动态误差演化和临界阈值曲线。
- `results/mesh_convergence/` : 网格收敛性验证图。

## 📄 许可证
本项目采用MIT许可证 - 详情请见[LICENSE](LICENSE)文件。

