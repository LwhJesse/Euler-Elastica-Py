### 🇰🇷 Korean (`README.ko.md`)

# 오일러 탄성 곡선 & 비선형 보 해석기
[English](README.md) | [Français](README.fr.md) |  [Español](README.es.md) |[繁體中文](README.zh-TW.md) | [简体中文](README.zh-CN.md) | [Deutsch](README.de.md) | [日本語](README.ja.md) | [한국어](README.ko.md)
<div align="center">
  <img src="results/batch_runs/3d_renders/3D_Stress_FEM_Case8_F-20.0N.jpg" width="800" alt="비선형 보의 3D 등각 투영 렌더">
  <p><em>극심한 중앙점 하중을 받는 단순 지지 보의 대변형 (Corotational FEM).</em></p>
</div>

**대변형 기하학적 비선형 보(오일러 탄성 곡선)** 문제를 해결, 시뮬레이션 및 검증하기 위해 설계된 매우 견고한 Python 기반 계산 역학 프레임워크입니다.

세 가지 다른 수학적 차원(해석적, 룽게-쿠타 슈팅 방법, Corotational FEM)을 교차 검증함으로써 이 도구는 선형 가정과 비선형 현실 사이의 정확한 수치적 경계를 명확하게 매핑합니다.

## 🧮 이론적 기반

전통적인 선형 역학(오일러-베르누이 보 이론)은 미소 회전을 가정하여($w' \approx 0$) 정확한 곡률 방정식을 단순화하며, 이는 극심한 하중 하에서 처짐을 심각하게 과대평가하게 만듭니다. 이 프레임워크는 완전한 기하학적 비선형성을 유지함으로써 **오일러 탄성 곡선** 문제를 근본적으로 해결합니다:

$$ \kappa = \frac{M(x)}{EI} = \frac{w''}{(1 + (w')^2)^{3/2}} $$

이 고도로 결합된 비선형 미분 방정식은 대부분의 복잡한 하중 사례에 대해 닫힌 형태의 해석해가 없으므로, 이 프로젝트에서 구현된 고급 수치 해석기가 필요합니다.

## ⚙️ 기술적 세부사항 및 알고리즘

이 프레임워크는 세 가지 독립적인 해석 엔진과 고급 공간 매핑을 통해 높은 정확도의 검증을 달성합니다:

### 1. 3중 수학적 교차 검증 및 라그랑주 좌표 정렬
* **선형 기준선:** 고전적인 오일러-베르누이 미소 처짐 이론에 기반한 닫힌 형태의 해. 이는 비선형 발산의 임계점을 정량화하기 위한 비교 기준으로 사용됩니다.
* **정확한 비선형 엔진:** 오일러 탄성 곡선의 정확한 미분 방정식에 대한 룽게-쿠타 적분.
* **FEM 골드 스탠다드:** `Corotational` 기하학적 비선형 공식을 활용하는 OpenSeesPy.
* **원자 대 원자 라그랑주 정렬:** 극심한 굽힘으로 인해 보의 수평 투영이 급격히 줄어듭니다. 모든 오차 평가와 시각화는 엄격하게 **라그랑주 좌표계**(초기 호 길이 $s$ 기반, `scipy.interpolate.interp1d` 사용)에서 수행됩니다. 이는 오일러 좌표계에서 극심한 기하학적 변형에 내재된 공간적 불일치 문제를 효과적으로 해결합니다.

### 2. RK 슈팅 방법 엔진
* **상태 공간 공식화:** 고차 미분 방정식을 1차 상미분 방정식(ODE) 시스템으로 변환하여 상태 벡터 $\mathbf{y} =[w, \theta, M, V, N]^T$를 설정합니다.
* **경계값 문제(BVP)를 초기값 문제(IVP)로 변환:** 슈팅 방법을 통해 경계값 문제(BVP)를 해결합니다. **레벤버그-마쿼트(`lm`) 알고리즘**(`scipy.optimize.least_squares`)을 사용하여 초기 추측값을 반복적으로 개선함으로써 깊은 굽힘 영역에서 자코비 행렬의 특이성을 효과적으로 방지합니다.
* **불연속성 처리:** 내부 힘/모멘트의 점프를 부드럽게 해결하기 위해 조각별 연속 경계 일치 조건을 구현합니다.

### 3. 자동화된 경계 탐색 엔진
* **근 찾기 알고리즘:** **브렌트 방법**(`scipy.optimize.brentq`)을 통합하여 선형 해석 모델과 비선형 FEM 모델 간의 상대 오차가 엄격한 **5% 임계값**에 도달하는 정확한 적용 하중(또는 하중 위치)을 동적으로 찾아냅니다.

### 4. 제로 종속성 고화질 3D 렌더링
* VTK, Mayavi, ParaView와 같은 무거운 3D 과학 시각화 라이브러리를 사용하지 않습니다. 1D 보 요소를 수학적으로 3D 물리적 고체로 재구성하고 등가 응력 텐서를 표면에 매핑하여, **순수 Matplotlib**만을 사용하여 출판 가능한 3D 등각 투영 응력 등고선을 생성합니다.

## 📦 환경 설정

이 프로젝트는 계산용 Python 스크립트 모음입니다. 종속성이 충족되면 바로 실행할 수 있습니다.

### 옵션 A: 표준 설정 (Windows / macOS / Debian 기반 / Red Hat 기반)
표준 OS 환경에서는 가상 환경 사용이 선택 사항이지만 종속성 충돌을 피하기 위해 권장됩니다.
```bash
# 1. 리포지토리 복제
git clone https://github.com/LwhJesse/Euler-Elastica-Py.git
cd Euler-Elastica-Py

# 2. (선택 사항) 가상 환경 생성 및 활성화
python3 -m venv venv
source venv/bin/activate  # Windows에서는 다음을 사용: venv\Scripts\activate

# 3. 핵심 종속성 설치
pip install --upgrade pip
pip install -r requirements.txt
```
> **⚠️ macOS Apple Silicon (M1/M2/M3) 사용자 참고:** `openseespy`는 C++로 래핑된 프레임워크입니다. `pip`이 호환되는 ARM64 휠을 찾지 못하면, Rosetta 2를 사용하여 터미널을 실행하거나 [OpenSeesPy 공식 설치 문서](https://openseespydoc.readthedocs.io/en/latest/src/installation.html)를 참조해야 할 수 있습니다.

### 옵션 B: Arch 기반 Linux (AUR)
Arch 기반 Linux를 사용 중인 경우, 전역 `pip` 설치는 외부에서 관리됩니다(PEP 668). 가상 환경을 안전하게 건너뛰고 시스템 패키지 관리자를 통해 직접 종속성을 설치할 수 있습니다.
*(참고: `python-openseespy` AUR 패키지는 이 리포지토리의 저자 [@LwhJesse](https://aur.archlinux.org/packages/python-openseespy)가 공식적으로 관리합니다).*
```bash
# 1. 리포지토리 복제
git clone https://github.com/LwhJesse/Euler-Elastica-Py.git
cd Euler-Elastica-Py

# 2. pacman과 AUR 헬퍼(예: yay)를 통해 종속성 설치
sudo pacman -S python-numpy python-pandas python-scipy python-matplotlib python-rich
yay -S python-openseespy
```

## 🚀 사용 가이드

이 프로젝트는 `core/config.py`에 위치한 단일 파라미터 구성 소스를 따릅니다.

### 1. 단일 사례 분석
단일 사례를 실행하여 2D 라그랑주 비교 플롯을 생성합니다:
```bash
python run_single.py
```
현재 구성에 대한 3D 등각 투영 응력 등고선을 생성합니다:
```bash
python run_3d_render.py
```

### 2. 일괄 실행 및 벤치마킹
10개의 사전 정의된 벤치마크 사례를 모두 실행하여 포괄적인 검증 스위트를 생성합니다:
```bash
python run_batch.py
python run_batch_3d.py
```
<details>
<summary><b>클릭하여 보기: 정확한 RK 대 FEM 다중 물리 정렬 (사례 8)</b></summary>
<br>
<div align="center">
  <img src="results/batch_runs/2d_plots/case8-1_comparison.jpg" width="800" alt="2D 비교 플롯">
  <p><em>선형 해석해(녹색)가 급격히 벗어나는 반면, RK(빨간색)와 FEM(검은색)은 완벽하게 일치하는 것을 확인하세요.</em></p>
</div>
</details>

### 3. 임계 경계 탐색 (라이브 클러스터)
`rich` 기반의 라이브 터미널 대시보드를 사용하여 모든 CPU 코어에서 동시에 5% 비선형 오차 상 경계를 매핑합니다:
```bash
python run_multiprocess.py
```
<details>
<summary><b>클릭하여 보기: 이변수 안전/위험 구역 맵 (사례 3)</b></summary>
<br>
<div align="center">
  <img src="results/boundary_analysis/pure_critical_boundary_fem_case3.jpg" width="600" alt="경계 포락선 다이어그램">
  <p><em>빨간색 임계 경계는 기하학적 비선형 모델을 채택해야 할 시점을 명시적으로 나타냅니다.</em></p>
</div>
</details>

## 🛠️ 지원되는 벤치마크 하중 사례

| 사례 ID | 경계 조건 | 하중 유형 | 활성 변수 |
| :---: | :--- | :--- | :--- |
| **1** | 캔틸레버 | 끝단 굽힘 모멘트 | $M_e$ |
| **2** | 캔틸레버 | 끝단 집중 하중 | $F$ |
| **3** | 캔틸레버 | 중간 집중 하중 | $F, a$ |
| **4** | 캔틸레버 | 등분포 하중 | $q$ |
| **5** | 단순 지지 | 왼쪽 끝 모멘트 | $M_e$ |
| **6** | 단순 지지 | 오른쪽 끝 모멘트 | $M_e$ |
| **7** | 단순 지지 | 중간 모멘트 | $M_e, a$ |
| **8** | 단순 지지 | 중앙점 집중 하중 | $F$ |
| **9** | 단순 지지 | 중간 집중 하중 | $F, a$ |
| **10** | 단순 지지 | 등분포 하중 | $q$ |

## 📁 출력 디렉터리 구조
결과는 소스 코드와 완벽하게 분리되어 자동으로 정리됩니다:
- `results/single_runs/` : 수동 구성의 스냅샷 및 렌더링.
- `results/batch_runs/` : 10개 벤치마크 사례의 전체 검증 스위트.
- `results/boundary_analysis/` : 동적 오차 변화 및 임계 임계값 곡선.
- `results/mesh_convergence/` : 메쉬 독립성 검증 플롯.

## 📄 라이선스
이 프로젝트는 MIT 라이선스에 따라 라이선스가 부여됩니다 - 자세한 내용은 [LICENSE](LICENSE) 파일을 참조하십시오.