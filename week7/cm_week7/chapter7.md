# 第七章 常微分方程初值问题的数值解法
## §7.1 引言
本章着重讨论一阶常微分方程初值问题
$
\begin{cases}
y'=\dfrac{dy}{dx}=f(x, y) \\
y(a)=y_0
\end{cases}
$
在区间 $[a, b]$ 上的数值解法。

这类问题多数无法求出解析解，只能用近似方法求解，近似方法分为两类：
1. **近似解析法**：级数解法、逐次逼近法等；
2. **数值解法**：给出解在一些离散点上的近似值。

### 初值问题解的存在唯一性
设 $f(x, y)$ 在区域 $D=\{a \leq x \leq b, y \in \mathbb{R}\}$ 上连续，且关于 $y$ 满足**李普希兹（Lipschitz）条件**：存在常数 $L$，使得对 $D$ 内任意两个 $y_1,y_2$，有
$
|f(x, y_1)-f(x, y_2)| \leq L|y_1-y_2|
$
则初值问题**存在唯一解**，且解连续可微。

**简便验证方法**：若 $f(x, y)$ 对 $y$ 可微，且偏导数 $\dfrac{\partial f}{\partial y}$ 有界，则Lipschitz条件成立，其中
$
L=\max\left|\dfrac{\partial f(x, y)}{\partial y}\right|
$

### 数值解法基本概念
- 离散点：$x_i = x_{i-1}+h_i$，$h$ 为**步长**，通常取等步长；
- 符号：$y(x_n)$ 为解析解，$y_n$ 为数值解，$f_n=f(x_n,y_n)$，$y'(x_n)=f(x_n,y(x_n))$；
- 求解方式：**逐步计算**，由 $y_n$ 计算 $y_{n+1}$。

### 方法分类
1. **单步法**：计算 $y_{n+1}$ 仅用到 $x_{n+1},x_n,y_n$（前一步值）；
2. **多步法**：计算 $y_{n+1}$ 需用到 $x_{n-p},y_{n-p}(p=1,2,\dots,k)$（前 $k$ 步值）。

单步法与多步法均分为**显式**和**隐式**：
- 显式：$y_{n+1}$ 可直接由已知值算出；
- 隐式：$y_{n+1}$ 含于方程中，需迭代求解。

---

## §7.2 欧拉法与梯形法
### 一、显式Euler公式
#### 公式形式
等步长 $x_n=x_0+nh$，显式Euler公式：
$
y_{n+1}=y_n+hf(x_n,y_n) \quad (n=0,1,2,\dots)
$
这是**最简单的显式单步法**。

#### 推导方法
1. **差商法**：用向前差商 $\dfrac{y_{n+1}-y_n}{h}$ 近似导数 $y'(x_n)$；
2. **Taylor展开法**：取 $y(x_{n+1})$ Taylor展开的线性部分；
3. **数值积分法**：对 $y'=f(x,y)$ 从 $x_n$ 到 $x_{n+1}$ 积分，用**左矩形公式**近似积分项。

#### 几何意义
Euler折线法：从初始点 $(x_0,y_0)$ 出发，以 $f(x_n,y_n)$ 为斜率作直线段，依次得到折线 $P_0P_1P_2\cdots$ 近似解曲线 $y(x)$。

### 二、隐式Euler公式（后退Euler公式）
#### 公式形式
$
y_{n+1}=y_n+hf(x_{n+1},y_{n+1}) \quad (n=0,1,2,\dots)
$
属于**隐式单步法**，需迭代求解。

#### 迭代求解格式
$
\begin{cases}
y_{n+1}^{(0)}=y_n+hf(x_n,y_n) \\
y_{n+1}^{(s+1)}=y_n+hf(x_{n+1},y_{n+1}^{(s)})
\end{cases} \quad s=0,1,2,\dots
$
当 $|y_{n+1}^{(s+1)}-y_{n+1}^{(s)}|<\varepsilon$（误差限）时停止迭代。

#### 收敛条件
$f(x,y)$ 满足Lipschitz条件，且 $hL<1$ 时迭代收敛。

### 三、梯形公式
#### 公式形式
对积分项用**梯形求积公式**近似，得：
$
y_{n+1}=y_n+\frac{h}{2}\left[f(x_n,y_n)+f(x_{n+1},y_{n+1})\right]
$
**隐式单步法**。

#### 迭代求解格式
$
\begin{cases}
y_{n+1}^{(0)}=y_n+hf(x_n,y_n) \\
y_{n+1}^{(s+1)}=y_n+\dfrac{h}{2}\left[f(x_n,y_n)+f(x_{n+1},y_{n+1}^{(s)})\right]
\end{cases} \quad s=0,1,2,\dots
$

#### 收敛条件
$\dfrac{hL}{2}<1$ 时迭代收敛。

### 四、Euler-梯形预估校正公式（改进Euler法）
实用中仅迭代一次，得到**显式格式**：
$
\begin{cases}
\text{预估：} y_{n+1}^{(0)}=y_n+hf(x_n,y_n) \\
\text{校正：} y_{n+1}=y_n+\dfrac{h}{2}\left[f(x_n,y_n)+f(x_{n+1},y_{n+1}^{(0)})\right]
\end{cases}
$

**K形式**：
$
\begin{cases}
y_{n+1}=y_n+\dfrac{h}{2}(K_1+K_2) \\
K_1=f(x_n,y_n) \\
K_2=f(x_n+h,y_n+hK_1)
\end{cases}
$

---

## 重点精讲7.3 单步法的局部截断误差和阶
### 基本定义
设单步法通式：
- 显式：$y_{n+1}=y_n+h\phi(x_n,y_n,h)$
- 隐式：$y_{n+1}=y_n+h\phi(x_n,y_n,y_{n+1},h)$

1. **整体截断误差**：$e_n=y(x_n)-y_n$（误差累计）；
2. **局部截断误差**：假设 $y_n=y(x_n)$，则 $R_{n+1}=y(x_{n+1})-y_{n+1}$；
3. **方法阶数**：若 $R_{n+1}=O(h^{p+1})$，称方法为**$p$ 阶方法**；
4. **主局部截断误差**：$R_{n+1}$ 按 $h$ 展开的首项 $\Psi(x_n,y(x_n))h^{p+1}$。

### 典型方法误差与阶数
1. **显式Euler**：$R_{n+1}=\dfrac{h^2}{2}y''(x_n)+O(h^3)$，**1阶**；
2. **隐式Euler**：$R_{n+1}=-\dfrac{h^2}{2}y''(x_n)+O(h^3)$，**1阶**；
3. **梯形公式**：$R_{n+1}=-\dfrac{h^3}{12}y'''(x_n)+O(h^4)$，**2阶**；
4. **预估校正法**：$R_{n+1}=O(h^3)$，**2阶**。

---

## §7.3 泰勒展开法与龙格-库塔（Runge-Kutta）方法
### 一、Taylor级数展开法
设解 $y(x)$ 与 $f(x,y)$ 充分光滑，$y(x)$ 在 $x_n$ 处Taylor展开：
$
y(x_{n+1})=y(x_n)+hy'(x_n)+\frac{h^2}{2!}y''(x_n)+\cdots+\frac{h^p}{p!}y^{(p)}(x_n)+O(h^{p+1})
$
略去余项得**$p$ 阶Taylor方法**。

- $p=1$：显式Euler公式；
- $p=2$：二阶Taylor公式
$
y_{n+1}=y_n+hf_n+\frac{h^2}{2}\left(f_x+ff_y\right)\big|_{(x_n,y_n)}
$

**缺陷**：需计算高阶导数，计算量大，极少单独使用。

### 二、龙格-库塔（R-K）方法
#### 基本思想
用**不同节点的函数值线性组合**构造高阶单步格式，避免高阶导数计算。

#### 一般形式
$
\begin{cases}
y_{n+1}=y_n+h\phi(x_n,y_n,h) \\
\phi(x_n,y_n,h)=\displaystyle\sum_{i=1}^r c_iK_i \\
K_1=f(x_n,y_n) \\
K_i=f\left(x_n+\alpha_ih,y_n+h\displaystyle\sum_{j=1}^{i-1}\beta_{ij}K_j\right),\ i=2,3,\dots,r
\end{cases}
$
$r$ 为级数，通过系数匹配使局部截断误差阶数最高。

#### 二级二阶R-K公式（无穷多组）
通式：
$
\begin{cases}
y_{n+1}=y_n+h(c_1K_1+c_2K_2) \\
K_1=f(x_n,y_n) \\
K_2=f(x_n+\alpha_2h,y_n+h\beta_{21}K_1)
\end{cases}
$
系数满足：
$
\begin{cases}
c_1+c_2=1 \\
c_2\alpha_2=\dfrac{1}{2} \\
c_2\beta_{21}=\dfrac{1}{2}
\end{cases}
$

**常用二级二阶格式**：
1. **预估校正法**：$c_1=c_2=\dfrac{1}{2},\alpha_2=\beta_{21}=1$；
2. **中点公式**：$c_1=0,c_2=1,\alpha_2=\beta_{21}=\dfrac{1}{2}$；
3. **二阶Heun方法**：$c_1=\dfrac{1}{4},c_2=\dfrac{3}{4},\alpha_2=\beta_{21}=\dfrac{2}{3}$。

#### 高阶经典R-K公式
1. **三阶Kutta公式**
$
y_{n+1}=y_n+\frac{h}{6}(K_1+4K_2+K_3)
$
2. **三阶Heun公式**
$
y_{n+1}=y_n+\frac{h}{4}(K_1+3K_3)
$
3. **标准四阶R-K公式（最常用）**
$
\begin{cases}
y_{n+1}=y_n+\dfrac{h}{6}(K_1+2K_2+2K_3+K_4) \\
K_1=f(x_n,y_n) \\
K_2=f\left(x_n+\dfrac{h}{2},y_n+\dfrac{h}{2}K_1\right) \\
K_3=f\left(x_n+\dfrac{h}{2},y_n+\dfrac{h}{2}K_2\right) \\
K_4=f(x_n+h,y_n+hK_3)
\end{cases}
$

#### R-K方法重要结论
1. 二级R-K最多**2阶**，三级最多**3阶**，四级最多**4阶**；
2. 四阶R-K精度高、计算量适中，工程最常用；
3. 更高阶R-K计算量激增，一般不采用。

### 三、单步法的收敛性与稳定性
#### 收敛性
若对固定 $x=x_0+nh$，$\lim\limits_{h \to 0}y_n=y(x)$，则方法**收敛**。
**定理**：局部截断误差 $O(h^{p+1})$ 且 $\phi$ 对 $y$ 满足Lipschitz条件 $\implies$ 方法收敛。

#### 稳定性（绝对稳定）
针对模型方程 $y'=\lambda y\ (\text{Re}(\lambda)<0)$，若计算误差不随步数放大，称**绝对稳定**。

**常用方法绝对稳定区间**：
1. 显式Euler：$h\lambda \in (-2,0)$；
2. 隐式Euler：$h\lambda \in (-\infty,0)$（无条件稳定）；
3. 梯形公式：$\text{Re}(h\lambda)<0$（无条件稳定）；
4. 四阶经典R-K：$h\lambda \in (-2.785,0)$。

---

## §7.4 线性多步法
### 基本思想
充分利用**前多步信息**提高精度，不显著增加计算量。构造方法：**数值积分法**、**Taylor展开法**。

### 线性多步法一般形式
$
y_{n+1}=\sum_{i=0}^r \alpha_i y_{n-i}+h\sum_{i=-1}^r \beta_i f_{n-i}
$
- $\beta_{-1}=0$：**显式**；$\beta_{-1} \neq 0$：**隐式**；
- 用到 $y_n,y_{n-1},\dots,y_{n-r}$，称为**线性 $r+1$ 步法**。

### 一、Adams方法
#### 1. Adams外插公式（显式四步）
$
y_{n+1}=y_n+\frac{h}{24}\left[55f_n-59f_{n-1}+37f_{n-2}-9f_{n-3}\right]
$
局部截断误差：$\dfrac{251}{720}h^5y^{(5)}(\eta)$，**4阶**。

#### 2. Adams内插公式（隐式三步）
$
y_{n+1}=y_n+\frac{h}{24}\left[9f_{n+1}+19f_n-5f_{n-1}+f_{n-2}\right]
$
**4阶隐式方法**。

#### 3. Adams预估-校正格式
$
\begin{cases}
\text{预估：} y_{n+1}^{(0)}=y_n+\dfrac{h}{24}\left[55f_n-59f_{n-1}+37f_{n-2}-9f_{n-3}\right] \\
\text{校正：} y_{n+1}=y_n+\dfrac{h}{24}\left[9f(x_{n+1},y_{n+1}^{(0)})+19f_n-5f_{n-1}+f_{n-2}\right]
\end{cases}
$

### 二、Taylor展开法构造线性多步法
以**两步法**为例，令局部截断误差首项系数为0，解得：
$
y_{n+1}=y_{n-1}+\frac{h}{3}\left(f_{n+1}+4f_n+f_{n-1}\right)
$
即**Simpson公式**，局部截断误差 $-\dfrac{1}{90}h^5y^{(5)}(x_n)$，**4阶**。

### 三、出发值计算
线性 $k$ 步法需要 $k$ 个初始值 $y_0,y_1,\dots,y_{k-1}$，仅 $y_0$ 已知，其余用**同阶单步法（如四阶R-K）**计算，保证精度匹配。

---

## §7.5 数值算例
**问题**：求解初值问题
$
\begin{cases}
\dfrac{dy}{dx}=y-\dfrac{2x}{y} \\
y(0)=1
\end{cases}
$
区间 $[0,1]$，步长 $h=0.1$，**理论解**：$y(x)=\sqrt{1+2x}$。

**方法对比**：显式Euler、Adams外插公式与理论解的计算结果表明，高阶方法精度显著优于低阶方法。

