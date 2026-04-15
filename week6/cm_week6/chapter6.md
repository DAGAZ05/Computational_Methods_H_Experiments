# 第六章 数值微分与数值积分
## §6.1 引言
积分学基本定理（Newton-Leibniz公式）：
$$\int_{a}^{b} f(x) d x=F(b)-F(a)$$

实际应用中需使用数值方法的情况：
1. $f(x)$ 的原函数无法用初等函数表示；
2. $f(x)$ 的原函数表达式过于复杂；
3. $f(x)$ 以表格形式给出。

---

## §6.2 数值微分公式
### 重点精讲6.1 数值微分
以离散点 $(x_k, f(x_k))(k=0,1,...,n)$ 近似表达 $y=f(x)$ 在节点 $x_k$ 处的微分，称为**数值微分**。

### 一、Taylor展开法
1. 一阶向前差分公式
$$f(x_k+h)=f(x_k)+hf'(x_k)+\frac{h^2}{2!}f''(\xi_1),\quad \xi_1\in(x_k,x_k+h)$$
$$f'(x_k)=\frac{f(x_k+h)-f(x_k)}{h}-\frac{h}{2}f'(\xi_1)$$

2. 一阶向后差分公式
$$f(x_k-h)=f(x_k)-hf'(x_k)+\frac{h^2}{2!}f''(\xi_2),\quad \xi_2\in(x_k-h,x_k)$$
$$f'(x_k)=\frac{f(x_k)-f(x_k-h)}{h}+\frac{h}{2}f'(\xi_2)$$

3. 中点公式（二阶精度）
$$f'(x_k)=\frac{f(x_k+h)-f(x_k-h)}{2h}-\frac{h^2}{6}f'''(\xi_3),\quad \xi_3\in(x_k-h,x_k+h)$$

4. 二阶导数中点公式
$$f''(x_k)=\frac{f(x_k+h)-2f(x_k)+f(x_k-h)}{h^2}-\frac{h^2}{12}f^{(4)}(\xi_4),\quad \xi_4\in(x_k-h,x_k+h)$$

### 二、插值法
对列表函数建立插值多项式 $p_n(x)$，以 $p_n'(x)$ 近似 $f'(x)$，称为**插值型求导公式**。

余项公式：
$$f'(x_k)-p_n'(x_k)=\frac{f^{(n+1)}(\xi)}{(n+1)!}\omega_{n+1}'(x_k),\quad \xi\in(x_0,x_n)$$

#### 1. 一阶两点公式
$$f'(x)=\frac{f(x+h)-f(x-h)}{2h}+O(h^2)$$

#### 2. 一阶三点公式（等距节点，步长$h$）
$$
\begin{align*}
f'(x_0)&=\frac{-3f(x_0)+4f(x_1)-f(x_2)}{2h}+\frac{h^2}{3}f''(\xi_1)\\
f'(x_1)&=\frac{-f(x_0)+f(x_2)}{2h}-\frac{h^2}{6}f''(\xi_2)\\
f'(x_2)&=\frac{f(x_0)-4f(x_1)+3f(x_2)}{2h}+\frac{h^2}{3}f''(\xi_3)
\end{align*}
$$

### 三、Richardson外推法
### 重点精讲6.2 Richardson外推法
设数值方法近似值 $S^*(h)$ 的截断误差：
$$S-S^*(h)=a_1h^{p_1}+a_2h^{p_2}+\cdots,\quad p_1<p_2<\cdots$$

以 $h/2$ 替换 $h$：
$$S-S^*\left(\frac{h}{2}\right)=a_1\left(\frac{h}{2}\right)^{p_1}+a_2\left(\frac{h}{2}\right)^{p_2}+\cdots$$

线性组合消去首项误差，得到更高精度公式：
$$f'(x)=\frac{4}{3}S^*\left(\frac{h}{2}\right)-\frac{1}{3}S^*(h)+O(h^4)$$

### 数值微分重要备注
1. 步长$h$并非越小精度越高，数值微分对舍入误差敏感，$h$过小会导致计算不稳定；
2. 插值多项式收敛到$f(x)$时，$p_n'(x)$不一定收敛到$f'(x)$；
3. 可使用样条插值函数的导函数避免上述问题。

---

## §6.3 Newton-Cotes求积公式
### 一、插值型求积公式
### 重点精讲6.3 插值型求积公式
用插值多项式 $p(x)$ 近似 $f(x)$，积分近似为：
$$\int_{a}^{b}f(x)dx\approx\sum_{k=0}^nA_kf(x_k)$$
其中：$x_k$ 为求积节点，$A_k$ 为求积系数，$A_k=\int_{a}^b l_k(x)dx$。

截断误差：
$$E[f]=\frac{1}{(n+1)!}\int_{a}^b f^{(n+1)}(\xi)\prod_{i=0}^n(x-x_i)dx$$

### 二、Newton-Cotes求积公式
将 $[a,b]$ $n$ 等分，步长 $h=\frac{b-a}{n}$，节点 $x_k=a+kh$。

求积系数：
$$A_k=(b-a)C_k^{(n)}$$
$C_k^{(n)}$ 为**柯特斯系数**。

#### 常用低阶公式
1. **梯形公式**（$n=1$）
$$\int_{a}^b f(x)dx\approx\frac{b-a}{2}[f(a)+f(b)]$$

2. **Simpson公式（抛物线公式）**（$n=2$）
$$\int_{a}^b f(x)dx\approx\frac{b-a}{6}\left[f(a)+4f\left(\frac{a+b}{2}\right)+f(b)\right]$$

3. **柯特斯公式**（$n=4$）
$$
\begin{align*}
\int_{a}^b f(x)dx&\approx(b-a)\left[\frac{7}{90}f(x_0)+\frac{32}{90}f(x_1)+\frac{12}{90}f(x_2)\right.\\
&\quad\left.+\frac{32}{90}f(x_3)+\frac{7}{90}f(x_4)\right]
\end{align*}
$$

### 三、代数精度
### 重点精讲6.4 代数精度
定义：若求积公式对所有不高于$m$次的多项式精确成立，对$m+1$次多项式不精确成立，则称公式具有**$m$次代数精度**。

判定：对 $f(x)=1,x,x^2,...,x^m$ 精确成立，对 $x^{m+1}$ 不精确成立。

常用公式代数精度：
- 梯形公式：1次
- Simpson公式：3次
- 柯特斯公式：5次

性质：
- Newton-Cotes公式代数精度**不低于$n$次**；
- $n$为偶数时，代数精度为$n+1$次；$n$为奇数时，为$n$次。

### 四、Newton-Cotes公式截断误差
### 重点精讲6.5 Newton-Cotes公式截断误差
1. **梯形公式误差**（$f(x)$二阶连续可导）
$$E_T[f]=-\frac{(b-a)^3}{12}f''(\eta),\quad a\leq\eta\leq b$$

2. **Simpson公式误差**（$f(x)$四阶连续可导）
$$E_S[f]=-\frac{(b-a)^5}{2880}f^{(4)}(\eta),\quad a\leq\eta\leq b$$

3. **柯特斯公式误差**（$f(x)$六阶连续可导）
$$E_C[f]=-\frac{2(b-a)}{945}\left(\frac{b-a}{4}\right)^6f^{(6)}(\eta),\quad a\leq\eta\leq b$$

### 五、稳定性与收敛性
1. **收敛性**：$\lim\limits_{n\to\infty}\sum_{k=0}^nA_kf(x_k)=\int_a^bf(x)dx$；
2. **稳定性**：求积系数$A_k>0$时，公式稳定；
3. $n\leq7$时柯特斯系数全正，$n\geq8$出现负系数，稳定性变差。

---

## §6.4 复化求积法
### 重点精讲6.6 复化求积法
将$[a,b]$等分为$n$个子区间，步长$h=\frac{b-a}{n}$，在每个子区间用低阶公式求和。

### 一、复化梯形公式
$$
\begin{align*}
T_n&=\frac{h}{2}\left[f(a)+2\sum_{k=1}^{n-1}f(x_k)+f(b)\right]\\
E_{T_n}&=-\frac{(b-a)}{12}h^2f''(\eta),\quad \eta\in[a,b]
\end{align*}
$$

误差估计：
$$|E_{T_n}|\leq\frac{b-a}{12}h^2M_2,\quad M_2=\max|f''(x)|$$

### 二、复化Simpson公式
$$
\begin{align*}
S_n&=\frac{h}{6}\left[f(a)+4\sum_{k=0}^{n-1}f(x_{k+\frac{1}{2}})+2\sum_{k=1}^{n-1}f(x_k)+f(b)\right]\\
E_{S_n}&=-\frac{b-a}{2880}h^4f^{(4)}(\eta),\quad \eta\in[a,b]
\end{align*}
$$

### 三、区间逐次分半求积法
### 重点精讲6.7 区间逐次分半求积法
1. **事后估计误差**：无需先验估计导数，通过区间分半，判断相邻两次结果差值是否满足精度。
2. **递推公式**
$$T_{2n}=\frac{1}{2}T_n+\frac{h}{2}\sum_{k=1}^nf\left(a+(2k-1)\frac{h}{2}\right)$$
3. **外推加速**
$$
\begin{align*}
I&\approx T_{2n}+\frac{1}{3}(T_{2n}-T_n)\\
I&\approx S_{2n}+\frac{1}{15}(S_{2n}-S_n)\\
I&\approx C_{2n}+\frac{1}{63}(C_{2n}-C_n)
\end{align*}
$$

---

## §6.5 Romberg求积法
### 一、外推构造高精度公式
由低精度序列组合得到高精度序列：
$$
\begin{align*}
S_n&=\frac{4}{3}T_{2n}-\frac{1}{3}T_n\\
C_n&=\frac{16}{15}S_{2n}-\frac{1}{15}S_n\\
R_n&=\frac{64}{63}C_{2n}-\frac{1}{63}C_n
\end{align*}
$$

- $T_n$：梯形序列
- $S_n$：Simpson序列
- $C_n$：柯特斯序列
- $R_n$：Romberg序列

每次外推误差阶**提高二阶**。

### 二、Romberg算法实现
构造T数表，逐次计算梯形、Simpson、柯特斯、Romberg序列，直到Romberg序列相邻项差值满足误差限。

优点：内存占用少、精度高，工程常用。

---

## §6.6 Gauss型求积公式
### 一、Gauss型求积公式定义
$n+1$个节点的插值型求积公式，具有**$2n+1$次代数精度**，称为**Gauss型求积公式**，节点为**Gauss点**。

Gauss点充要条件：以Gauss点为零点的$n+1$次多项式 $\omega_{n+1}(x)=\prod\limits_{j=0}^n(x-x_j)$ 与任意次数$\leq n$的多项式$P(x)$正交：
$$\int_a^b\omega_{n+1}(x)P(x)dx=0$$

### 二、Gauss-Legendre求积公式
区间$[-1,1]$上，以**Legendre多项式**零点为Gauss点。

Legendre多项式：
$$P_{n+1}(x)=\frac{1}{(n+1)!2^{n+1}}\frac{d^{n+1}}{dx^{n+1}}(x^2-1)^{n+1}$$

#### 常用公式
1. 两点公式
$$\int_{-1}^1f(x)dx\approx f\left(-\frac{\sqrt{3}}{3}\right)+f\left(\frac{\sqrt{3}}{3}\right)$$

2. 三点公式
$$\int_{-1}^1f(x)dx\approx\frac{5}{9}f(-\sqrt{0.6})+\frac{8}{9}f(0)+\frac{5}{9}f(\sqrt{0.6})$$

任意区间$[a,b]$变换：
$$x=\frac{a+b}{2}+\frac{b-a}{2}t,\quad \int_a^bf(x)dx=\frac{b-a}{2}\int_{-1}^1f\left(\frac{a+b}{2}+\frac{b-a}{2}t\right)dt$$

### 三、截断误差与稳定性
截断误差：
$$R[f]=\frac{f^{(2n+2)}(\xi)}{(2n+2)!}\int_a^b\omega_{n+1}^2(x)dx,\quad \xi\in(a,b)$$

优点：代数精度最高、收敛、稳定；缺点：节点增加时，旧节点函数值无法复用，可使用**复化Gauss求积**。

---

## 数值算例
1. 计算 $\int_0^1\frac{4}{1+x^2}dx$；
2. 用复化梯形、复化Simpson公式计算 $\int_2^8\frac{1}{x}dx$ 求$\ln2$，确定满足有效数字的最少节点数。
