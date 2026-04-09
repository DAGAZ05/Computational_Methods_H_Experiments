import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator

# ============================================================
# SCI 学术风格全局设置
# ============================================================
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman'],
    'mathtext.fontset': 'stix',
    'font.size': 10,
    'axes.linewidth': 0.8,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'xtick.minor.width': 0.4,
    'ytick.minor.width': 0.4,
    'xtick.major.size': 4,
    'ytick.major.size': 4,
    'xtick.minor.size': 2,
    'ytick.minor.size': 2,
    'legend.fontsize': 9,
    'legend.frameon': True,
    'legend.edgecolor': 'black',
    'legend.fancybox': False,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

# ============================================================
# 数据
# ============================================================
x = np.array([0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0])
y = np.array([0.302, 0.185, 0.106, 0.093, 0.240, 0.579, 0.561, 0.468, 0.302])
n = len(x) - 1  # 区间数 = 8

# 步长 h_i = x_i - x_{i-1}, i = 1,...,n
h = np.diff(x)

# 一阶差商 f[x_{i-1}, x_i]
f1 = np.diff(y) / h

# 右端项 d_i = 6 * f[x_{i-1}, x_i, x_{i+1}], i = 1,...,n-1
d = np.array([6 * (f1[i] - f1[i-1]) / (x[i+1] - x[i-1]) for i in range(1, n)])

# λ_i, μ_i
lam = np.array([h[i-1] / (h[i-1] + h[i]) for i in range(1, n)])
mu  = np.array([h[i]   / (h[i-1] + h[i]) for i in range(1, n)])

# ============================================================
# 边界条件名称
# ============================================================
bc_names = {
    1: "BC1: Clamped (slope = finite diff)",
    2: "BC2: Natural (M0=Mn=0)",
    3: "BC3: M0=M1, Mn=Mn-1",
    4: "BC4: Not-a-knot"
}

# ============================================================
# 求解三弯矩方程组
# ============================================================
def solve_spline(bc_type):
    A = np.zeros((n+1, n+1))
    b_vec = np.zeros(n+1)

    # 内部节点方程 i = 1,...,n-1
    for i in range(1, n):
        A[i, i-1] = lam[i-1]
        A[i, i]   = 2.0
        A[i, i+1] = mu[i-1]
        b_vec[i]  = d[i-1]

    if bc_type == 1:
        # 夹持边界: S'(x0)=f[x0,x1], S'(xn)=f[xn-1,xn]
        # => 2M0 + M1 = 6/h1*(f[x0,x1] - f'(x0)) = 0
        # => Mn-1 + 2Mn = 6/hn*(f'(xn) - f[xn-1,xn]) = 0
        A[0, 0] = 2.0;  A[0, 1] = 1.0;  b_vec[0] = 0.0
        A[n, n-1] = 1.0; A[n, n] = 2.0;  b_vec[n] = 0.0
    elif bc_type == 2:
        # 自然边界: M0 = 0, Mn = 0
        A[0, 0] = 1.0; b_vec[0] = 0.0
        A[n, n] = 1.0; b_vec[n] = 0.0
    elif bc_type == 3:
        # M0 = M1, Mn = Mn-1
        A[0, 0] = 1.0; A[0, 1] = -1.0; b_vec[0] = 0.0
        A[n, n] = 1.0; A[n, n-1] = -1.0; b_vec[n] = 0.0
    elif bc_type == 4:
        # (M1-M0)/h1 = (M2-M1)/h2  =>  -h2*M0 + (h1+h2)*M1 - h1*M2 = 0
        A[0, 0] = -h[1]; A[0, 1] = h[0]+h[1]; A[0, 2] = -h[0]; b_vec[0] = 0.0
        # (Mn-Mn-1)/hn = (Mn-1-Mn-2)/hn-1
        A[n, n-2] = -h[-2]; A[n, n-1] = h[-2]+h[-1]; A[n, n] = -h[-1]; b_vec[n] = 0.0

    M = np.linalg.solve(A, b_vec)

    # 样条系数: S_i(x) = a_i + b_i*(x-x_i) + c_i*(x-x_i)^2 + d_i*(x-x_i)^3
    a_c = y[:-1].copy()
    b_c = (y[1:] - y[:-1]) / h - h * (2*M[:-1] + M[1:]) / 6
    c_c = M[:-1] / 2
    d_c = (M[1:] - M[:-1]) / (6 * h)

    return M, a_c, b_c, c_c, d_c

def eval_spline(xv, a_c, b_c, c_c, d_c):
    """对标量或数组 xv 求值"""
    xv = np.atleast_1d(xv)
    result = np.empty_like(xv)
    for k, xval in enumerate(xv):
        i = np.searchsorted(x, xval) - 1
        i = np.clip(i, 0, n-1)
        dx = xval - x[i]
        result[k] = a_c[i] + b_c[i]*dx + c_c[i]*dx**2 + d_c[i]*dx**3
    return result

# ============================================================
# 求解四种边界条件
# ============================================================
results = {}
for bc in range(1, 5):
    M, a_c, b_c, c_c, d_c = solve_spline(bc)
    results[bc] = {'M': M, 'a': a_c, 'b': b_c, 'c': c_c, 'd': d_c}

# ============================================================
# 1. 输出每种边界条件下的三次样条函数（分段多项式系数）
# ============================================================
print("=" * 80)
print("Cubic Spline Functions: S_i(x) = a_i + b_i*(x-x_i) + c_i*(x-x_i)^2 + d_i*(x-x_i)^3")
print("=" * 80)
for bc in range(1, 5):
    r = results[bc]
    print(f"\n--- {bc_names[bc]} ---")
    print(f"{'Interval':<16} {'a_i':>12} {'b_i':>12} {'c_i':>12} {'d_i':>12}")
    print("-" * 68)
    for i in range(n):
        print(f"[{x[i]:.1f}, {x[i+1]:.1f}]       {r['a'][i]:12.6f} {r['b'][i]:12.6f} {r['c'][i]:12.6f} {r['d'][i]:12.6f}")

# ============================================================
# 2. 输出弯矩 M_i
# ============================================================
print("\n" + "=" * 80)
print("Moments M_i (second derivatives at nodes)")
print("=" * 80)
header = f"{'Node':>6}"
for bc in range(1, 5):
    header += f"  {'BC'+str(bc):>12}"
print(header)
print("-" * 60)
for i in range(n+1):
    row = f"x_{i} = {x[i]:.1f}"
    for bc in range(1, 5):
        row += f"  {results[bc]['M'][i]:12.6f}"
    print(row)

# ============================================================
# 3. 输出 x=0.05 处的样条值
# ============================================================
x_target = 0.05
print("\n" + "=" * 80)
print(f"Spline values at x = {x_target}")
print("=" * 80)
for bc in range(1, 5):
    r = results[bc]
    val = eval_spline(x_target, r['a'], r['b'], r['c'], r['d'])[0]
    print(f"  {bc_names[bc]}: S({x_target}) = {val:.6f}")

# ============================================================
# 4. 在 [0.2, 0.8] 区间上的比较 + 均方误差 & 最大偏差
# ============================================================
# 以 BC2（自然样条）为基准，计算其他 BC 相对于 BC2 的偏差
# 同时计算所有 BC 两两之间的偏差
x_full = np.linspace(x[0], x[-1], 500)
x_inner = np.linspace(0.2, 0.8, 300)

spline_full = {}
spline_inner = {}
for bc in range(1, 5):
    r = results[bc]
    spline_full[bc] = eval_spline(x_full, r['a'], r['b'], r['c'], r['d'])
    spline_inner[bc] = eval_spline(x_inner, r['a'], r['b'], r['c'], r['d'])

# 以四种 BC 的均值作为参考基准
ref_inner = np.mean([spline_inner[bc] for bc in range(1, 5)], axis=0)
ref_full  = np.mean([spline_full[bc]  for bc in range(1, 5)], axis=0)

print("\n" + "=" * 80)
print("Deviation analysis on [0.2, 0.8] (reference = mean of all 4 BCs)")
print("=" * 80)
print(f"{'BC':>6}  {'RMSD':>12}  {'Max |dev|':>12}  {'Mean dev':>12}")
print("-" * 50)
for bc in range(1, 5):
    diff = spline_inner[bc] - ref_inner
    rmsd = np.sqrt(np.mean(diff**2))
    maxd = np.max(np.abs(diff))
    meand = np.mean(diff)
    print(f"{'BC'+str(bc):>6}  {rmsd:12.6e}  {maxd:12.6e}  {meand:12.6e}")

# 两两之间的偏差
print("\n" + "=" * 80)
print("Pairwise comparison on [0.2, 0.8]")
print("=" * 80)
print(f"{'Pair':>10}  {'RMSD':>12}  {'Max |diff|':>12}")
print("-" * 40)
for i in range(1, 5):
    for j in range(i+1, 5):
        diff = spline_inner[i] - spline_inner[j]
        rmsd = np.sqrt(np.mean(diff**2))
        maxd = np.max(np.abs(diff))
        print(f"BC{i} vs BC{j}  {rmsd:12.6e}  {maxd:12.6e}")

# 全区间偏差
print("\n" + "=" * 80)
print("Pairwise comparison on full interval [0.0, 1.0]")
print("=" * 80)
print(f"{'Pair':>10}  {'RMSD':>12}  {'Max |diff|':>12}")
print("-" * 40)
for i in range(1, 5):
    for j in range(i+1, 5):
        diff = spline_full[i] - spline_full[j]
        rmsd = np.sqrt(np.mean(diff**2))
        maxd = np.max(np.abs(diff))
        print(f"BC{i} vs BC{j}  {rmsd:12.6e}  {maxd:12.6e}")

# ============================================================
# 5. 可视化 (SCI 学术风格, 多子图)
# ============================================================
colors = ['#1f77b4', '#d62728', '#2ca02c', '#ff7f0e']
linestyles = ['-', '--', '-.', ':']
markers_bc = ['s', '^', 'v', 'D']

fig = plt.figure(figsize=(7.2, 9.0))  # 单栏宽度 ~3.5in, 双栏 ~7.2in
gs = gridspec.GridSpec(3, 2, hspace=0.35, wspace=0.35)

# --- (a) 全区间样条曲线 ---
ax1 = fig.add_subplot(gs[0, :])
ax1.plot(x, y, 'ko', markersize=5, zorder=5, label='Data')
for idx, bc in enumerate(range(1, 5)):
    ax1.plot(x_full, spline_full[bc], color=colors[idx], ls=linestyles[idx],
             lw=1.2, label=f'BC {bc}')
ax1.axvline(0.05, color='gray', ls=':', lw=0.6, alpha=0.7)
ax1.annotate('$x=0.05$', xy=(0.05, 0.05), xycoords=('data', 'axes fraction'),
             fontsize=8, color='gray')
ax1.set_xlabel('$t$')
ax1.set_ylabel('Brightness')
ax1.set_title('(a) Cubic spline interpolation — full interval')
ax1.legend(loc='upper right', ncol=2)
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator())

# --- (b) [0.2, 0.8] 局部放大 ---
ax2 = fig.add_subplot(gs[1, 0])
for idx, bc in enumerate(range(1, 5)):
    ax2.plot(x_inner, spline_inner[bc], color=colors[idx], ls=linestyles[idx], lw=1.2,
             label=f'BC {bc}')
# 标注节点
mask = (x >= 0.2) & (x <= 0.8)
ax2.plot(x[mask], y[mask], 'ko', markersize=4, zorder=5)
ax2.set_xlabel('$t$')
ax2.set_ylabel('Brightness')
ax2.set_title('(b) Zoom-in: $t \\in [0.2, 0.8]$')
ax2.legend(loc='best', fontsize=8)
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())

# --- (c) 端点附近放大 [0, 0.2] ---
x_left = np.linspace(0.0, 0.25, 200)
ax3 = fig.add_subplot(gs[1, 1])
for idx, bc in enumerate(range(1, 5)):
    r = results[bc]
    ax3.plot(x_left, eval_spline(x_left, r['a'], r['b'], r['c'], r['d']),
             color=colors[idx], ls=linestyles[idx], lw=1.2, label=f'BC {bc}')
ax3.plot(x[:2], y[:2], 'ko', markersize=4, zorder=5)
ax3.axvline(0.05, color='gray', ls=':', lw=0.6, alpha=0.7)
ax3.set_xlabel('$t$')
ax3.set_ylabel('Brightness')
ax3.set_title('(c) Zoom-in: $t \\in [0, 0.25]$')
ax3.legend(loc='best', fontsize=8)
ax3.xaxis.set_minor_locator(AutoMinorLocator())
ax3.yaxis.set_minor_locator(AutoMinorLocator())

# --- (d) 差异曲线 (相对于均值基准) 全区间 ---
ax4 = fig.add_subplot(gs[2, 0])
for idx, bc in enumerate(range(1, 5)):
    diff = spline_full[bc] - ref_full
    ax4.plot(x_full, diff, color=colors[idx], ls=linestyles[idx], lw=1.0,
             label=f'BC {bc}')
ax4.axhline(0, color='k', lw=0.4)
ax4.set_xlabel('$t$')
ax4.set_ylabel('$S_{\\mathrm{BC}_i}(t) - \\bar{S}(t)$')
ax4.set_title('(d) Deviation from mean — full interval')
ax4.legend(loc='best', fontsize=8)
ax4.xaxis.set_minor_locator(AutoMinorLocator())
ax4.yaxis.set_minor_locator(AutoMinorLocator())

# --- (e) RMSD & Max deviation 柱状图 ---
ax5 = fig.add_subplot(gs[2, 1])
pairs = []
rmsd_vals = []
maxd_vals = []
for i in range(1, 5):
    for j in range(i+1, 5):
        diff = spline_full[i] - spline_full[j]
        pairs.append(f'{i}-{j}')
        rmsd_vals.append(np.sqrt(np.mean(diff**2)))
        maxd_vals.append(np.max(np.abs(diff)))

x_bar = np.arange(len(pairs))
width = 0.35
bars1 = ax5.bar(x_bar - width/2, rmsd_vals, width, color='#4c72b0', edgecolor='black',
                linewidth=0.5, label='RMSD')
bars2 = ax5.bar(x_bar + width/2, maxd_vals, width, color='#dd8452', edgecolor='black',
                linewidth=0.5, label='Max $|\\Delta|$')
ax5.set_xticks(x_bar)
ax5.set_xticklabels([f'BC{p}' for p in pairs], fontsize=8)
ax5.set_ylabel('Deviation')
ax5.set_title('(e) Pairwise deviation (full interval)')
ax5.legend(loc='upper left', fontsize=8)
ax5.yaxis.set_minor_locator(AutoMinorLocator())

plt.savefig('cubic_spline_comparison.png')
plt.show()
print("\n[Figure saved: cubic_spline_comparison.png]")
