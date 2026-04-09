import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'KaiTi']
matplotlib.rcParams['axes.unicode_minus'] = False

# ==================== 数据 ====================
years = np.array([1960, 1961, 1962, 1963, 1964, 1965, 1966, 1967, 1968], dtype=float)
population = np.array([29.72, 30.61, 31.51, 32.13, 32.34, 32.85, 33.56, 34.20, 34.83], dtype=float)

t = years - 1960  # t = year - 1960
y = population
n = len(t)

actual_2000 = 60.55
t_predict = 2000 - 1960


# ==================== 最小二乘法核心函数 ====================
def least_squares_fit(A, b):
    """通过法方程 (A^T A) x = A^T b 求解最小二乘问题"""
    ATA = A.T @ A
    ATb = A.T @ b
    x = np.linalg.solve(ATA, ATb)
    return x


def compute_errors(y_true, y_fit):
    """计算均方偏差(RMSE)和最大偏差"""
    residuals = y_true - y_fit
    rmse = np.sqrt(np.sum(residuals**2) / len(y_true))
    max_dev = np.max(np.abs(residuals))
    return rmse, max_dev


def gauss_newton_logistic(t, y, K0=100.0, b0=2.0, r0=0.05, max_iter=200, tol=1e-10):
    """
    Gauss-Newton 迭代求解 Logistic 模型: y = K / (1 + b * exp(-r * t))
    返回: (K, b, r)
    """
    K, b, r = K0, b0, r0
    for _ in range(max_iter):
        exp_rt = np.exp(-r * t)
        denom = 1 + b * exp_rt
        f = K / denom
        res = y - f
        # Jacobian
        J_K = 1.0 / denom
        J_b = -K * exp_rt / denom**2
        J_r = K * b * t * exp_rt / denom**2
        J = np.column_stack([J_K, J_b, J_r])
        delta = np.linalg.lstsq(J, res, rcond=None)[0]
        K += delta[0]
        b += delta[1]
        r += delta[2]
        if np.linalg.norm(delta) < tol:
            break
    return K, b, r


# ==================== 模型拟合 ====================
results = {}  # {name: (formula_str, coeffs_str, y_fit, y_pred, rmse, max_dev)}

# --- 模型1: 线性模型 y = a0 + a1*t ---
A1 = np.column_stack([np.ones(n), t])
c1 = least_squares_fit(A1, y)
y_fit1 = A1 @ c1
y_pred1 = c1[0] + c1[1] * t_predict
rmse1, maxd1 = compute_errors(y, y_fit1)
results['Linear'] = (
    'y = a0 + a1*t',
    f'a0={c1[0]:.6f}, a1={c1[1]:.6f}',
    y_fit1, y_pred1, rmse1, maxd1
)

# --- 模型2: 二次多项式 y = a0 + a1*t + a2*t^2 ---
A2 = np.column_stack([np.ones(n), t, t**2])
c2 = least_squares_fit(A2, y)
y_fit2 = A2 @ c2
y_pred2 = c2[0] + c2[1] * t_predict + c2[2] * t_predict**2
rmse2, maxd2 = compute_errors(y, y_fit2)
results['Quadratic'] = (
    'y = a0 + a1*t + a2*t^2',
    f'a0={c2[0]:.6f}, a1={c2[1]:.6f}, a2={c2[2]:.6f}',
    y_fit2, y_pred2, rmse2, maxd2
)

# --- 模型3: 指数模型 y = e^(a + b*t) ---
ln_y = np.log(y)
A3 = np.column_stack([np.ones(n), t])
c3 = least_squares_fit(A3, ln_y)
y_fit3 = np.exp(c3[0] + c3[1] * t)
y_pred3 = np.exp(c3[0] + c3[1] * t_predict)
rmse3, maxd3 = compute_errors(y, y_fit3)
results['Exponential'] = (
    'y = exp(a + b*t)',
    f'a={c3[0]:.6f}, b={c3[1]:.6f}',
    y_fit3, y_pred3, rmse3, maxd3
)

# --- 模型4: 幂函数模型 y = a*(t+1)^b ---
ln_t1 = np.log(t + 1)
A4 = np.column_stack([np.ones(n), ln_t1])
c4 = least_squares_fit(A4, ln_y)
a4_val = np.exp(c4[0])
y_fit4 = a4_val * (t + 1)**c4[1]
y_pred4 = a4_val * (t_predict + 1)**c4[1]
rmse4, maxd4 = compute_errors(y, y_fit4)
results['Power'] = (
    'y = a*(t+1)^b',
    f'a={a4_val:.6f}, b={c4[1]:.6f}',
    y_fit4, y_pred4, rmse4, maxd4
)

# --- 模型5: Logistic模型 y = K / (1 + b*exp(-r*t)) (Gauss-Newton) ---
# 尝试多组初始值, 选择残差最小的结果
best_logistic = None
best_logistic_rmse = np.inf
for K0 in [80, 120, 200, 500]:
    for r0 in [0.02, 0.05, 0.1]:
        try:
            b0 = K0 / y[0] - 1  # 由 y(0) = K/(1+b) 估计 b0
            Kt, bt, rt = gauss_newton_logistic(t, y, K0=K0, b0=b0, r0=r0)
            if Kt > 0 and bt > 0 and rt > 0:
                yf = Kt / (1 + bt * np.exp(-rt * t))
                rmse_t, _ = compute_errors(y, yf)
                if rmse_t < best_logistic_rmse:
                    best_logistic_rmse = rmse_t
                    best_logistic = (Kt, bt, rt)
        except Exception:
            continue

K5, b5, r5 = best_logistic
y_fit5 = K5 / (1 + b5 * np.exp(-r5 * t))
y_pred5 = K5 / (1 + b5 * np.exp(-r5 * t_predict))
rmse5, maxd5 = compute_errors(y, y_fit5)
results['Logistic'] = (
    'y = K / (1 + b*exp(-r*t))',
    f'K={K5:.6f}, b={b5:.6f}, r={r5:.6f}',
    y_fit5, y_pred5, rmse5, maxd5
)


# ==================== 数据输出 ====================
print("=" * 90)
print(f"{'Model':<14} {'Formula':<28} {'RMSE':<10} {'MaxDev':<10} {'Pred(2000)':<12} {'RelErr(%)':<10}")
print("=" * 90)
for name, (formula, coeffs, yf, yp, rmse, maxd) in results.items():
    rel_err = abs(yp - actual_2000) / actual_2000 * 100
    print(f"{name:<14} {formula:<28} {rmse:<10.6f} {maxd:<10.6f} {yp:<12.4f} {rel_err:<10.4f}")
print("=" * 90)
print(f"Actual population in 2000: {actual_2000} (unit: 100 million)")
print()

# 输出各模型拟合系数
print("-" * 60)
print("Fitted coefficients:")
print("-" * 60)
for name, (formula, coeffs, *_) in results.items():
    print(f"  {name:<14} {coeffs}")
print()

# 输出各数据点的拟合值与残差
print("-" * 90)
header = f"{'Year':<6}"
for name in results:
    header += f" {name+'_fit':<12} {name+'_res':<12}"
print(header)
print("-" * 90)
for i in range(n):
    row = f"{int(years[i]):<6}"
    for name, (_, _, yf, *__) in results.items():
        res_i = y[i] - yf[i]
        row += f" {yf[i]:<12.4f} {res_i:<12.4f}"
    print(row)
print("-" * 90)
print()

# 输出法方程矩阵 (以线性模型为例)
print("Normal equation (A^T A)x = A^T b for Linear model:")
print(f"  A^T A = \n{A1.T @ A1}")
print(f"  A^T b = {A1.T @ y}")
print(f"  x     = {c1}")
print()


# ==================== 学术风格可视化 ====================
plt.style.use('seaborn-v0_8-whitegrid')
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

t_plot = np.linspace(0, 42, 300)

# 各模型预测曲线
y_plots = {
    'Linear':      c1[0] + c1[1] * t_plot,
    'Quadratic':   c2[0] + c2[1] * t_plot + c2[2] * t_plot**2,
    'Exponential': np.exp(c3[0] + c3[1] * t_plot),
    'Power':       a4_val * (t_plot + 1)**c4[1],
    'Logistic':    K5 / (1 + b5 * np.exp(-r5 * t_plot)),
}

colors = {
    'Linear': '#1f77b4', 'Quadratic': '#ff7f0e',
    'Exponential': '#2ca02c', 'Power': '#9467bd', 'Logistic': '#d62728'
}
linestyles = {
    'Linear': '-', 'Quadratic': '--',
    'Exponential': '-.', 'Power': ':', 'Logistic': '-'
}

# --- 图(a): 拟合曲线与外推预测 ---
ax = axes[0]
ax.scatter(t, y, c='k', s=40, zorder=5, marker='o', label='Data (1960-1968)', edgecolors='k', linewidths=0.5)
for name, yp in y_plots.items():
    ax.plot(t_plot, yp, color=colors[name], ls=linestyles[name], lw=1.5, label=name)
ax.scatter([t_predict], [actual_2000], c='red', s=120, zorder=5, marker='*', label=f'Actual 2000: {actual_2000}')
# 标注各模型预测点
for name, (_, _, _, yp, *__) in results.items():
    ax.scatter([t_predict], [yp], c=colors[name], s=50, zorder=4, marker='s', edgecolors='k', linewidths=0.3)
ax.set_xlabel(r'$t$ (year $-$ 1960)', fontsize=11)
ax.set_ylabel('Population ($10^8$)', fontsize=11)
ax.set_title('(a) Least Squares Fitting & Prediction', fontsize=12)
ax.legend(fontsize=7.5, loc='upper left', framealpha=0.9, edgecolor='gray')
ax.set_xlim(-1, 44)
ax.tick_params(labelsize=9)

# --- 图(b): 拟合残差 ---
ax = axes[1]
for name, (_, _, yf, *__) in results.items():
    res = y - yf
    ax.plot(t, res, color=colors[name], marker='o', ms=4, lw=1.2,
            ls=linestyles[name], label=name)
ax.axhline(y=0, color='k', lw=0.8, ls='-')
ax.set_xlabel(r'$t$ (year $-$ 1960)', fontsize=11)
ax.set_ylabel(r'Residual ($y_i - \hat{y}_i$)', fontsize=11)
ax.set_title('(b) Fitting Residuals', fontsize=12)
ax.legend(fontsize=7.5, framealpha=0.9, edgecolor='gray')
ax.set_xticks(t)
ax.set_xticklabels([str(int(yr)) for yr in years], fontsize=8, rotation=45)
ax.tick_params(labelsize=9)

# --- 图(c): 模型误差对比柱状图 ---
ax = axes[2]
model_names = list(results.keys())
rmse_vals = [results[m][4] for m in model_names]
maxd_vals = [results[m][5] for m in model_names]
rel_errs = [abs(results[m][3] - actual_2000) / actual_2000 * 100 for m in model_names]

x_idx = np.arange(len(model_names))
w = 0.25
bars1 = ax.bar(x_idx - w, rmse_vals, w, label='RMSE', color='#4c72b0', edgecolor='k', linewidth=0.5)
bars2 = ax.bar(x_idx, maxd_vals, w, label='Max Dev.', color='#55a868', edgecolor='k', linewidth=0.5)

ax2 = ax.twinx()
bars3 = ax2.bar(x_idx + w, rel_errs, w, label='Rel. Err. (%)', color='#c44e52', edgecolor='k', linewidth=0.5)

ax.set_xlabel('Model', fontsize=11)
ax.set_ylabel('Error ($10^8$)', fontsize=11)
ax2.set_ylabel('Relative Error (%)', fontsize=11)
ax.set_title('(c) Model Comparison', fontsize=12)
ax.set_xticks(x_idx)
ax.set_xticklabels(model_names, fontsize=8.5, rotation=15)
ax.tick_params(labelsize=9)
ax2.tick_params(labelsize=9)

# 合并图例
lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines1 + lines2, labels1 + labels2, fontsize=7.5, loc='upper left',
          framealpha=0.9, edgecolor='gray')

plt.tight_layout(pad=2.0)
plt.savefig('population_prediction.png', dpi=300, bbox_inches='tight')
plt.show()

print("Figure saved: population_prediction.png")
