"""
第八章 乘幂法求矩阵主特征值及对应特征向量
Power Method for Dominant Eigenvalue and Eigenvector

矩阵 A = [[3, -4, 3], [-4, 6, 3], [3, 3, 1]]
收敛条件: |λ₁^(k+1) - λ₁^(k)| ≤ 1e-3
初始向量: V^(0) = (1, 1, 1)^T
每5步规范化一次
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# ============================================================
# 1. 乘幂法核心算法
# ============================================================

def power_method(A, v0, tol=1e-3, max_iter=200, norm_interval=5):
    """
    乘幂法（带规范化）求按模最大特征值及对应特征向量。

    Parameters
    ----------
    A : ndarray, shape (n, n)
        输入矩阵
    v0 : ndarray, shape (n,)
        初始向量
    tol : float
        收敛容限 |λ^(k+1) - λ^(k)| ≤ tol
    max_iter : int
        最大迭代次数
    norm_interval : int
        每隔多少步做一次规范化（按模最大分量归一）

    Returns
    -------
    eigenvalue : float
        主特征值近似
    eigenvector : ndarray
        对应特征向量近似
    history : list of dict
        每步迭代记录 {k, v, lam}
    """
    n = len(v0)
    v = v0.astype(float).copy()
    history = [{"k": 0, "v": v.copy(), "v_raw": None, "lam": None}]

    lam_prev = None
    for k in range(1, max_iter + 1):
        v_new = A @ v

        # 计算特征值近似: 取第一个分量的比值 (与教材一致)
        lam = v_new[0] / v[0]

        # 每 norm_interval 步规范化，同时保存规范化前的向量
        v_raw = None
        if k % norm_interval == 0:
            v_raw = v_new.copy()
            max_component = v_new[np.argmax(np.abs(v_new))]
            v_new = v_new / max_component

        history.append({"k": k, "v": v_new.copy(), "v_raw": v_raw, "lam": lam})
        v = v_new

        # 收敛判断
        if lam_prev is not None and abs(lam - lam_prev) <= tol:
            break
        lam_prev = lam

    return lam, v, history


# ============================================================
# 2. 终端表格输出
# ============================================================

def print_iteration_table(history):
    """打印迭代过程表格，规范化步骤同时输出规范化前后的向量"""
    header = f"{'k':>4s} | {'':>6s} | {'V1':>14s}  {'V2':>14s}  {'V3':>14s} | {'lambda_1':>14s}"
    sep = "-" * len(header)
    print("\n" + sep)
    print("  Power Method Iteration")
    print(sep)
    print(header)
    print(sep)
    for rec in history:
        k = rec["k"]
        v = rec["v"]
        v_raw = rec["v_raw"]
        lam = rec["lam"]
        lam_str = f"{lam:>14.7f}" if lam is not None else f"{'---':>14s}"

        if v_raw is not None:
            # 规范化前
            raw_str = f"{v_raw[0]:>14.4f}  {v_raw[1]:>14.4f}  {v_raw[2]:>14.4f}"
            print(f"{k:>4d} | {'raw':>6s} | {raw_str} | {lam_str}")
            # 规范化后
            norm_str = f"{v[0]:>14.7f}  {v[1]:>14.7f}  {v[2]:>14.7f}"
            print(f"{'':>4s} | {'norm':>6s} | {norm_str} | {'':>14s}")
        else:
            v_str = f"{v[0]:>14.7f}  {v[1]:>14.7f}  {v[2]:>14.7f}"
            print(f"{k:>4d} | {'':>6s} | {v_str} | {lam_str}")
    print(sep)


# ============================================================
# 3. SCI 学术风格可视化
# ============================================================

def plot_convergence(history, true_eigenvalue, save_path="power_method_convergence.png"):
    """绘制特征值收敛曲线和误差曲线（SCI学术风格）"""

    # SCI 风格全局设置
    rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman"],
        "mathtext.fontset": "stix",
        "font.size": 11,
        "axes.linewidth": 1.0,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "legend.fontsize": 10,
        "legend.frameon": True,
        "legend.edgecolor": "black",
        "figure.dpi": 150,
    })

    # 提取有效数据（跳过 k=0）
    ks = [rec["k"] for rec in history if rec["lam"] is not None]
    lams = [rec["lam"] for rec in history if rec["lam"] is not None]
    errors = [abs(l - true_eigenvalue) for l in lams]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    # --- 左图: 特征值逐步逼近 ---
    ax1.plot(ks, lams, "o-", color="#2166AC", markersize=4, linewidth=1.2,
             label=r"$\lambda_1^{(k)}$")
    ax1.axhline(y=true_eigenvalue, color="#B2182B", linestyle="--", linewidth=1.0,
                label=rf"Exact $\lambda_1 = {true_eigenvalue:.4f}$")
    ax1.set_xlabel("Iteration $k$")
    ax1.set_ylabel(r"Eigenvalue approximation $\lambda_1^{(k)}$")
    ax1.set_title("(a) Convergence of dominant eigenvalue")
    ax1.legend(loc="best")
    ax1.grid(True, linestyle=":", linewidth=0.5, alpha=0.7)

    # --- 右图: 误差半对数 ---
    ax2.semilogy(ks, errors, "s-", color="#D6604D", markersize=4, linewidth=1.2,
                 label=r"$|\lambda_1^{(k)} - \lambda_1|$")
    ax2.set_xlabel("Iteration $k$")
    ax2.set_ylabel("Absolute error")
    ax2.set_title("(b) Error decay (semi-log scale)")
    ax2.legend(loc="best")
    ax2.grid(True, linestyle=":", linewidth=0.5, alpha=0.7)

    fig.tight_layout(pad=2.0)
    fig.savefig(save_path, bbox_inches="tight", dpi=200)
    plt.show()
    print(f"\nFigure saved: {save_path}")


# ============================================================
# 4. 主程序
# ============================================================

if __name__ == "__main__":
    # 定义矩阵和初始向量
    A = np.array([
        [ 3, -4,  3],
        [-4,  6,  3],
        [ 3,  3,  1]
    ], dtype=float)
    v0 = np.array([1.0, 1.0, 1.0])

    # 用 numpy 求精确特征值作为参考
    exact_eigenvalues = np.linalg.eigvals(A)
    true_lam1 = max(exact_eigenvalues, key=abs)
    print(f"numpy 精确特征值: {sorted(exact_eigenvalues, key=abs, reverse=True)}")
    print(f"按模最大特征值 (参考): {true_lam1:.10f}")

    # 运行乘幂法
    lam, vec, history = power_method(A, v0, tol=1e-3, max_iter=200, norm_interval=5)

    # 打印迭代表格
    print_iteration_table(history)

    # 打印最终结果
    print(f"\n主特征值 lam_1 = {lam:.7f}")
    print(f"对应特征向量 V = ({vec[0]:.7f}, {vec[1]:.7f}, {vec[2]:.7f})^T")
    print(f"迭代次数: {history[-1]['k']}")

    # 可视化
    plot_convergence(history, true_lam1)
