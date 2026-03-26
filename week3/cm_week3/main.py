"""
计算方法H 第三章实践作业 - 主脚本
求解四元线性方程组，比较各方法的结果

用法:
    python main.py                  # 全部方法执行并比较（含可视化）
    python main.py gauss            # 仅执行Gauss顺序消去法
    python main.py gauss_full       # 仅执行全主元Gauss消去法
    python main.py doolittle        # 仅执行Doolittle分解法
    python main.py jacobi           # 仅执行Jacobi迭代法
    python main.py gauss_seidel     # 仅执行Gauss-Seidel迭代法
    python main.py sor              # 仅执行SOR逐次超松弛迭代法
"""
import sys
import time
import warnings
import numpy as np

warnings.filterwarnings('ignore', category=RuntimeWarning)

from gauss import gauss_elimination
from gauss_full_pivot import gauss_full_pivot
from doolittle import doolittle_decomposition
from jacobi import jacobi_iteration
from gauss_seidel import gauss_seidel_iteration
from sor import sor_iteration, find_optimal_omega, scan_omega

TOL = 1e-5
MAX_ITER = 10000

METHOD_NAMES = {
    'gauss': 'Gauss顺序消去法', 'gauss_full': '全主元Gauss消去法',
    'doolittle': 'Doolittle分解法', 'jacobi': 'Jacobi迭代法',
    'gauss_seidel': 'Gauss-Seidel迭代法', 'sor': 'SOR逐次超松弛迭代法',
}
METHOD_NAMES_EN = {
    'gauss': 'Gauss', 'gauss_full': 'Full-Pivot Gauss',
    'doolittle': 'Doolittle', 'jacobi': 'Jacobi',
    'gauss_seidel': 'Gauss-Seidel', 'sor': 'SOR',
}


def generate_system(n):
    """生成n阶同形方程组: 对角元=-n, 非对角元=1, b=全1, 精确解=全-1"""
    A = np.ones((n, n)) - (n + 1) * np.eye(n)
    b = np.ones(n)
    exact = -np.ones(n)
    return A, b, exact

# ========== 求解与单方法输出 ==========
A4, b4, EXACT4 = generate_system(4)


def solve(mk, A_mat, b_vec, tol=TOL, record_history=False, omega=None):
    """统一求解接口，返回 (x, ok, iters, elapsed, info, history)"""
    t0 = time.perf_counter()
    hist = []
    iters = 0
    if mk == 'gauss':
        x, ok = gauss_elimination(A_mat, b_vec)
    elif mk == 'gauss_full':
        x, ok = gauss_full_pivot(A_mat, b_vec)
    elif mk == 'doolittle':
        x, ok = doolittle_decomposition(A_mat, b_vec)
    elif mk == 'jacobi':
        x, ok, iters, hist = jacobi_iteration(
            A_mat, b_vec, tol=tol, max_iter=MAX_ITER,
            record_history=record_history)
    elif mk == 'gauss_seidel':
        x, ok, iters, hist = gauss_seidel_iteration(
            A_mat, b_vec, tol=tol, max_iter=MAX_ITER,
            record_history=record_history)
    elif mk == 'sor':
        if omega is None:
            omega, _ = find_optimal_omega(A_mat, b_vec, tol=tol,
                                          max_iter=MAX_ITER)
        x, ok, iters, hist = sor_iteration(
            A_mat, b_vec, omega=omega, tol=tol, max_iter=MAX_ITER,
            record_history=record_history)
    else:
        return None, False, 0, 0, "", []
    elapsed = time.perf_counter() - t0

    # 构造info
    if mk in ('jacobi', 'gauss_seidel'):
        info = f"迭代{iters}次" if ok else f"未收敛({iters}次)"
    elif mk == 'sor':
        info = f"ω={omega:.2f}, 迭代{iters}次"
    else:
        info = ""
    return x, ok, iters, elapsed, info, hist


def print_solution(mk):
    """单方法详细输出（含计时、迭代逐步结果）"""
    need_hist = mk in ('jacobi', 'gauss_seidel', 'sor')
    omega = None
    if mk == 'sor':
        omega, _ = find_optimal_omega(A4, b4, tol=TOL, max_iter=MAX_ITER)
    x, ok, iters, elapsed, info, hist = solve(
        mk, A4, b4, record_history=need_hist, omega=omega)

    name = METHOD_NAMES[mk]
    print(f"\n{'='*70}")
    print(f"  {name}")
    print(f"{'='*70}")

    if need_hist and hist:
        n = len(EXACT4)
        hdr = "  k |" + "".join(f"{'x'+str(i+1):>12s}" for i in range(n))
        hdr += " | 精度"
        print(hdr)
        print("  " + "-" * (len(hdr) - 2))
        for k, xk, acc in hist:
            row = f"  {k:>2d} |"
            row += "".join(f"{v:>12.8f}" for v in xk)
            row += f" | {acc:.6e}"
            print(row)

    if ok and x is not None:
        err = np.linalg.norm(x - EXACT4, ord=np.inf)
        print(f"\n  解向量: x = [{', '.join(f'{v:.10f}' for v in x)}]")
        print(f"  精确解: x*= [{', '.join(f'{v:.1f}' for v in EXACT4)}]")
        print(f"  误差(∞范数): {err:.6e}")
        print(f"  执行时间: {elapsed:.6f} s")
        if info:
            print(f"  {info}")
    else:
        print(f"  求解失败  {info}")

# ========== 全部执行 ==========
def run_all():
    """全部方法求解4阶方程组，打印比较表格"""
    methods = list(METHOD_NAMES.keys())
    results = {}
    for mk in methods:
        x, ok, iters, elapsed, info, _ = solve(mk, A4, b4)
        err = np.linalg.norm(x - EXACT4, ord=np.inf) if (ok and x is not None) else float('inf')
        results[mk] = dict(x=x, ok=ok, iters=iters, err=err,
                           elapsed=elapsed, info=info)

    print(f"\n{'='*95}")
    print("  各方法求解结果比较 (n=4)")
    print(f"{'='*95}")
    print(f"{'方法':<22s} | {'x1':>11s} {'x2':>11s} {'x3':>11s} {'x4':>11s}"
          f" | {'误差':>10s} | {'时间/s':>10s}")
    print("-" * 95)
    for mk in methods:
        r = results[mk]
        name = METHOD_NAMES[mk]
        if r['ok'] and r['x'] is not None:
            xs = ''.join(f"{v:>11.7f} " for v in r['x'])
            print(f"{name:<22s} | {xs}| {r['err']:>10.3e} | {r['elapsed']:>10.6f}")
        else:
            print(f"{name:<22s} | {'求解失败':>48s} |")
    print(f"{'精确解':<22s} | "
          f"{''.join(f'{v:>11.7f} ' for v in EXACT4)}|")

    print(f"\n--- 迭代法详细信息 ---")
    for mk in ('jacobi', 'gauss_seidel', 'sor'):
        r = results[mk]
        print(f"  {METHOD_NAMES[mk]}: {r['info']}, 耗时 {r['elapsed']:.6f}s")

    return results


# ========== 要求4: 高阶方程组扩展比较 ==========
N_SCALING = [4, 8, 16, 32, 64, 128]


def run_scaling():
    """对同形高阶方程组，比较各方法执行时间与迭代次数"""
    methods = list(METHOD_NAMES.keys())
    data = {mk: {'ns': [], 'times': [], 'iters': []} for mk in methods}

    print(f"\n{'='*95}")
    print("  要求4: 高阶同形方程组各方法性能比较")
    print(f"{'='*95}")
    hdr = f"{'n':>5s}"
    for mk in methods:
        hdr += f" | {METHOD_NAMES_EN[mk]+'/s':>14s}"
    print(hdr)
    print("-" * (5 + 17 * len(methods)))

    for n in N_SCALING:
        An, bn, _ = generate_system(n)
        row = f"{n:>5d}"
        for mk in methods:
            # SOR用理论最优omega加速
            omega = None
            if mk == 'sor':
                rho_J = (n - 1) / n
                omega = 2.0 / (1.0 + np.sqrt(1.0 - rho_J ** 2))
            x, ok, iters, elapsed, info, _ = solve(
                mk, An, bn, omega=omega)
            data[mk]['ns'].append(n)
            data[mk]['times'].append(elapsed)
            data[mk]['iters'].append(iters if ok else None)
            if ok:
                row += f" | {elapsed:>14.6f}"
            else:
                row += f" | {'fail':>14s}"
        print(row)

    # 迭代法迭代次数表
    iter_methods = ['jacobi', 'gauss_seidel', 'sor']
    print(f"\n--- 迭代次数 ---")
    hdr2 = f"{'n':>5s}"
    for mk in iter_methods:
        hdr2 += f" | {METHOD_NAMES_EN[mk]:>14s}"
    print(hdr2)
    print("-" * (5 + 17 * len(iter_methods)))
    for i, n in enumerate(N_SCALING):
        row = f"{n:>5d}"
        for mk in iter_methods:
            it = data[mk]['iters'][i]
            row += f" | {it:>14d}" if it else f" | {'fail':>14s}"
        print(row)

    return data

# ========== 可视化 ==========
def _setup_rc():
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Times New Roman', 'DejaVu Serif'],
        'mathtext.fontset': 'stix', 'font.size': 11,
        'axes.linewidth': 1.0, 'axes.labelsize': 13,
        'axes.titlesize': 13, 'axes.unicode_minus': False,
        'xtick.direction': 'in', 'ytick.direction': 'in',
        'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
        'xtick.minor.visible': True, 'ytick.minor.visible': True,
        'xtick.top': True, 'ytick.right': True,
        'legend.frameon': True, 'legend.framealpha': 0.9,
        'legend.edgecolor': '0.6', 'legend.fontsize': 9,
        'lines.linewidth': 1.4, 'lines.markersize': 5,
        'grid.alpha': 0.25, 'grid.linewidth': 0.5,
        'grid.linestyle': '--',
        'savefig.dpi': 300, 'figure.dpi': 150,
    })


def visualize(scaling_data):
    """
    三子图:
      (a) 迭代法收敛曲线 (精度 vs 迭代次数)
      (b) SOR ω-迭代次数
      (c) 高阶方程组执行时间 vs n
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    _setup_rc()

    iter_keys = ['jacobi', 'gauss_seidel', 'sor']
    iter_colors = {'jacobi': '#ff7f0e', 'gauss_seidel': '#2ca02c',
                   'sor': '#d62728'}
    iter_markers = {'jacobi': 'D', 'gauss_seidel': 's', 'sor': '^'}

    fig, axes = plt.subplots(1, 3, figsize=(17, 5))

    # ---- (a) 收敛曲线 ----
    ax = axes[0]
    best_w, _ = find_optimal_omega(A4, b4, tol=TOL, max_iter=MAX_ITER)
    for mk in iter_keys:
        omega = best_w if mk == 'sor' else None
        _, _, _, _, _, hist = solve(mk, A4, b4, record_history=True,
                                    omega=omega)
        if hist:
            ks = [h[0] for h in hist]
            accs = [h[2] for h in hist]
            lbl = METHOD_NAMES_EN[mk]
            if mk == 'sor':
                lbl += rf' ($\omega$={best_w:.2f})'
            ax.semilogy(ks, accs, marker=iter_markers[mk],
                        color=iter_colors[mk], label=lbl, markersize=4)
    ax.axhline(y=TOL, color='gray', linestyle=':', linewidth=1,
               label=rf'tol=$10^{{-5}}$')
    ax.set_xlabel('Iteration $k$')
    ax.set_ylabel(r'$\max|x_i^{(k+1)}-x_i^{(k)}|$')
    ax.set_title('(a)  Convergence history', loc='left', fontweight='bold')
    ax.legend(loc='upper right')
    ax.grid(True)

    # ---- (b) SOR ω 曲线 ----
    ax = axes[1]
    omegas, iters_list = scan_omega(A4, b4, tol=TOL, max_iter=MAX_ITER,
                                    step=0.02)
    pw = [w for w, it in zip(omegas, iters_list) if it is not None]
    pi = [it for it in iters_list if it is not None]
    ax.plot(pw, pi, '-', color='#1f77b4', linewidth=1.5, label='SOR')
    ax.plot(best_w, min(pi), 'r*', markersize=14, markeredgecolor='black',
            markeredgewidth=0.5, zorder=5,
            label=rf'$\omega_{{opt}}={best_w:.2f}$')
    _, gs_ok, gs_it, _ = gauss_seidel_iteration(A4, b4, tol=TOL,
                                                  max_iter=MAX_ITER)
    if gs_ok:
        ax.plot(1.0, gs_it, 'D', color='#2ca02c', markersize=7,
                markeredgecolor='black', markeredgewidth=0.5, zorder=5,
                label=f'G-S ($\\omega$=1, iter={gs_it})')
    ax.set_xlabel(r'Relaxation factor $\omega$')
    ax.set_ylabel('Iterations')
    ax.set_title(r'(b)  SOR: iterations vs. $\omega$',
                 loc='left', fontweight='bold')
    ax.legend(loc='upper right')
    ax.grid(True)

    # ---- (c) 高阶扩展: 执行时间 vs n ----
    ax = axes[2]
    all_methods = list(METHOD_NAMES.keys())
    styles = {
        'gauss': ('-', 'o', '#1f77b4'), 'gauss_full': ('-', 's', '#9467bd'),
        'doolittle': ('-', '^', '#8c564b'),
        'jacobi': ('--', 'D', '#ff7f0e'), 'gauss_seidel': ('--', 'v', '#2ca02c'),
        'sor': ('--', 'p', '#d62728'),
    }
    for mk in all_methods:
        ns = scaling_data[mk]['ns']
        ts = scaling_data[mk]['times']
        ls, m, c = styles[mk]
        ax.loglog(ns, ts, linestyle=ls, marker=m, color=c,
                  label=METHOD_NAMES_EN[mk], markeredgecolor='black',
                  markeredgewidth=0.3)
    ax.set_xlabel('Matrix order $n$')
    ax.set_ylabel('Execution time (s)')
    ax.set_title('(c)  Scaling comparison', loc='left', fontweight='bold')
    ax.legend(loc='upper left', fontsize=8)
    ax.grid(True)

    plt.tight_layout(w_pad=2.5)
    plt.savefig('comparison.png', bbox_inches='tight')
    plt.close(fig)
    print("\n可视化结果已保存至 comparison.png")


# ========== 入口 ==========
def main():
    if len(sys.argv) > 1:
        mk = sys.argv[1]
        if mk not in METHOD_NAMES:
            print(f"未知方法: {mk}")
            print(f"可选: {', '.join(METHOD_NAMES.keys())}")
            sys.exit(1)
        print_solution(mk)
    else:
        run_all()
        scaling_data = run_scaling()
        visualize(scaling_data)


if __name__ == '__main__':
    main()
