"""
计算方法H第7章实践作业
求解常微分方程初值问题: y' = -xy^2, y(0) = 2, x ∈ [0, 5]
使用三种方法:
1. Euler预估-校正方法（改进Euler法）
2. 经典四阶Runge-Kutta方法
3. 四阶Adams预估-校正方法
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# 设置学术风格的绘图参数
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 10
rcParams['axes.labelsize'] = 11
rcParams['axes.titlesize'] = 12
rcParams['xtick.labelsize'] = 9
rcParams['ytick.labelsize'] = 9
rcParams['legend.fontsize'] = 9
rcParams['figure.dpi'] = 300
rcParams['savefig.dpi'] = 300
rcParams['savefig.bbox'] = 'tight'


def f(x, y):
    """微分方程右端函数: y' = -xy^2"""
    return -x * y**2


def exact_solution(x):
    """
    理论解析解
    y' = -xy^2 可分离变量求解
    dy/y^2 = -x dx
    -1/y = -x^2/2 + C
    y(0) = 2 => C = -1/2
    y = 2/(1 + x^2)
    """
    return 2 / (1 + x**2)


def euler_predictor_corrector(f, x0, y0, x_end, h):
    """
    Euler预估-校正方法（改进Euler法）
    预估: y_{n+1}^(0) = y_n + h*f(x_n, y_n)
    校正: y_{n+1} = y_n + h/2 * [f(x_n, y_n) + f(x_{n+1}, y_{n+1}^(0))]
    """
    n = int((x_end - x0) / h)
    x = np.zeros(n + 1)
    y = np.zeros(n + 1)
    x[0], y[0] = x0, y0

    for i in range(n):
        x[i+1] = x[i] + h
        # 预估
        y_pred = y[i] + h * f(x[i], y[i])
        # 校正
        y[i+1] = y[i] + h/2 * (f(x[i], y[i]) + f(x[i+1], y_pred))

    return x, y


def runge_kutta_4(f, x0, y0, x_end, h):
    """
    经典四阶Runge-Kutta方法
    y_{n+1} = y_n + h/6 * (K1 + 2*K2 + 2*K3 + K4)
    K1 = f(x_n, y_n)
    K2 = f(x_n + h/2, y_n + h/2*K1)
    K3 = f(x_n + h/2, y_n + h/2*K2)
    K4 = f(x_n + h, y_n + h*K3)
    """
    n = int((x_end - x0) / h)
    x = np.zeros(n + 1)
    y = np.zeros(n + 1)
    x[0], y[0] = x0, y0

    for i in range(n):
        x[i+1] = x[i] + h
        K1 = f(x[i], y[i])
        K2 = f(x[i] + h/2, y[i] + h/2*K1)
        K3 = f(x[i] + h/2, y[i] + h/2*K2)
        K4 = f(x[i] + h, y[i] + h*K3)
        y[i+1] = y[i] + h/6 * (K1 + 2*K2 + 2*K3 + K4)

    return x, y


def adams_predictor_corrector(f, x0, y0, x_end, h):
    """
    四阶Adams预估-校正方法
    需要前4个点，使用RK4方法计算前3个点
    预估（Adams外插公式）:
    y_{n+1}^(0) = y_n + h/24 * [55*f_n - 59*f_{n-1} + 37*f_{n-2} - 9*f_{n-3}]
    校正（Adams内插公式）:
    y_{n+1} = y_n + h/24 * [9*f_{n+1}^(0) + 19*f_n - 5*f_{n-1} + f_{n-2}]
    """
    n = int((x_end - x0) / h)
    x = np.zeros(n + 1)
    y = np.zeros(n + 1)
    x[0], y[0] = x0, y0

    # 使用RK4计算前3个点
    for i in range(min(3, n)):
        x[i+1] = x[i] + h
        K1 = f(x[i], y[i])
        K2 = f(x[i] + h/2, y[i] + h/2*K1)
        K3 = f(x[i] + h/2, y[i] + h/2*K2)
        K4 = f(x[i] + h, y[i] + h*K3)
        y[i+1] = y[i] + h/6 * (K1 + 2*K2 + 2*K3 + K4)

    # Adams预估-校正方法
    for i in range(3, n):
        x[i+1] = x[i] + h
        # 预估
        f_vals = [f(x[i-j], y[i-j]) for j in range(4)]
        y_pred = y[i] + h/24 * (55*f_vals[0] - 59*f_vals[1] + 37*f_vals[2] - 9*f_vals[3])
        # 校正
        f_pred = f(x[i+1], y_pred)
        y[i+1] = y[i] + h/24 * (9*f_pred + 19*f_vals[0] - 5*f_vals[1] + f_vals[2])

    return x, y


def compute_errors(x, y_numerical, y_exact):
    """计算数值解的误差"""
    errors = np.abs(y_numerical - y_exact)
    max_error = np.max(errors)
    mean_error = np.mean(errors)
    return errors, max_error, mean_error


def main():
    """主程序"""
    x0, y0 = 0, 2
    x_end = 5
    step_sizes = [0.5, 0.2, 0.1, 0.05]
    results = {}

    print("=" * 70)
    print("常微分方程初值问题数值解法比较")
    print("问题: y' = -xy^2, y(0) = 2, x ∈ [0, 5]")
    print("解析解: y(x) = 2/(1 + x^2)")
    print("=" * 70)

    for h in step_sizes:
        print(f"\n步长 h = {h}")
        print("-" * 70)

        # 三种方法求解
        x_euler, y_euler = euler_predictor_corrector(f, x0, y0, x_end, h)
        x_rk4, y_rk4 = runge_kutta_4(f, x0, y0, x_end, h)
        x_adams, y_adams = adams_predictor_corrector(f, x0, y0, x_end, h)

        # 计算精确解
        y_exact_euler = exact_solution(x_euler)
        y_exact_rk4 = exact_solution(x_rk4)
        y_exact_adams = exact_solution(x_adams)

        # 输出数值结果表格
        print(f"\n{'n':<5} {'x_n':<10} {'y(x_n)':<12} {'Euler P-C':<12} {'RK4':<12} {'Adams P-C':<12}")
        print("-" * 70)
        for i in range(len(x_euler)):
            print(f"{i:<5} {x_euler[i]:<10.4f} {y_exact_euler[i]:<12.8f} {y_euler[i]:<12.8f} {y_rk4[i]:<12.8f} {y_adams[i]:<12.8f}")

        # 计算误差
        _, max_err_euler, mean_err_euler = compute_errors(x_euler, y_euler, y_exact_euler)
        _, max_err_rk4, mean_err_rk4 = compute_errors(x_rk4, y_rk4, y_exact_rk4)
        _, max_err_adams, mean_err_adams = compute_errors(x_adams, y_adams, y_exact_adams)

        # 输出误差统计
        print(f"\n{'方法':<25} {'最大误差':<15} {'平均误差':<15}")
        print(f"Euler预估-校正法        {max_err_euler:.6e}    {mean_err_euler:.6e}")
        print(f"四阶Runge-Kutta法       {max_err_rk4:.6e}    {mean_err_rk4:.6e}")
        print(f"四阶Adams预估-校正法    {max_err_adams:.6e}    {mean_err_adams:.6e}")

        # 存储结果用于绘图
        results[h] = {
            'euler': (x_euler, y_euler, y_exact_euler),
            'rk4': (x_rk4, y_rk4, y_exact_rk4),
            'adams': (x_adams, y_adams, y_exact_adams)
        }

    # 绘制比较图
    plot_comparison(results, step_sizes)
    plot_error_analysis(results, step_sizes)

    print("\n" + "=" * 70)
    print("计算完成！图像已保存为 PNG 文件。")
    print("=" * 70)


def plot_comparison(results, step_sizes):
    """绘制不同步长下三种方法的比较图"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for idx, h in enumerate(step_sizes):
        ax = axes[idx]
        x_euler, y_euler, y_exact_euler = results[h]['euler']
        x_rk4, y_rk4, y_exact_rk4 = results[h]['rk4']
        x_adams, y_adams, y_exact_adams = results[h]['adams']

        # 绘制精确解
        x_fine = np.linspace(0, 5, 1000)
        y_fine = exact_solution(x_fine)
        ax.plot(x_fine, y_fine, 'k-', linewidth=1.5, label='Exact solution', alpha=0.7)

        # 绘制数值解
        ax.plot(x_euler, y_euler, 'o-', markersize=4, linewidth=1, label='Euler P-C', alpha=0.8)
        ax.plot(x_rk4, y_rk4, 's-', markersize=4, linewidth=1, label='RK4', alpha=0.8)
        ax.plot(x_adams, y_adams, '^-', markersize=4, linewidth=1, label='Adams P-C', alpha=0.8)

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(f'h = {h}')
        ax.legend(loc='best', frameon=True, fancybox=False, edgecolor='black')
        ax.grid(True, linestyle='--', alpha=0.3)
        ax.set_xlim([0, 5])

    plt.tight_layout()
    plt.savefig('E:\\_Materials_in_NWPU_\\2026_Spring\\course\\Computational Methods H\\week7\\cm_week7\\comparison.png')
    print("\n数值解比较图已保存: comparison.png")
    plt.close()


def plot_error_analysis(results, step_sizes):
    """绘制误差分析图"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # 左图：不同步长下的最大误差
    max_errors = {'Euler P-C': [], 'RK4': [], 'Adams P-C': []}
    for h in step_sizes:
        x_euler, y_euler, y_exact_euler = results[h]['euler']
        x_rk4, y_rk4, y_exact_rk4 = results[h]['rk4']
        x_adams, y_adams, y_exact_adams = results[h]['adams']

        _, max_err_euler, _ = compute_errors(x_euler, y_euler, y_exact_euler)
        _, max_err_rk4, _ = compute_errors(x_rk4, y_rk4, y_exact_rk4)
        _, max_err_adams, _ = compute_errors(x_adams, y_adams, y_exact_adams)

        max_errors['Euler P-C'].append(max_err_euler)
        max_errors['RK4'].append(max_err_rk4)
        max_errors['Adams P-C'].append(max_err_adams)

    ax1 = axes[0]
    ax1.loglog(step_sizes, max_errors['Euler P-C'], 'o-', label='Euler P-C', linewidth=2, markersize=6)
    ax1.loglog(step_sizes, max_errors['RK4'], 's-', label='RK4', linewidth=2, markersize=6)
    ax1.loglog(step_sizes, max_errors['Adams P-C'], '^-', label='Adams P-C', linewidth=2, markersize=6)
    ax1.set_xlabel('Step size h')
    ax1.set_ylabel('Maximum error')
    ax1.set_title('Maximum Error vs Step Size')
    ax1.legend(loc='best', frameon=True, fancybox=False, edgecolor='black')
    ax1.grid(True, which='both', linestyle='--', alpha=0.3)

    # 右图：h=0.1时的误差分布
    h_selected = 0.1
    x_euler, y_euler, y_exact_euler = results[h_selected]['euler']
    x_rk4, y_rk4, y_exact_rk4 = results[h_selected]['rk4']
    x_adams, y_adams, y_exact_adams = results[h_selected]['adams']

    err_euler, _, _ = compute_errors(x_euler, y_euler, y_exact_euler)
    err_rk4, _, _ = compute_errors(x_rk4, y_rk4, y_exact_rk4)
    err_adams, _, _ = compute_errors(x_adams, y_adams, y_exact_adams)

    ax2 = axes[1]
    ax2.semilogy(x_euler, err_euler, 'o-', label='Euler P-C', linewidth=1.5, markersize=4)
    ax2.semilogy(x_rk4, err_rk4, 's-', label='RK4', linewidth=1.5, markersize=4)
    ax2.semilogy(x_adams, err_adams, '^-', label='Adams P-C', linewidth=1.5, markersize=4)
    ax2.set_xlabel('x')
    ax2.set_ylabel('Absolute error')
    ax2.set_title(f'Error Distribution (h = {h_selected})')
    ax2.legend(loc='best', frameon=True, fancybox=False, edgecolor='black')
    ax2.grid(True, which='both', linestyle='--', alpha=0.3)

    plt.tight_layout()
    plt.savefig('E:\\_Materials_in_NWPU_\\2026_Spring\\course\\Computational Methods H\\week7\\cm_week7\\error_analysis.png')
    print("误差分析图已保存: error_analysis.png")
    plt.close()


if __name__ == '__main__':
    main()
