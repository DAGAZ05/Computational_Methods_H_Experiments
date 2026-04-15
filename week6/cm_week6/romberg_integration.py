# -*- coding: utf-8 -*-
"""
Romberg积分算法实现
计算方法课程第六周作业
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import sys
import io

# 设置UTF-8编码输出
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

# 设置中文字体和学术风格
rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
rcParams['axes.unicode_minus'] = False  # 解决负号显示问题
rcParams['font.family'] = 'sans-serif'
plt.style.use('seaborn-v0_8-paper')

# 抑制matplotlib字体警告
import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')


def trapezoidal(f, a, b, n):
    """复化梯形公式"""
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    result = (f(a) + f(b)) / 2
    for i in range(1, n):
        result += f(x[i])
    return h * result


def romberg_integration(f, a, b, tol=1e-8, max_iter=20):
    """
    Romberg积分算法

    参数:
        f: 被积函数
        a, b: 积分区间
        tol: 误差容限
        max_iter: 最大迭代次数

    返回:
        result: 积分结果
        T_table: T数表
        sequences: 各序列值 (梯形、Simpson、柯特斯、Romberg)
    """
    # 初始化T数表
    T = np.zeros((max_iter, max_iter))

    # 第一列：梯形序列 T_n
    n = 1
    T[0, 0] = (b - a) * (f(a) + f(b)) / 2

    sequences = {
        'Trapezoidal': [T[0, 0]],
        'Simpson': [],
        'Cotes': [],
        'Romberg': []
    }

    for k in range(1, max_iter):
        # 区间逐次分半
        h = (b - a) / (2**k)
        sum_term = 0
        for i in range(1, 2**k, 2):
            sum_term += f(a + i * h)
        T[k, 0] = T[k-1, 0] / 2 + h * sum_term
        sequences['Trapezoidal'].append(T[k, 0])

        # 外推加速：逐次计算Simpson、柯特斯、Romberg序列
        for m in range(1, k + 1):
            T[k, m] = (4**m * T[k, m-1] - T[k-1, m-1]) / (4**m - 1)

            # 记录各序列
            if m == 1 and k >= 1:
                sequences['Simpson'].append(T[k, 1])
            elif m == 2 and k >= 2:
                sequences['Cotes'].append(T[k, 2])
            elif m == 3 and k >= 3:
                sequences['Romberg'].append(T[k, 3])

        # 检查收敛性（使用Romberg序列）
        if k >= 3 and abs(T[k, 3] - T[k-1, 3]) < tol:
            return T[k, 3], T[:k+1, :k+1], sequences

    return T[max_iter-1, 3], T, sequences


def print_T_table(T_table, title="T数表"):
    """打印T数表到终端"""
    print(f"\n{'='*80}")
    print(f"{title:^80}")
    print(f"{'='*80}")
    print(f"{'k':<5}", end="")
    for j in range(T_table.shape[1]):
        if j == 0:
            print(f"{'T(k,0) 梯形':<20}", end="")
        elif j == 1:
            print(f"{'T(k,1) Simpson':<20}", end="")
        elif j == 2:
            print(f"{'T(k,2) 柯特斯':<20}", end="")
        elif j == 3:
            print(f"{'T(k,3) Romberg':<20}", end="")
        else:
            print(f"{'T(k,' + str(j) + ')':<20}", end="")
    print()
    print("-" * 80)

    for i in range(T_table.shape[0]):
        print(f"{i:<5}", end="")
        for j in range(T_table.shape[1]):
            if j <= i:
                print(f"{T_table[i, j]:<20.12f}", end="")
            else:
                print(f"{'—':<20}", end="")
        print()
    print("=" * 80)


def plot_convergence(sequences_list, labels, exact_values, filename):
    """绘制收敛性比较图"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Romberg算法收敛性分析', fontsize=16, fontweight='bold')

    sequence_names = ['Trapezoidal', 'Simpson', 'Cotes', 'Romberg']
    sequence_labels = ['梯形序列', 'Simpson序列', '柯特斯序列', 'Romberg序列']

    for idx, (seq_name, seq_label) in enumerate(zip(sequence_names, sequence_labels)):
        ax = axes[idx // 2, idx % 2]

        for i, (sequences, label, exact) in enumerate(zip(sequences_list, labels, exact_values)):
            if seq_name in sequences and len(sequences[seq_name]) > 0:
                seq = sequences[seq_name]
                errors = [abs(val - exact) for val in seq]
                iterations = list(range(len(seq)))
                ax.semilogy(iterations, errors, marker='o', label=label, linewidth=2, markersize=6)

        ax.set_xlabel('迭代次数 k', fontsize=11)
        ax.set_ylabel('绝对误差 |I - I_k|', fontsize=11)
        ax.set_title(seq_label, fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"\n可视化图表已保存: {filename}")


# 定义被积函数
def f1(x):
    """积分(1): x^2 * e^(-x^2)"""
    return x**2 * np.exp(-x**2)


def f2(x):
    """积分(2): cot(x)"""
    return 1 / np.tan(x)


if __name__ == "__main__":
    print("\n" + "="*80)
    print("Romberg积分算法 - 计算方法第六周作业".center(80))
    print("="*80)

    # 积分(1): ∫₀² x²e^(-x²)dx
    print("\n\n【积分(1)】 I1 = ∫[0,2] x^2*e^(-x^2)dx")
    result1, T_table1, sequences1 = romberg_integration(f1, 0, 2, tol=1e-10)
    print_T_table(T_table1, "积分(1) T数表")
    print(f"\n最终结果: I1 = {result1:.12f}")

    # 积分(2): ∫_{π/2}^{3π/4} cot(x)dx
    print("\n\n【积分(2)】 I2 = ∫[π/2, 3π/4] cot(x)dx")
    a2, b2 = np.pi/2, 3*np.pi/4
    result2, T_table2, sequences2 = romberg_integration(f2, a2, b2, tol=1e-10)
    print_T_table(T_table2, "积分(2) T数表")
    print(f"\n最终结果: I2 = {result2:.12f}")

    # 理论值（用于误差分析）
    # 积分(1)的理论值可以通过高精度数值方法或符号计算得到
    from scipy import integrate
    exact1, _ = integrate.quad(f1, 0, 2)
    exact2, _ = integrate.quad(f2, np.pi/2, 3*np.pi/4)

    print(f"\n\n{'='*80}")
    print("结果验证（与scipy.integrate.quad比较）".center(80))
    print("="*80)
    print(f"积分(1): Romberg = {result1:.12f}, scipy = {exact1:.12f}, 误差 = {abs(result1-exact1):.2e}")
    print(f"积分(2): Romberg = {result2:.12f}, scipy = {exact2:.12f}, 误差 = {abs(result2-exact2):.2e}")
    print("="*80)

    # 生成可视化图表
    plot_convergence(
        [sequences1, sequences2],
        ['积分(1): ∫[0,2] x^2*e^(-x^2)dx', '积分(2): ∫[π/2, 3π/4] cot(x)dx'],
        [exact1, exact2],
        r'E:\_Materials_in_NWPU_\2026_Spring\course\Computational Methods H\week6\cm_week6\romberg_convergence.png'
    )

    print("\n程序执行完成！\n")
