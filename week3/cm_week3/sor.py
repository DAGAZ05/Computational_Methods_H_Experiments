"""
SOR逐次超松弛迭代法模块
实现SOR方法求解线性方程组
"""
import numpy as np


def sor_iteration(A, b, omega=1.5, tol=1e-10, max_iter=10000,
                  record_history=False):
    """
    SOR逐次超松弛迭代法求解线性方程组 Ax = b

    返回:
        (x, success, iterations, history)
    """
    n = len(b)
    A = A.astype(float)
    b = b.astype(float)

    for i in range(n):
        if abs(A[i, i]) < 1e-15:
            return None, False, 0, []

    x = np.zeros(n)
    history = []

    for k in range(max_iter):
        x_old = x.copy()
        for i in range(n):
            s1 = sum(A[i, j] * x[j] for j in range(i))
            s2 = sum(A[i, j] * x_old[j] for j in range(i, n))
            x[i] = x_old[i] + omega / A[i, i] * (b[i] - s1 - s2)

        if np.any(np.isnan(x)) or np.any(np.isinf(x)):
            return x.copy(), False, k + 1, history

        acc = np.linalg.norm(x - x_old, ord=np.inf)
        if record_history:
            history.append((k + 1, x.copy(), acc))

        if acc < tol:
            return x.copy(), True, k + 1, history

    return x.copy(), False, max_iter, history


def find_optimal_omega(A, b, tol=1e-5, max_iter=10000):
    """两轮搜索寻找最佳松弛因子: 粗搜(0.1) + 细搜(0.01)"""
    best_omega = 1.0
    min_iters = max_iter

    for omega_10 in range(1, 20):
        omega = omega_10 / 10.0
        _, success, iters, _ = sor_iteration(A, b, omega, tol, max_iter)
        if success and iters < min_iters:
            min_iters = iters
            best_omega = omega

    lo = max(0.01, best_omega - 0.15)
    hi = min(1.99, best_omega + 0.15)
    steps = int(round((hi - lo) / 0.01)) + 1
    for i in range(steps):
        omega = lo + i * 0.01
        _, success, iters, _ = sor_iteration(A, b, omega, tol, max_iter)
        if success and iters < min_iters:
            min_iters = iters
            best_omega = omega

    return best_omega, min_iters


def scan_omega(A, b, tol=1e-5, max_iter=10000, step=0.02):
    """全范围扫描omega与迭代次数的关系"""
    omegas, iters_list = [], []
    omega = step
    while omega < 2.0:
        _, success, iters, _ = sor_iteration(A, b, omega, tol, max_iter)
        omegas.append(round(omega, 4))
        iters_list.append(iters if success else None)
        omega += step
    return omegas, iters_list
