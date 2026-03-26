"""
Gauss-Seidel迭代法模块
实现G-S迭代法求解线性方程组
"""
import numpy as np


def gauss_seidel_iteration(A, b, tol=1e-10, max_iter=10000,
                           record_history=False):
    """
    Gauss-Seidel迭代法求解线性方程组 Ax = b

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
            s2 = sum(A[i, j] * x_old[j] for j in range(i + 1, n))
            x[i] = (b[i] - s1 - s2) / A[i, i]

        if np.any(np.isnan(x)) or np.any(np.isinf(x)):
            return x.copy(), False, k + 1, history

        acc = np.linalg.norm(x - x_old, ord=np.inf)
        if record_history:
            history.append((k + 1, x.copy(), acc))

        if acc < tol:
            return x.copy(), True, k + 1, history

    return x.copy(), False, max_iter, history
