"""
Jacobi迭代法模块
实现Jacobi迭代法求解线性方程组
"""
import numpy as np


def jacobi_iteration(A, b, tol=1e-10, max_iter=10000, record_history=False):
    """
    Jacobi迭代法求解线性方程组 Ax = b

    迭代格式: x_i^{(k+1)} = (b_i - sum_{j!=i} a_ij * x_j^{(k)}) / a_ii

    返回:
        (x, success, iterations, history)
        history: record_history=True时为[(k, x_copy, accuracy), ...]
    """
    n = len(b)
    A = A.astype(float)
    b = b.astype(float)

    for i in range(n):
        if abs(A[i, i]) < 1e-15:
            return None, False, 0, []

    x = np.zeros(n)
    x_new = np.zeros(n)
    history = []

    for k in range(max_iter):
        for i in range(n):
            s = sum(A[i, j] * x[j] for j in range(n) if j != i)
            x_new[i] = (b[i] - s) / A[i, i]

        if np.any(np.isnan(x_new)) or np.any(np.isinf(x_new)):
            return x_new.copy(), False, k + 1, history

        acc = np.linalg.norm(x_new - x, ord=np.inf)
        if record_history:
            history.append((k + 1, x_new.copy(), acc))

        if acc < tol:
            return x_new.copy(), True, k + 1, history

        x = x_new.copy()

    return x_new.copy(), False, max_iter, history
