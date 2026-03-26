"""
Doolittle分解法模块
实现Doolittle三角分解法求解线性方程组
"""
import numpy as np


def doolittle_decomposition(A, b):
    """
    Doolittle分解法求解线性方程组 Ax = b
    A = LU，其中L为单位下三角矩阵，U为上三角矩阵

    参数:
        A: 系数矩阵 (n×n)
        b: 右端向量 (n×1)

    返回:
        x: 解向量
        success: 是否成功求解
    """
    n = len(b)
    A = A.astype(float).copy()
    b = b.astype(float).copy()

    L = np.zeros((n, n))
    U = np.zeros((n, n))

    # Doolittle分解
    for k in range(n):
        # 计算U的第k行
        for j in range(k, n):
            U[k, j] = A[k, j] - sum(L[k, m] * U[m, j] for m in range(k))

        # 检查对角元素
        if abs(U[k, k]) < 1e-15:
            return None, False

        # L的对角元素为1
        L[k, k] = 1.0

        # 计算L的第k列
        for i in range(k + 1, n):
            L[i, k] = (A[i, k] - sum(L[i, m] * U[m, k]
                       for m in range(k))) / U[k, k]

    # 前推求解 Ly = b
    y = np.zeros(n)
    for k in range(n):
        y[k] = b[k] - sum(L[k, m] * y[m] for m in range(k))

    # 回代求解 Ux = y
    x = np.zeros(n)
    for k in range(n - 1, -1, -1):
        x[k] = (y[k] - sum(U[k, m] * x[m] for m in range(k + 1, n))) / U[k, k]

    return x, True
