"""
Gauss顺序消去法模块
实现Gauss消去法求解线性方程组
"""
import numpy as np


def gauss_elimination(A, b):
    """
    Gauss顺序消去法求解线性方程组 Ax = b

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

    # 消元过程
    for k in range(n - 1):
        # 检查主元是否为零
        if abs(A[k, k]) < 1e-15:
            return None, False

        # 计算消元因子并消元
        for i in range(k + 1, n):
            l_ik = A[i, k] / A[k, k]
            A[i, k:n] = A[i, k:n] - l_ik * A[k, k:n]
            b[i] = b[i] - l_ik * b[k]

    # 回代过程
    x = np.zeros(n)
    for k in range(n - 1, -1, -1):
        if abs(A[k, k]) < 1e-15:
            return None, False
        x[k] = (b[k] - np.dot(A[k, k+1:n], x[k+1:n])) / A[k, k]

    return x, True
