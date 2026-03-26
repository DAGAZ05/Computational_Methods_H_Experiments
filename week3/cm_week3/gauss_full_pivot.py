"""
全主元Gauss消去法模块
实现全主元Gauss消去法求解线性方程组
"""
import numpy as np


def gauss_full_pivot(A, b):
    """
    全主元Gauss消去法求解线性方程组 Ax = b

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

    # 记录列交换信息
    col_order = list(range(n))

    # 消元过程
    for k in range(n - 1):
        # 在k~n行、k~n列中选取绝对值最大的元素
        max_val = 0
        max_i, max_j = k, k
        for i in range(k, n):
            for j in range(k, n):
                if abs(A[i, j]) > max_val:
                    max_val = abs(A[i, j])
                    max_i, max_j = i, j

        if max_val < 1e-15:
            return None, False

        # 行交换
        if max_i != k:
            A[[k, max_i], :] = A[[max_i, k], :]
            b[k], b[max_i] = b[max_i], b[k]

        # 列交换
        if max_j != k:
            A[:, [k, max_j]] = A[:, [max_j, k]]
            col_order[k], col_order[max_j] = col_order[max_j], col_order[k]

        # 计算消元因子并消元
        for i in range(k + 1, n):
            l_ik = A[i, k] / A[k, k]
            A[i, k:n] = A[i, k:n] - l_ik * A[k, k:n]
            b[i] = b[i] - l_ik * b[k]

    # 回代过程
    y = np.zeros(n)
    for k in range(n - 1, -1, -1):
        if abs(A[k, k]) < 1e-15:
            return None, False
        y[k] = (b[k] - np.dot(A[k, k+1:n], y[k+1:n])) / A[k, k]

    # 恢复未知量次序
    x = np.zeros(n)
    for i in range(n):
        x[col_order[i]] = y[i]

    return x, True
