# 计算方法H 第三章实践作业

求解四元线性方程组 $A\vec{x}=\vec{b}$，分别使用6种方法并比较结果。

## 方程组

$$
\begin{bmatrix}
-4 & 1 & 1 & 1 \\
1 & -4 & 1 & 1 \\
1 & 1 & -4 & 1 \\
1 & 1 & 1 & -4
\end{bmatrix}
\begin{bmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{bmatrix}
=
\begin{bmatrix} 1 \\ 1 \\ 1 \\ 1 \end{bmatrix}
$$

精确解：$\vec{x}^* = [-1, -1, -1, -1]^T$

## 文件结构

| 文件 | 说明 |
|------|------|
| `main.py` | 主脚本（命令行入口、比较表格、可视化） |
| `gauss.py` | Gauss顺序消去法 |
| `gauss_full_pivot.py` | 全主元Gauss消去法 |
| `doolittle.py` | Doolittle分解法 |
| `jacobi.py` | Jacobi迭代法 |
| `gauss_seidel.py` | Gauss-Seidel迭代法 |
| `sor.py` | SOR逐次超松弛迭代法（含最佳ω搜索） |

## 用法

```bash
# 全部方法执行：比较表格 + 高阶扩展 + 可视化
python main.py

# 单独执行某个方法（直接法输出解和执行时间，迭代法输出每步迭代详情）
python main.py gauss            # Gauss顺序消去法
python main.py gauss_full       # 全主元Gauss消去法
python main.py doolittle        # Doolittle分解法
python main.py jacobi           # Jacobi迭代法
python main.py gauss_seidel     # Gauss-Seidel迭代法
python main.py sor              # SOR逐次超松弛迭代法
```

## 作业要求

1. 初始向量 $\vec{x}^{(0)} = [0, 0, 0, 0]^T$
2. 迭代终止条件：$\max_{1 \le i \le 4} |x_i^{(k+1)} - x_i^{(k)}| \le 10^{-5}$
3. SOR最佳松弛因子 $\omega_{opt}$ 通过试算法在 $0 < \omega < 2$ 范围内搜索
4. 高阶同形方程组（对角元 $-n$，非对角元 $1$，$n=4\sim128$）性能扩展比较

## 可视化输出

`comparison.png` 包含三个子图：
- (a) 迭代法收敛曲线（精度随迭代次数变化）
- (b) SOR松弛因子 $\omega$ 与迭代次数的关系
- (c) 高阶方程组各方法执行时间随阶数的变化

## 依赖

- Python 3
- NumPy
- Matplotlib（仅全部执行时需要）
