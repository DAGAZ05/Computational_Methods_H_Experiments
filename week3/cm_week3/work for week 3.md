分别编写Gauss消去法、全选主元Gauss消去法、Doolittle分解法、Jacobi迭代法、Gauss-Seidel迭代法、SOR（逐次超松弛）方法​程序，求解求解如下四元线性方程组：
$$
\begin{bmatrix}
-4 & 1 & 1 & 1 \\
1 & -4 & 1 & 1 \\
1 & 1 & -4 & 1 \\
1 & 1 & 1 & -4
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
x_3 \\
x_4
\end{bmatrix}
=
\begin{bmatrix}
1 \\
1 \\
1 \\
1
\end{bmatrix}
$$

**要求**：

1. 取初始值$\vec{x}^{(0)}=[0, 0, 0, 0]^T$；
2. 对于迭代法，迭代终止条件：$$\max\limits_{1\leq i\leq 4} |x_i^{(k+1)}-x_i^{(k)}|\le 10^{-5}$$；
3. 对于SOR方法，在$0<\omega<2$范围内，试找出使迭代收敛最快的最佳松弛因子​**$\omega_{opt}$**；
4. 对于和题中方程组形式相同的更高阶方程组，比较各方法执行时间或迭代次数随阶数增加而产生的差异。
