# 第一章 实践题目

1．设$S_n=\sum_{j=2}^{n}\frac{1}{j^2-1}$，其精确值为$S_n=\frac{1}{2}\left(\frac{3}{2}-\frac{1}{n}-\frac{1}{n+1}\right)$。  

(1) 编制计算$S_n=\frac{1}{2^2-1}+\frac{1}{3^2-1}+\dots+\frac{1}{n^2-1}$的通用程序。  

(2) 编制计算$S_n=\frac{1}{n^2-1}+\frac{1}{(n-1)^2-1}+\dots+\frac{1}{2^2-1}$的通用程序。  

(3) 按两种顺序分别计算$n=10^k, k=1,\dots,10$时，计算上述量的结果，并指出有效位数(用单精度浮点数计算)。  

(4) 通过本上机题，你明白了什么？

2．选取$\varepsilon=10^{-2k}, k=0,1,\dots,10$，完成如下问题：  

(1) 用Gauss消去法求解方程组：
$$
\begin{bmatrix}\varepsilon&1\\1&1\end{bmatrix}\begin{bmatrix}x_1\\x_2\end{bmatrix}=\begin{bmatrix}1+\varepsilon\\2\end{bmatrix}
$$

(2) 用Gauss消去法求解方程组：

$$
\begin{bmatrix}1&1\\\varepsilon&1\end{bmatrix}\begin{bmatrix}x_2\\x_1\end{bmatrix}=\begin{bmatrix}2\\1+\varepsilon\end{bmatrix}
$$

(3) 随着$\varepsilon$的逐渐减小，问题的结果发生了什么改变？你明白了什么？

3．编写程序，用无穷级数

$$
e^x=1+x+\frac{x^2}{2!}+\frac{x^3}{3!}+\dots
$$

计算以下问题：  

(1) 用$x=1,5,10,15,20$测试程序，观察程序运行结果。  

(2) 用$x=-1,-5,-10,-15,-20$测试程序，观察程序运行结果。  

(3) 当$x<0$时，改用数值计算公式：
$$
e^{-x}=\frac{1}{e^x}=\frac{1}{1+x+\frac{x^2}{2!}+\frac{x^3}{3!}+\dots}
$$

计算$x=-1,-5,-10,-15,-20$时的近似值，数值实验中又有怎样的现象出现？你明白了什么？

4．对$n=0,1,2,\dots,20$计算定积分$y_n=\int_{0}^{1}x^n e^{x-1}dx$的近似值。  

(1) 利用递推公式：$y_n=1-n y_{n-1}$，$n=0,1,2,\dots,20$进行计算，观察计算结果。  

(2) 利用递推公式：$y_{n-1}=\frac{1-y_n}{n}$，$n=20,19,18,\dots,0$，注意到$0<y_n<\int_{0}^{1}x^n dx=\frac{1}{n+1}$，不妨取$y_{20}=0$进行计算，观察计算结果。

‍
