# 计算方法第二周作业

## 问题描述

求解方程 **4cos(x) = e^x** 的根，精度要求 ε = 10^-4

## 实现方法

本程序实现了四种数值迭代方法：

1. **简单迭代法** (Simple Iteration)
   - 初值：x₀ = π/4
   - 迭代函数：使用 x = arccos(e^x / 4)

2. **斯蒂芬森迭代法** (Steffensen's Method)
   - 初值：x₀ = π/4
   - 加速收敛的迭代方法

3. **Newton迭代法** (Newton's Method)
   - 初值：x₀ = π/4
   - 使用函数及其导数

4. **双点弦截法** (Secant Method)
   - 初值：x₀ = π/4, x₁ = π/2
   - 不需要计算导数的割线法

## 使用方法

### 基本用法

```bash
# 运行单个方法
python solver.py simple      # 简单迭代法
python solver.py steffensen  # 斯蒂芬森迭代法
python solver.py newton      # Newton迭代法
python solver.py secant      # 双点弦截法

# 运行所有方法并生成可视化
python solver.py all
```

### 显示迭代详情

使用 `--verbose` 参数可以显示每次迭代的详细信息，包括：
- k：迭代次数
- x_k：当前迭代值
- |x_k - x_(k-1)|：相邻两次迭代的差值

```bash
# 显示单个方法的迭代详情
python solver.py newton --verbose

# 显示所有方法的迭代详情
python solver.py all --verbose
```

### 自定义精度

```bash
# 自定义精度
python solver.py all --epsilon 1e-6

# 自定义精度并显示迭代详情
python solver.py all --epsilon 1e-6 --verbose
```

### 查看帮助

```bash
python solver.py --help
```

## 输出结果

### 1. 控制台输出

运行 `python solver.py all` 会输出：
- 每种方法的根
- 迭代次数
- 计算时间
- 函数值误差 |f(x)|
- 效率比较（迭代次数比和时间比）

**示例输出：**
```
============================================================
Solving equation: 4cos(x) = e^x
Tolerance: epsilon = 0.0001
============================================================

--- NEWTON METHOD ---

Newton's Method:
  Root: x = 0.9047882181
  Iterations: 3
  Time: 0.0311 ms
  |f(x)|: 1.19e-09

============================================================
EFFICIENCY COMPARISON:
============================================================
Simple Iteration    : Iter ratio = 11.00x, Time ratio = 4.74x
Steffensen's Method : Iter ratio = 1.00x, Time ratio = 15.99x
Newton's Method     : Iter ratio = 1.00x, Time ratio = 1.00x
Secant Method       : Iter ratio = 1.33x, Time ratio = 1.36x
```

### 2. 迭代详情输出（使用 --verbose）

使用 `--verbose` 参数时，会显示每次迭代的详细信息：

```
--- NEWTON METHOD ---
k     x_k             |x_k - x_(k-1)|
----------------------------------------
0     0.7853981634    -
1     0.9118784721    1.2648030875e-01
2     0.9048101871    7.0682850008e-03
3     0.9047882181    2.1969061469e-05
```

### 3. 可视化图表（results.png）

仅在运行 `python solver.py all` 时生成，包含：
- 函数图像和根的位置
- 迭代次数比较
- 计算时间比较
- 收敛过程（误差变化，对数坐标）
- 迭代点演化

## 依赖库

```bash
pip install numpy matplotlib
```

## 结果分析

根据测试结果（ε = 10^-4）：

| 方法 | 迭代次数 | 根 x | 相对效率 |
|------|---------|------|---------|
| 简单迭代法 | 33 | 0.9048280823 | 最慢（11倍） |
| 斯蒂芬森迭代法 | 3 | 0.9047882179 | 快速 |
| Newton迭代法 | 3 | 0.9047882181 | 最快 |
| 双点弦截法 | 4 | 0.9047880040 | 较快 |

**主要结论：**
- **Newton迭代法**：收敛最快，计算效率最高（3次迭代）
- **斯蒂芬森迭代法**：与Newton法迭代次数相同，但每次迭代需要计算两次φ(x)，时间略长
- **双点弦截法**：不需要导数，收敛速度略慢于Newton法（4次迭代）
- **简单迭代法**：收敛最慢，需要33次迭代，效率最低

**收敛条件：** 所有方法均使用 |x_k - x_(k-1)| < ε 作为终止条件

## 文件说明

- `solver.py` - 主程序文件（全英文代码和注释）
- `results.png` - 可视化结果图（运行 `python solver.py all` 后生成）
- `README.md` - 本说明文件

## 程序特性

1. **命令行参数支持**：可选择单个方法或全部执行
2. **迭代详情输出**：使用 `--verbose` 显示每次迭代的 k, x_k, |x_k - x_(k-1)|
3. **科学风格可视化**：符合学术论文标准的图表样式
4. **精度可调**：通过 `--epsilon` 参数自定义收敛精度
5. **效率比较**：自动计算并显示各方法的相对效率

## 技术说明

- **迭代终止条件**：|x_k - x_(k-1)| < ε
- **简单迭代函数**：x = arccos(e^x / 4)
- **可视化字体**：已修复负号显示问题（设置 `axes.unicode_minus: False`）
- **代码语言**：Python脚本全部使用英文
