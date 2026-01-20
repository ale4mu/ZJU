# Help 目录说明

`help` 目录包含了本研究项目中使用的各种 MATLAB 函数和脚本，主要用于实现不确定性量化 (UQ) 方法，特别是随机配点法 (Stochastic Collocation, SC) 和随机伽辽金法 (Stochastic Galerkin, SG)，以及一些通用的数值工具。

## 目录结构

```
help/
├── SC/             # 随机配点法相关函数
├── SG/             # gPC-SG法相关函数
└── *.m             # 通用工具函数
```

## SC 子目录

此子目录包含用于实现随机配点法的函数，主要用于计算参考解。

*   `Euler_forward.m`: 欧拉法求解器，替代 RK3 求解器，主要用于 Example 4.5。
*   `RK3_SC.m`: 三阶 Runge-Kutta 方法求解器。
*   `WENO.m`: WENO 格式实现(全局 Lax 分裂)。
*   `WENO_LLF.m`: WENO 格式实现(局部 Lax 分裂)。
*   `colloaction41.m`: 针对 Example 4.1 的随机配点法求解器。
*   `collocation42.m`: 针对 Example 4.2 的随机配点法求解器。
*   `collocation43.m`: 针对 Example 4.3 的随机配点法求解器。
*   `collocation45.m`: 针对 Example 4.5 的随机配点法求解器。
*   `solve_2D_45.m`: 针对 Example 4.5 的二维问题求解器(用于 collocation45.m,使用欧拉法求解器)。
*   `solve_SC.m`: 随机配点法的通用求解框架(使用RK3求解器)。
*   `upwind_SC43.m`: 针对 Example 4.3 的迎风格式求解器(用于 collocation43.m)。
*   `upwind_SC45.m`: 针对 Example 4.5 的迎风格式求解器(用于 Euler_forward.m)。

## SG 子目录

此子目录包含用于实现 gPC-SG 的函数。

*   `RHS_SG.m`: 关键函数，论文公式 (2.8) 右端项的计算。
*   `RK3_SG.m`: 三阶 Runge-Kutta 方法求解器。
*   `WENO_gPC.m`: WENO 格式实现(无限制器版本,用于 Example 4.1)。
*   `WENO_gPC43.m`: WENO 格式实现(有限制器版本,用于 Example 4.3)。
*   `WENO_gPC_2D.m`: 针对二维 gPC-SG 系统的 WENO 格式实现(有限制器版本,用于 Example 4.5)。
*   `gPC_SG4_2.m`: 针对 Example 4.2 的 gPC-SG 实现。
*   `gpc_SG4_1.m`: 针对 Example 4.1 的 gPC-SG 实现。
*   `gpc_SG4_3.m`: 针对 Example 4.3 的 gPC-SG 实现。
*   `gpc_SG_2D.m`: 针对二维问题的 gPC-SG 实现(Example 4.5)。

## 通用工具函数

这些函数直接位于 `help` 目录下。

*   `Euler_Eigen.m`: 欧拉方程的一维特征值、Jacobian 矩阵和左特征向量计算。
*   `Euler_Eigen_2D.m`: 欧拉方程的二维特征值计算。
*   `Euler_Eigen_Batch.m`: 批量计算欧拉方程的一维特征值。
*   `Euler_Eigen_Batch_2D.m`: 批量计算欧拉方程的二维特征值。
*   `L2Error.m`: 计算 L2 误差。
*   `Limiter.m`: 限制器Limiter，用于数值格式的稳定性。
*   `LinfError.m`: 计算 L 无穷误差。
*   `WENO_neg.m`: 5阶精度 WENO 重构, WENO 格式中的负通量计算。
*   `WENO_pos.m`: 5阶精度 WENO 重构, WENO 格式中的正通量计算。
*   `build_2D_gPC.m`: 构建两个随机参数 gPC 系统所需的 Gauss 积分点和权重等，使用总次数基。
*   `check_Euler.m`: 检查一维欧拉方程的物理量是否在容许状态集。
*   `check_Euler_2D.m`: 检查二维欧拉方程的物理量是否在容许状态集。
*   `f_direction.m`: 用于二维欧拉方程的算子分裂。
*   `hyper.m`: 一维双曲性检查。
*   `hyper_2D.m`: 二维双曲性检查。
*   `lgwt.m`: Gauss - Legendre 权重和节点计算。
*   `my_legendre.m`: Legendre 多项式计算。
*   `order.m`: 计算收敛阶。
