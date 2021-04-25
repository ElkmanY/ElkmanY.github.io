---
layout: post
title: "粒子群优化及其Matlab实现"
categories: techo
# author: elkman
---

#### 1 粒子群优化简介

考虑最优化问题
$$ 
\min_{\boldsymbol{x}\in \Omega}{J(\boldsymbol{x})} \tag{1}
$$ 
满足约束
$$ 
\underline{\boldsymbol{x}} \leqslant \boldsymbol{x} \leqslant \overline{\boldsymbol{x}} \tag{2}
$$ 
其中，$\underline{\boldsymbol{x}}$ 和 $\overline{\boldsymbol{x}}$ 分别为设计变量 $\boldsymbol{x} \in \mathbb{R}^d$ 的下界和上界，$\Omega$ 为搜素空间
$$\Omega = \{ \boldsymbol{x}|\underline{\boldsymbol{x}} \leqslant \boldsymbol{x} \leqslant \overline{\boldsymbol{x}} \}. \tag{3}$$

在粒子群优化（PSO）中，每一个可行解被称为 “粒子”，粒子群中的每个粒子代表了 $d$ 维搜索空间中的一个点。
在一个有 $n$ 个粒子的粒子群中，第 $i$ 个例子的位置为
$$ \boldsymbol{x}_i = [ x_{i1}, x_{i2}, ..., x_{id}]^\text{T} \tag{4}$$ 

粒子通过不断迭代更新自己的位置来搜索最优解，粒子的轨迹满足运动学方程：
$$ \boldsymbol{x}_i(t+1) = \boldsymbol{x}_i(t) + \boldsymbol{v}_i(t+1) \tag{5}$$式中，$t$ 为算法当前的迭代次数，而 $\boldsymbol{v}_i \in \mathbb{R}^d$ 是第 $i$ 个粒子的速度，向量的元素是该粒子在 $d$ 个维度上速度的分量。

速度矢量控制着粒子在搜索空间中移动的方式，第 $i$ 个粒子的速度定义为:
$$
 \boldsymbol{v}_i(t+1) = \omega \boldsymbol{v}_i(t) + c_1( \boldsymbol{p}_i -\boldsymbol{x}_i(t) ) \odot \boldsymbol{r}_1 + c_2( \boldsymbol{g} -\boldsymbol{x}_i(t) ) \odot \boldsymbol{r}_2 \tag{6}
$$ 速度由三个项组成:第一项定义了惯性或动量，通过保持先前的流动方向 $\boldsymbol{v}_i(t)$，防止粒子剧烈改变方向;第二项称为自我认知部分，它代表了粒子有回到自己曾经找到的最佳位置 $\boldsymbol{p}_i \in \mathbb{R}^d$ 的倾向;最后一项称为社会部分，它代表一个粒子向整个粒子群的当前最佳位置移动 $\boldsymbol{g} \in \mathbb{R}^d$ 的倾向。
另外，加速常数 $c_1$ 和 $c_2$ 为 $[0, 4]$ 之间的实数，数值越大，代表该项占的比重越大。
惯性权重 $\omega $ 为$[0.4, 0.8]$ 之间的实数，数值越大，粒子改变方向的趋势越小。
$\boldsymbol{r}_1 \in \mathbb{R}^d$ 和 $\boldsymbol{r}_2 \in \mathbb{R}^d$ 为服从均匀分布的随机数向量，而 $\odot$ 为向量或矩阵对应位置元素相乘的乘法运算。

#### 2 Matlab代码实现
本节将一步一步地阐述粒子群优化算法的Maltab代码实现。

在Matlab中进行矩阵化的运算而非for循环，这样会节省很多执行时间。于是以下操作均是以矩阵为单位。
根据式 $(4)$ 定义粒子群
$$ 
\boldsymbol{X}=[\boldsymbol{x}_1, \boldsymbol{x}_2, ..., \boldsymbol{x}_n] = 
\left[ 
\begin{matrix}
x_{11} & x_{21} & \ldots & x_{n1} \\ 
x_{12} & x_{22} & \ldots & x_{n2} \\ 
\vdots & \vdots & \ddots & \vdots  \\ 
x_{1d} & x_{2d} & \ldots & x_{nd}  
\end{matrix} \right]
$$ 矩阵中的每一列为一个粒子，该列元素为该粒子在搜索空间 $\Omega$ 中的坐标。
同样地，粒子定义粒子群的速度矩阵
$$ 
\boldsymbol{V}=[\boldsymbol{v}_1, \boldsymbol{v}_2, ..., \boldsymbol{v}_n] = 
\left[ 
\begin{matrix}
v_{11} & v_{21} & \ldots & v_{n1} \\ 
v_{12} & v_{22} & \ldots & v_{n2} \\ 
\vdots & \vdots & \ddots & \vdots  \\ 
v_{1d} & v_{2d} & \ldots & v_{nd}  
\end{matrix} \right]
$$ 
所有粒子个体的曾找到的最优位置
$$ 
\boldsymbol{P}=[\boldsymbol{p}_1, \boldsymbol{p}_2, ..., \boldsymbol{p}_n] = 
\left[ 
\begin{matrix}
p_{11} & p_{21} & \ldots & p_{n1} \\ 
p_{12} & p_{22} & \ldots & p_{n2} \\ 
\vdots & \vdots & \ddots & \vdots  \\ 
p_{1d} & p_{2d} & \ldots & p_{nd}  
\end{matrix} \right]
$$ 

##### Step 1
在开搜索之前，对粒子群各粒子的位置和速度进行初始化操作
$$ 
\boldsymbol{X}(0) = \text{rand}(d,n) \odot ( \overline{\boldsymbol{x}} - \underline{\boldsymbol{x}} ) + \underline{\boldsymbol{x}} \\
\boldsymbol{V}(0) = \text{rand}(d,n) \odot ( \overline{\boldsymbol{v}} - \underline{\boldsymbol{v}} ) + \underline{\boldsymbol{v}} 
$$
```matlab
X = rand(d,n).*(xp-xb) + xb;  
V = rand(d,n).*(vp-vb) + vb;  
```
以便后续的迭代，初始化所有粒子个体的曾找到的最优位置，即粒子初始位置
$$
\boldsymbol{P}(0) = \boldsymbol{X}(0)
$$
```matlab
P = X;
F = J(X);
Fp = F;
```
同时初始化粒子群当前的全局最优位置，即在初始位置中的一个使目标函数值最小的位置
$$
\boldsymbol{g}(0) = \argmin_{\boldsymbol{x}} J(\boldsymbol{X}(0))
$$
```matlab
F = J(X);
[fg(1),ig] = min(F);
g(:,1) = X(:,ig);
```
##### Step 2
粒子群优化算法的迭代用离散时间序列 $[0,1,...,T]$ 来表示，并通过for循环来实现，即从 $t=0$ 时刻开始一直到 $t=T$ 时刻结束：
```matlab
for t = 1:T
    % 粒子群优化算法的迭代搜索
end
```

在迭代过程中，首先根据式 $(6)、(5)$ 计算下一时刻粒子群的速度和位置
$$
 \boldsymbol{V}(t+1) = 
 \omega \boldsymbol{V}(t) + c_1( \boldsymbol{P}(t) -\boldsymbol{X}(t) ) \odot \boldsymbol{R}_1 + c_2( \boldsymbol{g}(t) -\boldsymbol{X}(t) ) \odot \boldsymbol{R}_2 \\
 \boldsymbol{X}(t+1) = \boldsymbol{X}(t) + \boldsymbol{V}(t+1)
$$
```matlab
V_ = V + w*c1*(P-X).*rand(d,n) + c2*(g(:,t)-X).*rand(d,n);
X_ = X + V_;  
```
##### Step 3
随后根据粒子群新的位置，更新个体的曾找到的最优位置。
即与之前的个体最优位置相比较，如果当前位置较优秀，则个体的当前位置替换个体的历史最优位置，否则不替换：
$$
\boldsymbol{P}(t+1) = \argmin_{\boldsymbol{p},\boldsymbol{x}} \{J(\boldsymbol{P}(t)),J(\boldsymbol{X}(t+1))\}
$$
```matlab
F_ = J(X_);
flag = F_ < Fp;
P_ = X_.*flag + P.*(~flag);
Fp_ = F_.*flag + Pp.*(~flag);
```
随后更新粒子群当前的全局最优位置，即在更新后的个体最优位置中找出一个最优解，并与全局最优解相比得到下一时刻的全局最优解：
$$
\boldsymbol{g}(t+1) = \argmin_{\boldsymbol{p}} \{J(\boldsymbol{P}(t+1)),\boldsymbol{g}(t)\}
$$
```matlab
Fpa = [Fp_,fg(t)];
Pa = [P_,g(:,t)];
[fg(t+1),ig] = min(Fpa);
g(:,t+1) = Pa(:,ig);
```
迭代的最后完成各变量的更新
```matlab
V = V_;
X = X_;
P = P_;
Fp = Fp_;
```

#### 3 结语
文中代码来源：*[ElkmanY/pso](https://github.com/ElkmanY/pso)*
本文参考文献：*[Particle swarm optimization (PSO). A tutorial](https://www.sciencedirect.com/science/article/pii/S0169743915002117)*
