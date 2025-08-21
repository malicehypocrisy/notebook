# 一、分子血缘相关矩阵$\mathbf{A}$及其逆矩阵$\mathbf{A}^{-1}$的计算

## 1.$\mathbf{A}$的计算方法

### &#9312;通径系数

$$
a_{ii}=1+f_{i}\\
a_{ij}=\sum{(\frac{1}{2})^{n+n^{`}}(1+f_{A})}
$$

其中$f_{i}$为个体$i$的近交系数（为父母近交系数的一半），$n$和$n^`$是个体$i$和$j$经过不同通径到共同祖先$A$的时代数，$f_A$为$A$的近交系数。

### &#9313;循环公式

$$
a_{ii} = 1 + \frac{1}{2} a_{\text{父母}} \quad \text{（对角线元素）}
\\
a_{ij} = \frac{1}{2} \left( a_{i(j\text{父})} + a_{i(j\text{母})} \right) = a_{ji}
$$

> [!NOTE]
>
> **若个体$j$的父亲或母亲或双亲未知时，则$a_{ii}=1，a_{ij}=0$**

| 个体 |  父  |  母  |
| :--: | :--: | :--: |
|  1   |  -   |  -   |
|  2   |  -   |  -   |
|  3   |  1   |      |
|  4   |  1   |  2   |
|  5   |  3   |  4   |
|  6   |  1   |  4   |
|  7   |  5   |  6   |

列表注意：

&#9332;在个体一列中应包含所有在父和母列出现过的个体

&#9333;在个体一列保证后代不会出现在亲代前，一般按出生日期排序，先出生在前

&#9334;个体编号从1开始练习编号

## 2.$\mathbf{A}^{-1}$的计算方法

### &#9312;HENDERSON法

$$
\mathbf{A}^{-1}=(\mathbf{L}\mathbf{L}^{T})^{-1}=(\mathbf{TDD}\mathbf{T}^{T})^{-1}=\mathbf{T}^{T^{-1}}\mathbf{D}^{-2}\mathbf{T}^{-1}
$$

$\mathbf{T}^{-1}$（**下三角矩阵**）的对角线元素全为$1$，在其第$i$行个体的每个亲本对应的元素为-0.5，其余元素为零。

$\mathbf{D}^{-2}$是对角线矩阵$\delta_i$是其对角线上第$i$个元素
$$
\delta_i = \frac{1}{l_{ii}^2}\tag{1}
$$
**$l_{ii}$是$\mathbf{L}$阵的第$i$个对角线元素，$f$近交系数$f_i=a_{ii}-1$,$a_{ii}$是$\mathbf{A}$的对角线元素。**
$$
\\
l_{ii}=
\begin{cases}
[0.5 - 0.25(f_s + f_d)]^{0.5} & \text{当个体 i 的双亲 s 和 d 已知} \\
(0.75 - 0.25f_p)^{0.5} & \text{当个体 i 的一个亲体 p 已知} \\
1 & \text{当个体 i 双亲未知}\tag{2}
\end{cases}
$$

### &#9313;QUAAS法（简化$l_{ii}$的求解）

#### &#9372;**求解$l_{ii}$的步骤**：

&#9332;将所有个体按照计算$\mathbf{A}$阵时方式排列。

&#9333;建立两个阶数为n的零向量$\mathbf{v}$和$\mathbf{u}$,$\mathbf{v}$存放$l_{ii}$，$\mathbf{u}$用于存放$\mathbf{u}_i$,i=1,……，n,n为个体总数。

&#9334;对于i=1,……，n,计算

* $\mathbf{v}_i=l_{ii}(使用公式3)$

$$
u_i = \sum_{k=1}^i l_{ik}^2=a_{ii}
\\
l_{ii} =
\begin{cases}
\sqrt{1 - 0.25(u_s + u_d)} & \text{双亲 } s \text{ 和 } d \text{ 已知} \\
\sqrt{1 - 0.25u_p} & \text{单亲 } p \text{ 已知} \\
1 & \text{双亲未知}
\end{cases}\tag{3}
\\

\text{其中} (u_s, u_d, u_p) \text{为父母个体的累计值，反映遗传贡献}。
$$

* $\mathbf{u}_i=u_i+v_i^2$
* 对于i=1,……，n计算,

$$
l_{ki} =
\begin{cases}
\frac{l_{si} + l_{di}}{2} & \text{如个体 \(k\) 的双亲 \(s\) 和 \(d\) 排在 \(i\) 以后} \\
\frac{l_{pi}}{2} & \text{如个体 \(k\) 的一个亲体 \(p\) 排在 \(i\) 以后} \\
0 & \text{如个体 \(k\) 的双亲均排在 \(i\) 之前}
\end{cases}

\\

u_k = u_k + l_{ki}^2
$$

#### &#9373;求$\mathbf{A}^{-1}$的步骤

&#9332;将所有个按照$\mathbf{A}$求阵时的方式排列，

&#9333;将$\mathbf{A}^{-1}$中的所有元素置为零

&#9334;对于i=1,……，n,计算

* 按上述方法计算$l_{ii}$
* 如i的双亲s 和 d 已知，将下列数值加到$\mathbf{A}^{-1}$中：

$$
\begin{aligned}
\delta_i &\rightarrow a_{ii} \\
-\delta_i / 2 &\rightarrow a_{is}, a_{si}, a_{id}, a_{di} \\
\delta_i / 4 &\rightarrow a_{ss}, a_{dd}, a_{sd}, a_{ds}
\end{aligned}
$$

* 如i的一个亲本p已知，则
	$$
	\begin{aligned}
	\delta_i &\rightarrow a_{ii} \\
	-\delta_i / 2 &\rightarrow a_{ip}, a_{pi} \\
	\delta_i / 4 &\rightarrow a_{pp}
	\end{aligned}
	$$
	
* 如i的双亲均未知，则
	$$
	\delta_i \rightarrow a_{ii}
	$$


其中 $\delta_i$ 由式 (1) 计算，但当**群体为非近交群体时**，HENDERSON 证明：

$$
\delta_i =
\begin{cases}
2 & \text{当个体 } i \text{ 的双亲已知} \\
\frac{4}{3} & \text{当个体 } i \text{ 的一个亲体已知} \\
1 & \text{当个体 } i \text{ 双亲未知}
\end{cases}
$$

## 3.非基础群血缘关系逆矩阵

$$
\mathbf{A} =
\begin{bmatrix}
\mathbf{A}_{11} & \mathbf{A}_{12} \\
\mathbf{A}_{12}^\prime & \mathbf{A}_{22}
\end{bmatrix}
\\
\mathbf{A}^{-1}  =
\begin{bmatrix}
\mathbf{B}_{11} & \mathbf{B}_{12} \\
\mathbf{B}_{12}^\prime & \mathbf{B}_{22}
\end{bmatrix}
\\
$$

其中$\mathbf{A}_{11}$ 和 $\mathbf{A}_{22}$ 分别对应的是基础群和非基础群的个体的血缘相关子矩阵，我们需要的是 $\mathbf{A}_{22}^{-1}$。
$$
\mathbf{A}_{22}^{-1}=\mathbf{B}_{22}-\mathbf{B}_{12}^{T}\mathbf{B}_{11}^{-1}\mathbf{B}_{12}
$$

# 二、吸收法

模型：$\mathbf{y}=\mathbf{X}_{1}\mathbf{h}+\mathbf{X}_{2}\mathbf{g}+\mathbf{Zs}+\mathbf{e}$

将$\mathbf{h}$所对应的方程吸收到$\mathbf{g}$和$\mathbf{s}$方程中的步骤

## 1.建立忽略$\mathbf{g}$方程与$\sigma_{e}^{2}\mathbf{G}^{-1}$

## 2.将$\mathbf{h}$方程吸收到$\mathbf{s}$方程中

$$
\mathbf{C}\hat{s}=\mathbf{r}\tag{6.21}\\
\text{当}\mathbf{X_{1}^{T}X_{1}}\text{和}\mathbf{Z_{1}Z_{1}}\text{为对角矩阵}
\\
c_{jk} =
\begin{cases}
n_{\cdot k} - \sum\limits_{i} \frac{n_{ik}^2}{n_{i \cdot}} & \text{（当 } j = k \text{ 时，即对角线元素）} \\
-\sum\limits_{i} \frac{(n_{ik} n_{ij})}{n_{i \cdot}} & \text{（当 } j \neq k \text{ 时，即非对角线元素）}
\end{cases}
\\
r_k = y_{\cdot k} - \sum_i \frac{n_{ik} y_{i\cdot}}{n_{i\cdot}}
\\
\mathbf{C} \text{阵和} \mathbf{r} \text{向量中的元素分别具有如下性质：}

\sum_j c_{jk} = 0, \quad \sum_k r_k = 0
$$

> [!NOTE]
>
> ==上述$i$最多取到$n$==

## 3.将$\mathbf{g}$方程加入方程组

* 建立$\mathbf{L}$阵的方法：**对角线方向上是一系列的单位向量，其向量个数等于公牛组数q,每一单位向量的元素个数与对应的公牛所含公牛头数相等**

$$
\mathbf{L} =
\begin{bmatrix}
1_1 &       &        &        \\
    & 1_2   &        &        \\
    &       & \ddots &        \\
    &       &        & 1_q
\end{bmatrix}
\\

\mathbf{D}_2 = \mathbf{CL}, \quad \mathbf{D}_1 = \mathbf{L}'\mathbf{D}_2, \quad \mathbf{q} = \mathbf{L}'\mathbf{r}
\\
\begin{bmatrix}
\mathbf{D}_1 & \mathbf{D}_2' \\
\mathbf{D}_2 & \mathbf{C}
\end{bmatrix}
\begin{bmatrix}
\mathbf{g} \\
\mathbf{s}
\end{bmatrix}
=
\begin{bmatrix}
\mathbf{q} \\
\mathbf{r}
\end{bmatrix}
$$

## 4.将$\sigma_{e}^{2}\mathbf{G}^{-1}$加入到系数矩阵中

$$
\begin{bmatrix}
\mathbf{D}_1 & \mathbf{D}_2' \\
\mathbf{D}_2 & \mathbf{C+\sigma_{e}^{2}G^{-1}}
\end{bmatrix}
\begin{bmatrix}
\mathbf{s} \\
\mathbf{g}
\end{bmatrix}
=
\begin{bmatrix}
\mathbf{r} \\
\mathbf{q}
\end{bmatrix}
$$

# 三、混合模型方程组的迭代求解

## 1.经典迭代方法

设有混合模型方程组：$\mathbf{A} \mathbf{b} = \mathbf{r} $

* 其中 $\mathbf{A}$ 为系数矩阵，$\mathbf{b}$为未知量向量，$\mathbf{r}$为右手向量。

### &#9312;高斯—赛德尔

$$
b_i^{(t)} = b_i^{(t-1)} + \left( r_i - \sum_{j=1}^{i-1} a_{ij} b_j^{(t)} - \sum_{j=i}^{p} a_{ij} b_j^{(t-1)} \right) / a_{ii}
$$

设定初值$b_{i}^{(0)}$进行迭代，至到收敛
$$
\frac{\sum\limits_{i=1}^p \left(b_i^{(t)} - b_i^{(t-1)}\right)^2}{\sum\limits_{i=1}^p \left(b_i^{(t)}\right)^2} < \varepsilon
$$


### &#9313;超松弛法

$$
b_i^{(t)} = b_i^{(t-1)} +\alpha \left( r_i - \sum_{j=1}^{i-1} a_{ij} b_j^{(t)} - \sum_{j=i}^{p} a_{ij} b_j^{(t-1)} \right) / a_{ii}
$$

其中$\alpha$为松弛因子通常取值$\alpha\in（1,2）$，当$\alpha\in(0,1)$为超松弛迭代法

### &#9314;雅可比法

$$
b_i^{(t)} = b_i^{(t-1)} +\alpha \left( r_i  - \sum_{j=i}^{p} a_{ij} b_j^{(t-1)} \right) / a_{ii}
$$

## 2.间接解法
