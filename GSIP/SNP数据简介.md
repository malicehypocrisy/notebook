# 一、SNP数据格式

## &#9352;UGA格式

| 代码 | 基因型   |
| ---- | -------- |
| 0    | AA       |
| 1    | AB or BA |
| 2    | BB       |
| 5    | 缺失     |

## &#9353;标记信息基础检查

### &#9312;检出率

* 观察到的基因型数量
* 典型阈值控制在90%或95%，低于此水平数据将被废弃。

### &#9313;等位基因频率

* 参照等位基因出现频率

### &#9314;次等位基因频率（MAF）

* 两个等位基因频率中的较低值

## &#9353;哈迪—温伯格平衡

| 基因型 | 0      | 1      | 2      |
| ------ | ------ | ------ | ------ |
| 观测值 | $n_0$  | $n_1$  | $n_2$  |
| 预测值 | $nq^2$ | $n2pq$ | $np^2$ |


$$
\chi^2 = \sum \frac{(\text{Observed}_i - \text{Expected}_i)^2}{\text{Expected}_i}
\\

Expected(0:2)= n (q^2, 2pq, p^2)
$$
另一种直接从计数数据获取相同统计量的方法 (Emigh 1980):
$$
\chi^2 = 16n \frac{\left(n_0 n_2 - \frac{n_1^2}{4}\right)^2}{(2n_0 + n_1)^2 (n_1 + 2n_2)^2}
$$
当观测到杂合子数量显著偏离预期值时剔除该标记位点：
$$
\bigg|\frac{n_1}{n}-2 p q\bigg|>t
$$
**t通常取值0.15**

### &#9312;杂交群体的基因频率

若品种A与B的等位基因频率分别为$p_A$和$p_B$，则$Expected(0:2)=n(q_Aq_B,p_Aq_B+q_Ap_B,p_Ap_B)$,用于第一代数据校验

# 二、连锁不平衡

不同位点上等位基因的非随机关联现象

## &#9352;通过配子基因型或个体基因型量化连锁不平衡

### &#9312;D衡量观测分布与预期分布的偏离程度

$$
D=freq(AB)-freq(A)freq(B)
$$

### &#9313;基因含量赋值

* {$A,a$}对应{$0,1$}
* {$B,b$}对应{$0,1$}

用$\mathbf{X}$代表第一等位基因**A**；用$\mathbf{Y}$代表第一等位基因**B**，并计算皮尔逊相关系数。
$$
r=\frac{D}{\sqrt{p_Aq_Ap_Bq_B}}
$$
其中：
$$
p_A=1-q_A=freq(A)
$$

### &#9314;缺少相位计算基因含量

假设存在5个个体，其基因型如下：
$$
\{AB,AB;ab,aB;ab,ab;Ab,AB;Ab,AB\}
$$
转为基因含量：

| X    | 2    | 0    | 0    | 2    | 2    |
| ---- | ---- | ---- | ---- | ---- | ---- |
| Y    | 2    | 1    | 0    | 1    | 1    |

# 三、标记的数量遗传学

**基因含量的取值范围：$z_i\{0,1,2\}$**

* **基因含量的均值**：若以A为参照等位基因时，$E(z)$等于A的出现次数，即$E(z)=2p$
* **基因含量的方差**：在哈迪—温伯格平衡下，$\sigma^2=Var(z)=E(z^2)-E(z)^2=2pq$
* **基因含量遗传力**：1

## &#9352;两个个体加基因含量的协方差

个体$i$和$j$在标记位点拥有两个等位基因拷贝。若从$i$抽取一个拷贝，从$j$抽取另一个拷贝，这两个等位基因“相同”（即同源，identical by descent, IBD）的概率被定义为$\theta$，即亲缘系数。<!--这里还存在疑惑-->
$$
\theta=\frac{A_{ij}}{2}
$$
因此，
$$
\text{Cov}(z_i, z_j) = E(z_i z_j) - E(z_i) E(z_j)=A_{ij}2pq
$$
**推广群体**
$$
E(\mathbf{z})=\mathbf{2}p
\\
Var(\mathbf{z})=\mathbf{A}2pq
\\
\text{其中}\mathbf{A}\text{为经典分子亲缘关系矩阵}
$$

## &#9353;利用基因含量进行质量控制

* ==将基因含量视为数量性状，估计其遗传力==
* 仅需系谱文件与数据文件

$$
\mathbf{z}=\mathbf{1}\mu+\mathbf{Wu}+e
$$

$\mathbf{W}$由基因含量、总体均值、动物编号构成，不存在基因型个体标记为零。

### &#9312;Gengler法估计基础群体确实基因型及等位基因频率

* 为处理缺失基因型的基因含量提供了分析工具
* 可估计未基因分型基础群体的等位基因频率。

## &#9354;基因含量多形状BLUP

$$
\mathbf{y} = \mathbf{X}_y \mathbf{b}_y + \mathbf{W}_y \mathbf{u}_y + \mathbf{e}_y
\\
\mathbf{z} = \mathbf{X}_z \mathbf{b}_z + \mathbf{W}_z \mathbf{u}_z + \mathbf{e}_z

\\
\text{如上所述性状间存在遗传协方差时，}
\mathbf{G}_0 = 
\begin{pmatrix}
\sigma_{u_y}^2 & \sigma_{u_{z,y}} \\
\sigma_{u_{z,y}} & \sigma_{u_z}^2
\end{pmatrix}
 ,\sigma_{u_{z,y}} = 2 p q a
$$

# 四、基因填充

## &#9352;经典填充法

两类信息源：

* 若某个个体染色体可追溯至亲本四条染色体之一，则认为完备。
* 标记位点呈现特定模式，且与已知模式高度匹配

### &#9312;快速近似填充法

**适用范围：**对于快速研究或原型设计，或缺失基因型数量极少的情况(例如个体和 动物检出率≈0.99)可能有用。

* 直接分配最常见基因型(”AA”、”Aa”或其他)。

* 根据概率分布$p^2,2pq,q^2$随机分配基因型

* ```R
  z=sample(c(0,1,2),1,prob=c(p2,2*p*q,q2))
  ```

==不建议使用==

### &#9313;线性填充



# 五、贝叶斯推断

**条件概率：**
$$
p(A,B)=p(A \mid B)p(B)=p(B \mid A)p(A)
$$


**贝叶斯定理表明:**
$$
p(B \mid A)=\frac{p(A \mid B)p(B)}{p(A)}
$$

* $p(A \mid B)$，即似然度。
* $p(B)$ ，即"先验"概率。
*  $p(A)$可以通过全概率公式计算得到，称为边缘概率（==A发生的所有的概率之和==）

### &#9312;贝叶斯公式推导

$$
\begin{aligned}
p(B \mid A) &= \frac{p(AB)}{p(B)} \\
&= \frac{p(A \mid B)p(B)}{p(A \mid B)p(B) + p(A \mid \overline{B})p(\overline{B})}
\end{aligned}
$$

### &#9313;一元正态分布密度函数

对于随机变量$X \sim N(\mu, \sigma^2)$，其概率密度函数为：
$$
f(x)=\frac{1}{\sigma \sqrt{2\pi}} \exp\left(-\frac{(x-\mu)^2}{2\sigma^2}\right)
$$

* 中心极限定理

### &#9314;多元正态分布

$$
\text{对于} d \text{维随机向量} \mathbf{x} = (x_1, x_2, \ldots, x_d)^T \text{，服从多元正态分布} \mathbf{x} \sim \mathcal{N}_d(\boldsymbol{\mu}, \boldsymbol{\Sigma}) \text{的概率密度函数为：}
\\
f(\mathbf{x}) = \frac{1}{(2\pi)^{d/2} |\boldsymbol{\Sigma}|^{1/2}} \exp\left( -\frac{1}{2} (\mathbf{x} - \boldsymbol{\mu})^T \boldsymbol{\Sigma}^{-1} (\mathbf{x} - \boldsymbol{\mu}) \right)\tag{1}
$$

## &#9352;贝叶斯推断示例

&#9312;若$\mu$值取值有限可以使用公式1。

&#9313;若均值分布是连续的
$$
p(\mu \mid \mathbf{y}) = N\left( \hat{\mu}, l h s^{-1} \right)
\\

\text{其中}

\begin{aligned}
lhs &= \frac{\mathbf{1}' \mathbf{1}}{\sigma_e^2} + \frac{1}{\sigma_{\mu}^2} \\
\hat{\mu} &= \left( lhs^{-1} \right) \mathbf{1}' \mathbf{y} \sigma_e^2
\end{aligned}
$$

## &#9353;吉布斯采样器

## &#9354;吉布斯采样后分析
