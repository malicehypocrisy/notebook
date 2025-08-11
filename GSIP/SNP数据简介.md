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
* **典型阈值控制在90%或95%**，低于此水平数据将被废弃。

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
另一种直接从计数数据获取相同统计量的方法 :
$$
\chi^2 = 16n \frac{\left(n_0 n_2 - \frac{n_1^2}{4}\right)^2}{(2n_0 + n_1)^2 (n_1 + 2n_2)^2}
$$
当观测到杂合子数量**显著偏离预期值**时剔除该标记位点：
$$
\bigg|\frac{n_1}{n}-2 p q\bigg|>t
$$
* t通常取值**$0.15$**

### &#9312;杂交群体的基因频率

若品种A与B的等位基因频率分别为$p_A$和$p_B$，则$Expected(0:2)=n(q_Aq_B,p_Aq_B+q_Ap_B,p_Ap_B)$==，用于第一代数据校验==

## &#9355;连锁不平衡

* 指群体中不同位点的等位基因组合频率偏离随机预期，随机关联的现象。

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

# 二、标记的数量遗传学

**基因含量的取值范围：$z_i\{0,1,2\}$**

* **基因含量的均值**：若以A为参照等位基因时，$E(z)$等于A的出现次数，即$E(z)=2p$
* **基因含量的方差**：在哈迪—温伯格平衡下，$\sigma^2=Var(z)=E(z^2)-E(z)^2=2pq$
* $\text{基因含量遗传力}=1$

## &#9352;基因含量的协方差

**亲缘系数的定义**：个体$i$和$j$在标记位点拥有两个等位基因拷贝。若从$i$抽取一个拷贝，从$j$抽取另一个拷贝，这两个等位基因“相同”（即同源）的概率被定义为$\theta$，即。
$$
\theta=\frac{A_{ij}}{2}
\\
\text{Cov}(z_i, z_j) = E(z_i z_j) - E(z_i) E(z_j)=A_{ij}2pq
$$


### &#9312;**推广到群体**

$$
E(\mathbf{z})=\mathbf{2}p 
\\
Var(\mathbf{z})=\mathbf{A}2pq
$$

* $\mathbf{A}$为经典分子亲缘关系矩阵

## &#9353;利用基因含量遗传力进行质量控制

* ==将基因含量视为数量性状，估计其遗传力==
* 仅需系谱文件与数据文件

$$
\mathbf{z}=\mathbf{1}\mu+\mathbf{Wu}+e
$$

* $\mathbf{W}$由基因含量、总体均值、动物编号构成，不存在基因型个体标记为零。
* 若$\hat{h}^2<0.98$时，则基因型或系谱数据存在问题

## &#9354;Gengler法

### &#9312;作用

* 为处理缺失基因型的基因含量提供了分析工具
* 可估计未基因分型基础群体的等位基因频率。

### &#9313;线性模型

$$
\mathbf{z}=1\mu+\mathbf{Wu}+\mathbf{e}
$$

* $\mathbf{W}$为关联矩阵，基因型个体标记为1，其余为0
* 使用0.99的遗传力
* $\mathbf{\hat{u}}$包含所有个体的基因含量估计值
* $\hat{\mu}$实际包含$2\hat{p}$

### &#9314;不足之出

* 主效基因上的相似度高于实际值

## &#9355;基因含量多性状BLUP

$$
\mathbf{y} = \mathbf{X}_y \mathbf{b}_y + \mathbf{W}_y \mathbf{u}_y + \mathbf{e}_y
\\
\mathbf{z} = \mathbf{X}_z \mathbf{b}_z + \mathbf{W}_z \mathbf{u}_z + \mathbf{e}_z
$$

* $y$常规性状
* $z$基因含量
* 性状间协方差：

$$
\mathbf{G}_0 = 
\begin{pmatrix}
\sigma_{u_y}^2 & \sigma_{u_{z,y}} \\
\sigma_{u_{z,y}} & \sigma_{u_z}^2
\end{pmatrix}
\\
 \sigma_{u_{z,y}} = 2 p q a
$$

# 三、基因填充

$$
\text{经典填充法}=
\begin{cases}
\text{快速近似填充法}\rightarrow{根据概率分布}p^2,2pq,q^2\text{随机分配基因型}
\\
\text{线性填充}\rightarrow ssGBLUP
\end{cases}
$$

**填充依据**：

* 若某个个体染色体可追溯至亲本四条染色体之一，则认为完备。
* 标记位点呈现特定模式，且与已知模式高度匹配

# 四、贝叶斯推断

**条件概率：**
$$
p(A,B)=p(A \mid B)p(B)=p(B \mid A)p(A)
$$

**贝叶斯定理解读:**
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

$$
\mathbf{y}=\mathbf{1}\mu+\mathbf{e}
\\
\text{其中，}Var(\mathbf{e})=\mathbf{R}=\mathbf{I\sigma^2_e}\text{和}\sigma^2_e=1
$$

&#9312;若$\mu$值取值**有限**，使用`公式 1`。

&#9313;若均值分布是**连续的**，则
$$
p(\mu \mid \mathbf{y}) = N\left( \hat{\mu}, l h s^{-1} \right)
\\
\begin{aligned}
lhs &= \frac{\mathbf{1}' \mathbf{1}}{\sigma_e^2} + \frac{1}{\sigma_{\mu}^2} \\
\hat{\mu} &= \left( lhs^{-1} \right) \mathbf{1}' \mathbf{y} /\sigma_e^2
\end{aligned}
\\
\text{其中，}\mathbf{1}\text{：是一个}n\times1\text{的列向量}
$$

## &#9353;吉布斯采样器

* 假设初始值$\mu = 0$和$\sigma_e^2=1$，从一下分布抽取新样本$\mu$

$$
p(\mu|y,\sigma_e^2) = N(\hat{\mu},lhs^{-1})
$$

* 从以下分布抽取$\sigma_e^2$

$$
p(\sigma_e^2|\mathbf{y},\mu)=(\mathbf{y-1}\mu)^`(\mathbf{y-1}\mu)\chi^{-2}_k
$$

## &#9354;吉布斯采样后分析

### &#9312;**预处理**

* 丢弃前n次迭代，每5次采样保留1个样本。

### &#9313;**收敛诊断（确保采样有效）**

* **R-hat统计量**：通过比较多条链的方差判断收敛，值越接近1（<1.01)越好。若R-hat过大，则需要增加迭代次数。
* **自相关函数(ACF)**：MCMC样本通常存在自相关，ACF应快速衰减到0，若自相关强（如滞后10仍>0.5)，可通过每隔n个样本保留1个来减少相关性。
* **有效样本数（ESS）**：衡量样本中独立信息的数量，ESS = 总样本量、（1+2x自然相和）。ESS越大，厚颜越可靠（一般要求>1000).

### &#9314;后验统计量与不确定性量化

* **点估计**：常用后验均值（兼顾信息）或中位数（对抗极端值）作为参数估计。
* **可信区间**：分位数法（非参数，适用于任何分布）比正态近似更稳健。

### &#9315;可视化分析

* **密度图**：展示后验分布形状
* **分段密度对比**：将后验样本分成多端，若密度曲线重合，说明链已平稳。
* **累计分布函数（CDF）**：可只管判断参数落在某个区间的概率

### &#9316;模型验证

* 从后验预测分布生成”模拟数据“与真实数据进行比较分析



