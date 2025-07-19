# 一、等位基因编码

以A为参考等位基因的$i$位点标记效应加性编码

| 基因型 | 101编码 | 012编码 | 中心化编码      |
| ------ | ------- | ------- | --------------- |
| aa     | $-a_i$  | 0       | $-2p_{i}a_i$    |
| Aa     | 0       | $a_i$   | $(1-2p_{i})a_i$ |
| AA     | $-a_i$  | $2a_i$  | $(2-2p_{i})a_i$ |

* 简化为单位点单效应

# 二、先验信息对编辑效应估计的影响

* 在重复实验中，$\sigma^2_a$取较小值可使均值误差最小化

# 三、标记解释的遗传方差

单个标记解释的方差为$2pqa_i^{2}$，具有==中等频率的标记将解释大部分遗传变异==，这是忽略等位基因频率标记的原因之一。

## &#9352;标记解释的总遗传方差

$$
\sigma_{a0}^2 = \frac{\sigma_u^2}{2 \sum\limits_i^{nsnp} p_i q_i}\tag{1}
$$

* $2 \sum\limits_i^{nsnp} p_i q_i\hat{a}^2_i$会**低估**总遗传方差。

# 四、RR-SNP与SNP-BLUP

SNP-BLUP模型通常表示为：
$$
\mathbf{u}=\mathbf{Xb}+\mathbf{Za}+\mathbf{e}
$$

* $\mathbf{a}$为标记效应

**假设标记见互相独立**
$$
p(\mathbf{a}) = MVN(0, \mathbf{D}) \quad ; \quad \text{Var}(\mathbf{a}) = \mathbf{D} = \mathbf{I} \sigma_{a0}^2
$$

## &#9352;MME

$$
\begin{pmatrix}
\mathbf{X}'\mathbf{R}^{-1}\mathbf{X} & \mathbf{X}'\mathbf{R}^{-1}\mathbf{Z} \\
\mathbf{Z}'\mathbf{R}^{-1}\mathbf{X} & \mathbf{Z}'\mathbf{R}^{-1}\mathbf{Z} + \mathbf{D}^{-1}
\end{pmatrix}
\begin{pmatrix}
\widehat{\mathbf{b}} \\
\widehat{\mathbf{a}}
\end{pmatrix}
=
\begin{pmatrix}
\mathbf{X}'\mathbf{R}^{-1}\mathbf{y} \\
\mathbf{Z}'\mathbf{R}^{-1}\mathbf{y}
\end{pmatrix}
$$

### &#9312;单性状模型

若 $\text{Var}(\mathbf{a}) = \mathbf{D} = \mathbf{I} \sigma_{a0}^2;\text{Var}(\mathbf{e}) = \mathbf{R} = \mathbf{I} \sigma_e^2，$则可简化为：
$$
\begin{pmatrix}
\mathbf{X}'\mathbf{X} & \mathbf{X}'\mathbf{Z} \\
\mathbf{Z}'\mathbf{X} & \mathbf{Z}'\mathbf{Z} + \mathbf{I} \lambda
\end{pmatrix}
\begin{pmatrix}
\widehat{\mathbf{b}} \\
\widehat{\mathbf{a}}
\end{pmatrix}
=
\begin{pmatrix}
\mathbf{X}'\mathbf{y} \\
\mathbf{Z}'\mathbf{y}
\end{pmatrix}
$$
其中：$ \lambda=\frac{\sigma_{e}^2}{\sigma_{ao}^2}$

可简记为：
$$
lhs\begin{pmatrix}
\mathbf{\hat{b}} \\
\mathbf{\hat{a}}
\end{pmatrix}=rhs
$$
特点：

* 维度=$(固定效应数+标记数)^2$,与个体数量无关
* $\mathbf{Z^{T}Z}$完全稠密且非稀疏

### &#9313;多性状

$$
\mathbf{R} = \mathbf{I} \otimes \mathbf{R}_0 \quad \text{和} \quad \mathbf{D} = \mathbf{I} \otimes \mathbf{S}_{a0}
$$

### &#9314;BLUP-SNP中的方差组分

MME假设已知两个方差组分的值：$\sigma_a^2$与$\sigma_e^2$。

* 通过`[公式1]`求标记方差，其中$\sigma_u^2$来自早期系谱研究，$p$来自遗传方差估计的群体。
* 多性状：$\mathbf{G_0}$代替性状间遗传协方差矩阵。

### &#9315;标记效应求解

#### &#9332; **GSRU（高斯-塞德尔残差更新法）**

$$
\begin{aligned}
(\mathbf{z}_i' \mathbf{z}_i + \lambda) \widehat{a}_i^{l+1} &= \mathbf{z}_i' \bigg( \mathbf{y} - \mathbf{X} \widehat{\mathbf{b}} - \mathbf{Z} \widehat{\mathbf{a}} + \mathbf{z}_i \widehat{a}_i^l \bigg ) \\
\text{故方程组为：} \quad (\mathbf{z}_i' \mathbf{z}_i + \lambda) \widehat{a}_i^{l+1} &= \mathbf{z}_i' \big( \widehat{\mathbf{e}}^l + \mathbf{z}_i \widehat{a}_i^l \big) = \mathbf{z}_i' \widehat{\mathbf{e}}^l + \mathbf{z}_i' \mathbf{z}_i \widehat{a}_i^l \\
\text{每次标记效应更新后需用下式修正误差项：} \quad \widehat{\mathbf{e}}^{l+1} &= \widehat{\mathbf{e}}^l - \mathbf{z}_i \big( \widehat{a}_i^{l+1} - \widehat{a}_i^l \big)
\end{aligned}
$$

* 只需极少改动即可转换为吉布斯采样器

#### &#9333;**PCG（预处理共轭梯度法）**

特点：**采用通用求解器并通过连续计算向量积**
$$
\left( \begin{matrix} \mathbf{X'X} & \mathbf{X'Z} \\ \mathbf{Z'X} & \mathbf{Z'Z} + \mathbf{I}\lambda \end{matrix} \right) \left( \begin{matrix} \widehat{\mathbf{b}}^l \\ \widehat{\mathbf{a}}^l \end{matrix} \right)  =  \left( \begin{matrix} \mathbf{X'} \\ \mathbf{Z'} \end{matrix} \right) \left( \left( \begin{matrix} \mathbf{X} & \mathbf{Z} \end{matrix} \right) \left( \begin{matrix} \widehat{\mathbf{b}}^l \\ \widehat{\mathbf{a}}^l \end{matrix} \right) \right)  +  \left( \begin{matrix} \mathbf{0} \\ \mathbf{I}\lambda \widehat{\mathbf{a}}^l \end{matrix} \right)
$$

## &#9353;基于标记效应：BayesC($P_i=0$)

* 通过贝叶斯推断可简便实现，并获得方差$\sigma_{a}^2$和$\sigma_{e}^2$的后验估计
* 该算法需要方差的初始值及其先验信息。

### &#9312;先验分布公式

$$
\text{方差} \sim \frac{S^2_a \cdot \nu}{\chi_\nu^2}
$$

- $S^2$：尺度参数（决定先验期望）
- $\nu$：自由度（决定先验强度），==通常选择 `4`==
- $\chi_\nu^2$：卡方分布抽样

### &#9313;先验参数设定

$$
S_{e}^{2}=\sigma_{e}^{2} \nu_{e};E(\sigma^2)=\frac{S^{2}}{\nu}
\\
S_{a}=\sigma_{a 0}^{2} \nu_{a} ; \sigma_{a 0}^{2}=\frac{\sigma_{u}^{2}}{2\sum_{i}^{\text{nsnp}} p_{i} q_{i}}
$$

## &#9354;将标记方差转化为遗传方差

注意：相关论文证实使用`公式1`的误差不大于10%

若某群体同时具备系谱记录与分子标记数据，在哈迪-温伯格平衡且无近交条件下，两种方差（**系谱与标记**）理论预期值应相同。

## &#9355;标记物的差异方差

**强效应 QTLs（主效基因）**：被过度收缩（估计值 $\widehat{a}_i$ 远小于真实值）

# 五、BayesA

因为 **t 分布**的方差是 $\nu / (\nu - 2)$，这相当于对标记效应采用缩放 t 分布作为先验:

$$
p\left(a_{i} \mid \sigma_{a 0}^{2}, \nu_{a}\right)=\sigma_{a 0} t\left(0, \nu_{a}\right)
$$
特点：**该分布具有”厚尾”特性，意味着大标记效应在先验中并非小概率事件。**

## &#9352;**标记效应条件分布**（正态先验）

$$
p\left(a_{i} \mid \sigma_{a i}^{2}\right)=N\left(0, \sigma_{a i}^{2}\right)
$$

## &#9353;**方差先验分布**（逆卡方分布）

$$
p\left(\sigma_{a i}^{2} \mid S_{a}, \nu_{a}\right)=S_{a} \chi_{\nu_{a}}^{-2}
$$

* 若 $(X \sim \chi_\nu^2)$，则 $\sigma^2 = S^2 / X$服从逆卡方分布
* 期望：$E(\sigma_{ai}^2) = \frac{S_a^2}{\nu_a - 2}$（要求 $\nu_a > 2$）

## &#9354;先验参数设定

$$
\sigma_{a 0}^{2}=\frac{\nu-2}{\nu} \frac{\sigma_{u}^{2}}{2 \sum_{i} p_{i} q_{i}}
$$

# 六、BayesB

产生原因：**多数标记与 QTL 无关，其真实效应为零**（$\sigma^2_{ai}$）。

**注意：**

* ==这在BayesA中不可能出现==
* <u>该方法对$S^2_a,ν_a$和$π$的 先验值较为敏感。</u>

$$
p\left(a_{i} \mid \sigma_{a i}^{2}\right)=N\left(0, \sigma_{a i}^{2}\right)
\\
\left\{\begin{aligned}
p\left(\sigma_{a i}^{2} \mid S_{a}, \nu_{a}\right) &= S_{a}\chi_{\nu_{a}}^{-2} &\text{with probability } 1-\pi \\
p\left(\sigma_{a i}^{2} \mid S_{a}, \nu_{a}\right) &= 0 &\text{with probability } \pi
\end{aligned}\right.
$$

# 七、BayesC

特点：

* 将有效应的标 记赋予”共同”方差$σ^2_{a0}$。
* 最终(数据拟合后)$σ^2_{a0}$的取值几乎不依赖于先验设定，因此即**使先验设定有误，模型仍可能正确。**

$$
p\left(a_{i} \mid \delta_{i}\right)=\begin{cases} 
N\left(0, \sigma_{ai}^{2}\right) & \text{if } \delta_{i}=1 \\
0 & \text{otherwise}
\end{cases}
\\
p\left(\sigma_{a 0}^{2} \mid S_{a}, \nu_{a}\right)=S_{a} \chi_{\nu_{a}}^{-2}
\\
p\left(\delta_{i}=1\right)=1-\pi
$$



其中 $S_{a}$ 可设置为类似 $S_{a}^{2}=\sigma_{a 0}^{2} \nu_{a 0}$ 的形式

$$
\sigma_{a 0}^{2}=\frac{\sigma_{u}^{2}}{(1-\pi) 2 \sum p_{i} q_{i}}
$$

## &#9352;与性状关联的标记

使用 **Beta分布** 描述 $π$ 的不确定性。

注意：BayesC的输出结果是$\delta_{i}$，即后验均值。该值不会是非0 即1，而是介于两者之间，故**不能用于筛选”控制**。

### &#9312;判断显著性：贝叶斯因子

$$
BF=\frac{\frac{p(\text{SNP in the model} \mid data)}{p(\text{SNP not in the model} \mid data)}}{\frac{p(\text{SNP in the model})}{p(\text{SNP not in the model})}}
$$

本例中计算公式为：

$$
BF_{i}=\frac{(1-\pi)}{\pi} \frac{p\left(\delta_{i}=1 \mid \mathbf{y}\right)}{1 - p\left(\delta_{i}=1 \mid \mathbf{y}\right)}
$$

*  BF=3-20 为” 提示性”证据 

* BF=20-150 为” 强”证据 

* BF>150 为”极强”证据

* 注意：由于所有单核苷酸多态性(SNP)都是同时引入的，且先验分布已对其估计值进行了惩罚，因此==无需进行多重检验校正==。

# 八、贝叶斯锁套法

特点：

* **将所有标记效应设为接近零的极小值**，通过连续压缩实现稳健建模。
* 用 **双指数先验（拉普拉斯分布）**，对效应施加**L1正则化**（绝对收缩），迫使大部分标记效应趋近零但不严格等于零。

## &#9352;先验分布公式

$$
p(a_i|\lambda)=\frac{\lambda}{2 \sigma} \exp \left( -\frac{\lambda |a_i|}{\sigma} \right)
\\
p\left(a_{i}\mid \sigma_{ai}^{2}\right)=N\left(0, \sigma_{ai}^{2}\right)
\\
p\left(\sigma_{ai}^{2}\mid \lambda\right)=\frac{\lambda^{2}}{2}\exp\left(-\frac{\lambda^{2}}{2}\frac{\sigma_{ai}^{2}}{\sigma^{2}}\right)
$$

*  此处的$\lambda$与BLUP-SNP中的$\lambda$无关

## &#9353;参数化

### &#9312;设定$σ^2$=1

$$
\text{此时}Var(a_i|\lambda)=\frac{2}{λ^2}
\\
\frac{2}{\lambda^{2}}=\frac{\sigma_{u}^{2}}{2 \sum p_{i} q_{i}}
$$

### &#9313;设定$\sigma^2=\sigma^2_e$

$$
\frac{2}{\lambda^{2}}=\frac{\sigma_{u}^{2}}{\sigma_{e}^{2} 2 \sum p_{i} q_{i}}
$$

# 九、VanRaden 非线性方法

$$

\sigma_{a i}^{2}=\sigma_{a 0}^{2}\left(c^\left(\frac{|\widehat{a}_{i}|}{\text{sd}\left(\widehat{a}_{1},\dots,\widehat{a}_{n}\right)}-2\right)\right)
$$

## &#9352;**“曲率”参数 c 的取值范围**

$$
c\in[1,1.25]
$$

## &#9353;**数值稳定性措施**

- **标准化限制**：对 $\frac{|\widehat{a}_{i}|}{\text{sd}()}$ 设置上限为5（避免指数爆炸）。

  **小数据集推荐**：使用 $c=1.12$ （平衡效率与稳定性）。

# 十、标记模型的可靠性

## &#9352; 通过 MCMC贝叶斯方法获得的标准误

### &#9312;MCMC迭代抽样

在贝叶斯框架下，MCMC（马尔科夫链蒙特卡洛）通过选代生成标记效应分布的样本：

- 第$t$次迭代：获得标记效应向量样本 $\tilde{a}_{(t)}$
- 转化为育种值：
  $$
  \tilde{u}^{(t)} = Z \tilde{a}^{(t)}
  \\
  \text{其中}Z\text{是基因型矩阵（行=个体，列=标记）}
  $$

### &#9313;后验估计计算

MCMC链收敛后（通常舍前 $n$ 次预热迭代）：

- 育种值估计 \(\tilde{u}
  $$
  _i\)
  = 个体 \(i\) 所有选代育种值的后验均值：
  \[
  \tilde{u}_i = \frac{1}{T} \sum_{t=1}^{T} \tilde{u}_{i}^{(t)}
  $$
  \]
  （\(T\) 为有效选代次数）
  
- 后验方差 \(\text{Var}(\tilde{u}_i)\)
  = 个体 \(i\) 育种值样本的方差：
  \[
  $$
  \text{Var}(\tilde{u}_i) = \frac{1}{T-1} \sum_{t=1}^{T} \left( \tilde{u}_{i}^{(t)} - \tilde{u}_i \right)^2
  $$

### &#9314;标准误差与置信区间

- 标准误差 (SE)：
  $$
  \text{SE}(\tilde{u}_i) = \sqrt{\text{Var}(\tilde{u}_i)}
  意义：衡量育种值估计的不确定性
  $$

- 95%置信区间：
  $$
  \tilde{u}_i \pm 2 \times \text{SE}(\tilde{u}_i)
  $$

## &#9353;可靠性

注意：但$z_iz^′_i$依赖参数化选 择——即便获得完全相同的育种值，==不同编码方式会导致可靠性估值差异==。

### &#9312;**方法1：基于MCMC的估计方差**

$$
Rel_i = 1 - \frac{\text{Var}(\widehat{u}_i)}{z_i z_i' \sigma_{a0}^2}
$$



- $Var(\widehat{u}_i)$：从MCMC（马尔科夫链蒙特卡洛）抽样获得的育种值方差（后验方差）。

- $z_i z_i' \sigma_{a0}^2$：个体i的基因组方差（遗传方差）


### &#9313;方法2：基于后验协方差矩阵 ($C^{aa}$\)

#### &#9332;获取标记效应后验协方差矩阵：


$$
\text{Var}(a|y)=C^{aa}
$$

- 通过MCMC或SNP-BLUP方程矩阵求逆计算（需解析方法）。

$$
\text{Var}(u|y)=Z C^{aa} Z'
$$

#### &#9333;可靠性公式：


$$
\text{Rel}_i = 1 - \frac{z_i C^{aa} z_i'}{z_i z_i' \sigma_{a0}^2}
$$
