# 一、等位基因编码

以A为参考等位基因的$i$位点标记效应加性编码

| 基因型 | 101编码 | 012编码 | 中心化编码      |
| ------ | ------- | ------- | --------------- |
| aa     | $-a_i$  | 0       | $-2p_{i}a_i$    |
| Aa     | 0       | $a_i$   | $(1-2p_{i})a_i$ |
| AA     | $-a_i$  | $2a_i$  | $(2-2p_{i})a_i$ |

* 简化为单位点单效应

# 二、先验信息对编辑效应估计的影响

&#9312;在重复实验中，$\sigma^2_a$取较小值可使均值误差最小化

# 三、标记解释的遗传方差

&#9312;单个标记解释的方差为$2pqa_i^{2}$，具有==中等频率的标记将解释大部分遗传变异==，这是忽略忽略等位基因频率标记的原因之一。

## &#9352;标记解释的总遗传方差

$$
\sigma_{a0}^2 = \frac{\sigma_u^2}{2 \sum\limits_i^{nsnp} p_i q_i}\tag{1}
$$

* $2 \sum\limits_i^{nsnp} p_i q_i\hat{a}^2_i$会**低估**总遗传方差。

# 三、RR-SNP与SNP-BLUP

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

* 通过公式1求标记方差，其中$\sigma_u^2$来自早期系谱研究，$p$来自遗传方差估计的群体。
* 多性状：$\mathbf{G_0}$代替性状间遗传协方差矩阵。

### &#9315;标记效应求解

