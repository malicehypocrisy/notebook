> [!NOTE]
>
> 解决的问题：**特定群体中仅有少数个体被基因分型。**
>
> 解决的方案：**将系谱关系与基因组关系合并，并将其作为MME的协方差结构。**

# 一、作为改进关系的SSGBLUP

$$
\begin{align}
\mathbf{H} & = \begin{pmatrix}
\operatorname{var}\left(\mathbf{u}_1\right) & \operatorname{cov}\left(\mathbf{u}_1, \mathbf{u}_2\right) \\
\operatorname{cov}\left(\mathbf{u}_2, \mathbf{u}_1\right) & \operatorname{var}\left(\mathbf{u}_2\right)
\end{pmatrix}
\\
& = \begin{pmatrix}
\mathbf{A}_{11} + \mathbf{A}_{12} \mathbf{A}_{22}^{-1} (\mathbf{G} - \mathbf{A}_{22}) \mathbf{A}_{22}^{-1} \mathbf{A}_{21} & \mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{G} \\
\mathbf{G}\mathbf{A}_{22}^{-1}\mathbf{A}_{21} & \mathbf{G}
\end{pmatrix}
\\
&=\mathbf{A} + 
\begin{bmatrix}
\mathbf{A}_{12}\mathbf{A}_{22}^{-1} & \mathbf{0} \\
\mathbf{0} & \mathbf{I}
\end{bmatrix}
\begin{bmatrix}
\mathbf{I} \\
\mathbf{I}
\end{bmatrix}
[\mathbf{G} - \mathbf{A}_{22}]
\begin{bmatrix}
\mathbf{I} & \mathbf{I}
\end{bmatrix}
\begin{bmatrix}
\mathbf{A}_{22}^{-1}\mathbf{A}_{21} & \mathbf{0} \\
\mathbf{0} & \mathbf{I}
\end{bmatrix}
\end{align}
$$

* 使用$\mathbf{A}$建模亲缘关系时，$\mathbf{G}$应具有“加性”特征

$$
\mathbf{H}^{-1}=\mathbf{A}^{-1}+\begin{bmatrix}
\mathbf{0} & \mathbf{0} \\
\mathbf{0} & \mathbf{G^{-1}-A_{22}^{-1}}
\end{bmatrix}
$$

# 二、作为线性插值的SSGBLUP

## &#9352;将基因型数据$\mathbf{z}$建模为**数量性状**

$$
\mathbf{z}=\mathbf{1}\mu+\mathbf{Wu}+\mathbf{e}
$$

* 其中$\mu$是基础群体等位基因频率$(\mu=2p)$，$\mathbf{u}$是随机遗传效应
* $\mathbf{\hat{z}}_1=\mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{z}_2$

## &#9353;**基因组关系矩阵构建**

$$
\hat{\mathbf{G}}=\frac{\hat{\mathbf{Z}}_1\hat{\mathbf{Z}}_1^`}{\sum{2p_jq_j}}
$$

**插补误差问题**
$$
\text{条件期望：}E(\mathbf{Z_1|Z_2})=\mathbf{\hat{Z}}_1
\\
\text{条件方差:}\operatorname{Var}(\widehat{\mathbf{Z}}_1 \mid \mathbf{Z}_2) = (\mathbf{A}_{11} - \mathbf{A}_{12} \mathbf{A}_{22}^{-1} \mathbf{A}_{21}) \mathbf{V}
\\
\text{其中，}\mathbf{V}\text{的对角线包含}p_kq_k
$$


## &#9354; **SSGBLUP 最终方差结构**

$$
\operatorname{Var}(\mathbf{u}_1) = \sigma_u^2 \left( \frac{\widehat{\mathbf{Z}}_1 \widehat{\mathbf{Z}}_1'}{2 \sum p_k q_k} + \mathbf{A}_{11} - \mathbf{A}_{12} \mathbf{A}_{22}^{-1} \mathbf{A}_{21} \right)
$$



# 三、混合方程

$$
\mathbf{y} = \mathbf{Xb} + \mathbf{Wu} + \mathbf{e}
\\
Var(\mathbf{u})=\mathbf{H}\sigma_{u}^2;Var(\mathbf{e})=\mathbf{I}\sigma_e^2
$$

## &#9352;对于单一性状

$$
\left( \begin{array}{cc}
\mathbf{X}'\mathbf{X}\sigma_e^{-2} & \mathbf{X}'\mathbf{W}\sigma_e^{-2} \\
\mathbf{W}'\mathbf{X}\sigma_e^{-2} & \mathbf{W}'\mathbf{W}\sigma_e^{-2} + \mathbf{H}^{-1}\sigma_u^{-2}
\end{array} \right)
\left( \begin{array}{c}
\widehat{\mathbf{b}} \\
\widehat{\mathbf{u}}
\end{array} \right)
=
\left( \begin{array}{c}
\mathbf{X}'\mathbf{y}\sigma_e^{-2} \\
\mathbf{W}'\mathbf{y}\sigma_e^{-2}
\end{array} \right)
$$
## &#9353;对于多性状情况

$$
\left( \begin{array}{cc}
\mathbf{X}'\mathbf{R}^{-1}\mathbf{X} & \mathbf{X}'\mathbf{R}^{-1}\mathbf{W} \\
\mathbf{W}'\mathbf{R}^{-1}\mathbf{X} & \mathbf{W}'\mathbf{R}^{-1}\mathbf{W} + \mathbf{H}^{-1} \otimes \mathbf{G}_0
\end{array} \right)
\left( \begin{array}{c}
\widehat{\mathbf{b}} \\
\widehat{\mathbf{u}}
\end{array} \right)
=
\left( \begin{array}{c}
\mathbf{X}'\mathbf{R}^{-1}\mathbf{y} \\
\mathbf{W}'\mathbf{R}^{-1}\mathbf{y}
\end{array} \right)
$$

* $\mathbf{G}_0$是性状间遗传协方差矩阵，通常$\mathbf{R}=\mathbf{I}\otimes\mathbf{R}_{0}$，$\mathbf{R}_0$为协方差矩阵
* $\mathbf{H}$中的近交系数:

$$
F_i=H_{ii}-1
$$

# 三、$\mathbf{G}$和$\mathbf{A}$矩阵的融合

**兼容性**：让$\mathbf{G}$和$\mathbf{A}$具有相同尺度

**混合**：通过技术手段将部分遗传方差分配给系谱（非标记）,确保$\mathbf{G}$可逆

## &#9352;混合处理

在单步法中，“标记基础”部分对应亲缘矩阵$\mathbf{H}$，而“系谱基础“部分对应亲缘矩阵$\mathbf{A}$。
$$
Var(\mathbf{u})=((1-\alpha)\mathbf{H}+\alpha\mathbf{A})\sigma_u^2\text{其中，}\alpha=\sigma^2_{u,m}/\sigma^2_{u}
$$

* **在实践中，$\alpha$取值较低且对预测结果影响甚微**。

**使用REML方法**
$$
\text{拟合两个独立随机效应：}
\\
\mathbf{u}_m∼N(0,\mathbf{H}σ_{u,m}^2)

\\
\mathbf{u}_p∼N(0,\mathbf{A}σ_{u,p}^2)
$$
**优势：**

* 精确估计$\mathbf{H}σ_{u,m}^2$和$\mathbf{A}σ_{u,p}^2$而非固定

### &#9312;使用技巧实现可逆性

$$
\mathbf{G}_w\leftarrow\mathbf{G}+0.01\mathbf{I}
$$

* 因此，$\mathbf{G}_w\approx\frac{\mathbf{ZZ^`}}{\sum{(2p_jq_j)}}$

$$
\mathbf{G}_w\leftarrow(1-\alpha)\mathbf{G}+\alpha\mathbf{A}_{22}
$$

* 与$\mathbf{G}$的相似度降低，对SNP效应的反向求解产生影响

## &#9353;兼容性

### &#9312;将$\mathbf{G}$矩阵适配至$\mathbf{A}$矩阵

要求：$\mathbf{G}$必须**基于基础等位基因频率构建**。

然而常采用观测群体的等位基因频率构建$\mathbf{G}$，合理的做法是为群体设置均值：
$$
\mathbf{u}_2\sim N(\mathbf{1}\mu,\mathbf{G}\sigma_u^2)
\\
\text{将}\mu\text{作为随机效应：}\mathbf{G^*}=\mathbf{G}+\mathbf{11^`}a\rightarrow a=\mathbf{\overline{A}_{22}}-\mathbf{\overline{G}}
$$
故：使用$a$可得到均值为0的$\mathbf{u}_2\sim N(0,\mathbf{G^*}\sigma_u^2)$

### &#9313;拟合$\mathbf{A}$到$\mathbf{G}$

目的：使$\mathbf{A}_{22}$趋近于$\mathbf{G}$

问题提出：任何等位基因频率都存在不确定性
$$
\mathbf{A}_{22}=\gamma+\mathbf{I}(1-\gamma/2)
\\
\mathbf{G}=\mathbf{ZZ^`}/s
\\
s=\frac{m}{2}\\
\gamma = 8var(p_{base})
$$
步骤：

* 使用Gengler法估算$p_{base}$。
* 计算$\hat{\gamma}=8Var(\hat{p}_{base})$

### &#9314;未知亲本群组

伪亲本并非真实个体，而是被构想为无线提供后代的基因池；他们后代既无近交也无亲缘关系。

UPG的应用通过特殊矩阵$\mathbf{A^*}$实现，该矩阵在BlUP中承担$\mathbf{A^{-1}}$的功能：
$$
\mathbf{A}^* = \begin{pmatrix}
\mathbf{Q}'\mathbf{A}^{-1}\mathbf{Q} & -\mathbf{Q}'\mathbf{A}^{-1} \\
-\mathbf{A}^{-1}\mathbf{Q} & \mathbf{A}^{-1}
\end{pmatrix}
\\
p(\mathbf{u}_2)=N(\mathbf{Qg},\mathbf{A}\sigma_u^2)
\\
\text{其中，}\mathbf{g}\text{未知亲本群的在“育种值”，}\mathbf{Q}\text{包含血统成分}
$$

#### &#9332;截断系谱与数据

* 删除”陈旧“数据后逆向追溯三代系谱，此时无需使用未知亲本
* 由于(截断后)系谱的基准世代与基因型个体高度接近，$\mathbf{A}$与$\mathbf{G}$几乎自动匹配。

#### &#9333;近似处理

$$
\mathbf{H^*}=\mathbf{A^*}+\begin{bmatrix}
\mathbf{0} & \mathbf{0} \\
\mathbf{0} & \mathbf{G^{-1}-A_{22}^{-1}}
\end{bmatrix}
\\
\text{其中，}\mathbf{A}_{22}\text{忽略位置亲本的方式构建}
$$

#### &#9334;元组群体

* 假设所有等位基因频率均为$0.5$。
* 使用元组群体的伪动物替代UPGs，这些群体名为$\Gamma$的特殊亲缘系数。

$$
\mathbf{H}^{\Gamma^{-1}} = \mathbf{A}^{\Gamma^{-1}} + 
\begin{pmatrix}
\mathbf{0} & \mathbf{0} \\
\mathbf{0} & \mathbf{G}^{-1} - \mathbf{A}_{22}^{\Gamma^{-1}}
\end{pmatrix}
$$

