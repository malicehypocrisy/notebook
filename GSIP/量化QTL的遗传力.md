# 一，QTL的基础

QTL遗传力的定义：衡量的是该特定QTL位点上的基因型差异所揭示的那一部分表型方差。

# 二、Wald检验

## &#9352;检验目标

$$
\mathbf{y}=\mathbf{X}\beta+\mathbf{Z}\gamma+\epsilon
$$

* 检验遗传效应$\gamma$是否显著
* **零假设**$H_0$：$\gamma = 0$（所有$q$个遗传效应均为0）
* **备择假设**$H_1$:$\gamma \neq 0$（至少一个效应非零）

## &#9353;检验统计量构建

### &#9312;原公式：

$$
W=\frac{(\hat\theta-\hat\theta_0)^2}{Var(\hat\theta)}
$$

### &#9313;带入变量后：

$$
W=\hat\gamma^T[Var(\hat\gamma)]^{-1}\hat\gamma^T
$$

* 本质是效应估计量与其他估计误差的加权平方和

### &#9354;统计量的分布与$p$值

* 零假设下：$W$服从自由度为$q$的卡方分布（$\chi^2_q$）
* 其中，$q=dim(\gamma)$
* $p$值计算公式

$$
p=1-F_{\chi^2_q}(W)
$$

* $F_{\chi^2_q}(.)$：自由度为$q$的卡方累计分布函数
* 软件实现

```R
p <- 1 - pchisq(W,q)
```

