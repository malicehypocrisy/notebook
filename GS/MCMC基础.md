# 一、MCMC基础

## &#9352;拒绝性采样

### &#9312;原理

通过使用工具分布`g(x)`从任意概率分布函数`f(x)`中模拟样本，唯一的限制是`f(x)< M g(x)`，其中`M>1`对`f(x)/g(x)`的适当界限。
$$
p(x|accepted)=\frac{p(accepted|x_i)g(x_i)}{\sum(accepted|x_i)g(x_i)}=\frac{[f(x_i)/cg(x_i)]g(x_i)}{\sum{\{[f(x_i)/cg(x_i)]g(x_i)\}}}=f(x_i)
$$

### &#9313;算法

#### &#9332;需求

* 由已知分布的概率密度函数$f(x)$，产生服从此分布的样本X。

#### &#9333;准备工作：

* 需要一个辅助的“建议分布G”来产生候选样本。
* 一个辅助的均匀分布$U(0,1)$。
* 计算常数值`c`。满足不等式$c*g(x)\geq f(x)$的最小值c（接近于1）

#### &#9334;流程

* 步骤1：从建议分布$G$抽样，得到样本$Y$.

* 步骤2：从分布$U(0,1)$抽样，得到样本$U$

* 步骤3：如果$U\leq\frac{f(Y)}{c*g(Y)}$，则令$X=Y$（接受Y），否则继续执行步骤1（拒绝）

## &#9353;Metropolis-Hastings算法

### &#9312;流程

* 初始化时间`t=1`

* 设置`u`的值，并初始化初始状态$\theta^{(t)} =u$

* 重复一下进程

  * 令$t=t+1$

  * 从已知分布中$q(\theta|\theta^{(t-1)})$中生成一个候选状态$\theta^{(*)}$

  * 计算接受这个的概率：

  * $$
    \alpha=min(1,\frac{p(\theta^{(*)})}{p(\theta^{(t-1)}}\frac{p(\theta^{(t-1)}|\theta^{(*)})}{\theta^{(*)})|p(\theta^{(t-1)}})
    $$

  * 从分布$U(0,1)$生成一个随机数`a`

  * 如果$a\leq\alpha$,接受新生成的值：$\theta^{(t)}=\theta^{(*)}$，否则：$\theta^{(t)}=\theta^{(t-1)}$。

* 直到$t=T$

## &#9354;Metropolis 算法

特点：它是Metropolis-Hastings采样器的一种特殊情况
$$
g(X^*|X_t)=g(X_t|X^*)
\\
\alpha=min(1,\frac{f(X^*)}{f(X)})
$$

## &#9355;原理

$$
\mathbf{\pi}^{(n+1)}=\mathbf{\pi^{(n)}}P
$$

* 转移矩阵$\mathbf{P}$是一个$N\times N$矩阵，其中元素$P_{ij}=P(i,j)$

### &#9312;转移概率（核）

* 从状态$s_k$到状态$s_j$的转移概率为：

$$
P(k,j)=pr(X_{n+1}=s_j|X_n=s_k)
$$

### &#9313;边际概率

* 在时间n的边际概率为：

$$
\pi^{(n)}_j=pr(X_n=s_j)
$$

## &#9356;吉布斯采样

### &#9312;算法流程

* **先验超参数**：在循环外（固定假设）
* **先验对后验的影响**：在循环内的条件后验分布公式中（每次抽样都要融合下先验和数据）

### &#9313;具有非信息先验的正态分布

* 无信息先验的正态分布

$$
p(\theta, \sigma^2) = p(\theta) p(\sigma^2) = C \times \frac{1}{\sigma^2} \propto \frac{1}{\sigma^2}
\\
p(\theta) = \text{Uniform}(\min, \max) \quad p(\log(\sigma^2)) = \text{Uniform}(0, 1)
$$

* 条件后验分布

$$
\sigma^2|\theta,\mathbf{y}\sim\chi^{-1}(n,\frac{\sum{(y_i-\theta)^2}}{n})
\\
\theta|\sigma^2,\mathbf{y}\sim N(\overline{y},\frac{\sigma^2}{n})
$$

### &#9314;正态分布与半共轭先验

* 先验分布

$$
\theta|\sigma^2 \sim N\left(\mu_0, \tau_0^2\right)
\\ \sigma^2 \sim Inv-\chi^2\left(\nu_0, S_0^2\right)
$$

* 后验分布

$$
\theta|\sigma^2, y \sim N\left(\mu_n, \sigma_n^2\right)
\\\sigma^2|\theta, y \sim \chi^{-2}\left(\nu_0 + n, \frac{\nu_0 S_0^2 + \sum_{i=1}^n (y_i - \theta)^2}{\nu_0 + n}\right)
\\
\text{where:}

\mu_n = \left(\frac{1}{\tau_0^2} + \frac{n}{\sigma^2}\right)^{-1} \left(\frac{1}{\tau_0^2} \mu_0 + \frac{n}{\sigma^2} \overline{y}\right)
\\
\frac{1}{\sigma_n^2} = \frac{1}{\tau_0^2} + \frac{n}{\sigma^2}
$$



### &#9315;具有共轭先验的正态数据

* 先验分布

$$
\theta | \sigma^2 \sim N\left( \mu_0, \frac{\sigma^2}{\kappa_0} \right), \quad \sigma^2 \sim \chi^{-2}\left( \nu_0, \tau_0^2 \right)
$$

* 后验分布

$$
\sigma^2 | \theta, y \sim \chi^{-2}\left( n + \nu_0 + 1, \frac{\sum_{i=1}^n (y_i - \theta)^2 + \nu_0 \tau_0^2 + \kappa_0 (\theta - \mu_0)^2}{n + \nu_0 + 1} \right)
\\
\theta | \sigma^2, y \sim N\left( \mu_n, \frac{\sigma^2}{\kappa_n} \right)
\\
\kappa_n = \kappa_0 + n, \quad \mu_n = \frac{\kappa_0}{\kappa_0 + n} \mu_0 + \frac{n}{\kappa_0 + n} \overline{y} 
$$

# 二、统计基础

## &#9352;二项分布（binomial distribution)

$$
X\sim B(n, p)\\
Pr(y|n,p) = \frac{n!}{y!(n-y)!}p^y(1-p)^{n-y}
\\
E[y|n, p] = np 
\\
Var[y|n, p] = np(1-p)
$$

* 当满足条件$nq\geq5$，二项分布接近正态分布$N(np,np(1-p))$

## &#9353;多项式分布

* 随机实验结果不再是两种状态，而是k种互斥的离散状态，每种状态出现的概率是$p_i$
* $p_1 +P_2 +p_3+...+p_k=1$

$$
\mathbf{y}|n,\mathbf{p}\sim Multin(n,p_1,p_2,...p_k)
\\
Pr(\mathbf{y}|n,\mathbf{p})=\frac{n!}{y_1!y_2!...y_K!}p_1^{y_1}p_2^{y_2}...p_k^{y_k}
\\E[y_i|n, p] = np_i 
\\
Var[y_i|n, p] = np_i(1-p_i)
\\
Cov[y_i,y_j|n,p]=-np_ip_j
$$

## &#9354;泊松分布

* 描述小概率事件在一定时间或空间范围内的发生次数的概率分布。

$$
Pr(y|\mathbf{\mu,\sigma^2})=\frac{\lambda^ye^{-\lambda}}{y!}
\\
E[y|\lambda]=Var[y|\lambda]=\lambda
$$

## &#9355;正态分布

$$
p(y|\mu,\sigma^2)=\frac{1}{\sqrt{2\pi\sigma^2}}exp[-\frac{1}{2\sigma^2}(y-\mu)^2]
\\
E(y|\mu,\sigma^2)=\mu,Var(y|\mu,\sigma^2)=\sigma^2
$$

### &#9312;多元正态分布

$$
y|\mu,\sum\sim N_p(\mu,\sum)\begin{cases}{-\infty} <\mu< {+\infty}\\
\sum & \text{正定}\\
{1\infty}<y{+\infty}
\end{cases}
\\
p(y|\mu,\sum)=(2\pi)^{-p/2}|\sum|^{-1/2}exp[-\frac{1}{2}(y-\mu)^T\sum{^{-1}}(y-\mu)]

\\
\Sigma = 
\begin{bmatrix}
\Sigma_{11} & \Sigma_{12} \\
\Sigma_{21} & \Sigma_{22}
\end{bmatrix}
$$

* 条件分布的均值向量为 $\mu_1 + \Sigma_{12} \Sigma_{22}^{-1} (a - \mu_2)$
* 条件分布的协方差矩阵为 $\Sigma_{11} - \Sigma_{12} \Sigma_{22}^{-1} \Sigma_{21}$

## &#9356;均匀分布

$$
p(y|a,b)=\frac{1}{b-a}
\\
E(y|a,b)=\frac{a+b}{2};Var(y|a,b)=\frac{(b-a)^2}{12}
$$

## &#9357;贝塔分布

* 关于连续变量$x \subseteq[0,1]$的概率分布，由两个参数$a>0$和$b>0$确定：

$$
p(y|a,b)=\frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)}y^{a-1}(1-y)^{b-1}
\\
E[y|a,b]=\frac{a}{a+b};Var[y|a,b]=\frac{ab}{(a+b)^2(a+b+1)}
$$

### &#9312;狄利克雷分布

* 贝塔分布向多变量的推广，常用于对**概率向量**进行建模。

* 若随机变量：$\mathbf{X}=(X_1,X_2,…,X_3)$满足：

$$
\begin{cases}
X_i&\geq 0
\\
\sum^K_{i=1}X_i&=1
\end{cases}
$$

*  则称$\boldsymbol{X}$ 服从参数为$\boldsymbol{\alpha} = (\alpha_1, \alpha_2, \dots, \alpha_K)$的狄利克雷分布

$$
\mathbf{X}\sim Dirichlet(\alpha)
$$

* **先验参数 + 本轮观测到的计数**

## &#9358;伽马分布

$$
p(y|a,b)=\frac{b^a}{\Gamma(a)}y^{a-1}exp(-by)
\\
E(y|a,b)=\frac{a}{b};Var(y|a,b)=\frac{a}{b^2}
$$

### &#9312;伽马函数

$$
\Gamma(x)=\int_0^{+\infty}u^{x-1}e^{-u}du
$$

### &#9313;逆伽马分布

$$
f(x;\alpha,\beta)=\frac{\beta^\alpha}{\Gamma(\alpha)}x^{-(\alpha+1)}exp(-\frac{\beta}{x})
$$

##  &#9359;卡方分布

$$
y|\phi \sim \chi^2_{\phi}\begin{cases}\phi>0\\
y>0
\end{cases}
\\
p(y|\phi)=\frac{2^{-\phi/2}}{\Gamma(\phi/2)}y^{\phi/2-1}exp{(-y/2)}
\\
E(y|\phi)=\phi;Var(y|\phi)=2\phi
$$

* 卡方分布式伽马分布的特列

$$
\chi^2_{\phi}=Gamma(\frac{\phi}{2},\frac{1}{2})
$$

## &#9360;指数分布

$$
p(y|\lambda)=\lambda e^{-\lambda y}
\\
E(y|\lambda)=\lambda^{-1};E(y|\lambda)=\lambda^{-2}
$$

* 是伽马分布的特列 

$$
Exp(\lambda)=Gamma(1,\lambda)
$$

## &#9361;Wishart 分布

$$
W|\nu,S\sim Wishart_\mu(S)\begin{cases}\nu:\text{自由度}\\
S:\text{正定对称矩阵};(k\times k)
\end{cases}
\\
p(W \mid \nu, S) \propto W^{-(\nu - k - 1)/2} \exp\left\{ -\frac{1}{2} \operatorname{tr}(S^{-1} W) \right\}
\\
E(W \mid \nu, S) \propto \nu S
$$

## &#9362;Inverse Wishart 分布

$$
W \mid \nu, S \sim \text{Inv-Wishart}_\nu(S^{-1})\begin{cases}\nu:\text{自由度}\\
S:\text{正定对称矩阵};(k\times k)
\end{cases}
\\
p(W \mid \nu, S) \propto W^{-(\nu + k + 1)/2} \exp\left\{ -\frac{1}{2} \operatorname{tr}(S W^{-1}) \right\}
\\
E(W \mid \nu, S) \propto (\nu - k - 1)^{-1} S
$$

# 三、贝叶斯数据分析

$$
p(\theta|y)=\frac{p(\theta,y)}{p(y)}
$$

## &#9352;模型评估与预测

### &#9312;点预测

$$
MSE=\frac{1}{2}\sum^n_{i=1}{(y_i-E(y_i|\theta))^2}
$$

### &#9313;加强版本

$$
MSE=\frac{1}{2}\sum^n_{i=1}{\frac{(y_i-E(y_i|\theta))^2}{Var(y_i|\theta)}}
$$

* 对不确定预测给予更大的容忍度

### &#9314;调整后样本内预测精度

#### &#9332;AIC

$$
\text{AIC} = -2 \log p(y|\hat{\theta}_{\text{MLE}}) + 2k
$$

* 惩罚项：$k$是模型中参数的个数。参数越多，模型越复杂，惩罚越大。

#### &#9333;DIC

$$
\text{DIC} = -2 \log p(y|\hat{\theta}{\text{Bayes}}) + 2p{\text{DIC}}
$$

* $+ 2p_{\text{DIC}}$：**有效参数数量（惩罚项）**，用于衡量模型的复杂度。

## &#9353;贝叶斯因子

$$
\frac{p(M_2 \mid y)}{p(M_1 \mid y)} = \frac{p(M_2)}{p(M_1)} \times \text{BF}(M_2, M_1)
\\
\text{BF}(M_2, M_1) = \frac{p(y \mid M_2)}{p(y \mid M_1)}

= \frac{\int p(\theta_2 \mid M_2) p(y \mid \theta_2, M_2) d\theta_2}{\int p(\theta_1 \mid M_1) p(y \mid \theta_1, M_1) d\theta_1}
$$

