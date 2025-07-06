# 一、引入和约束条件

```R
# 定义矩阵 X 和向量 y
X <- matrix(c(
  1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 1, 1, 1,
  1, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 1, 0, 1, 1, 1, 0,
  0, 0, 0, 1, 0, 0, 0, 1
), nrow = 8)

y <- matrix(c(129, 127, 174, 196, 160, 154, 181, 190), nrow = 8)

# 定义矩阵 H
H <- matrix(c(
  0, 1, 1, 0, 0, 0,
  0, 0, 0, 1, 1, 1
), nrow = 2, byrow = TRUE)

# 提取 X_r（假设 X_r 是 X 的某些列）
X_r <- X[, c(-2, -4)]  # 去掉第 2 列和第 4 列

# 计算 beta
X_r_transpose <- t(X_r)          # X_r 的转置
X_r_X <- X_r_transpose %*% X   # 正确地计算 X_r' * X_r

# 确保 H 的列数与 X_r_X 的列数一致
if (ncol(H) != ncol(X_r_X)) {
  stop("H 的列数必须与 X_r_X 的列数一致")
}

# 构建 combined_matrix
combined_matrix <- rbind(X_r_X, H)  # 将 X_r_X 和 H 按行绑定

# 检查是否可逆并求逆
if (det(combined_matrix) == 0) {
  stop("combined_matrix 不可逆！请检查输入矩阵。")
}
inverse_matrix <- solve(combined_matrix)

# 构建 [X_r'; 0]
zero_matrix <- matrix(0, nrow = nrow(H), ncol = ncol(X_r_transpose))
X_r_transpose_zero <- rbind(X_r_transpose, zero_matrix)

# 计算 T = inverse_matrix %*% X_r_transpose_zero
T <- inverse_matrix %*% X_r_transpose_zero

# 计算 beta_hat = T %*% y
beta_hat <- T %*% y

# 计算 V_hat_beta = T %*% t(T) * sigma^2
sigma_squared <- 14.675
V_hat_beta <- T %*% t(T) * sigma_squared
print(V_hat_beta)
# 标准化残差矩阵
S_hat_beta <- sqrt(diag(V_hat_beta))  # 开平方得到标准误差
S_beta_inv <- diag(1 / S_hat_beta)    # 对角线元素取倒数
R_beta <- S_beta_inv %*% V_hat_beta %*% S_beta_inv
print(R_beta)
# 输出结果
list(
  beta_hat = beta_hat,
  V_hat_beta = V_hat_beta,
  R_beta = R_beta
)
```

# 二、Harvey线性约束条件下

```R
# 定义矩阵 X 和向量 y
X <- matrix(c(
  1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 1, 1, 1,
  1, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 1, 0, 1, 1, 1, 0,
  0, 0, 0, 1, 0, 0, 0, 1
), nrow = 8)

y <- matrix(c(129, 127, 174, 196, 160, 154, 181, 190), nrow = 8)

last_col <- X[, ncol(X)] * -1

# 将最后一列加到倒数第二列和倒数第三列
X[, ncol(X) - 1] <- X[, ncol(X) - 1] + last_col
X[, ncol(X) - 2] <- X[, ncol(X) - 2] + last_col

# 删除最后一列
X <- X[, -ncol(X)]

#提取第三列
third_col <- X[, 3] * -1

# 将第三列加到第二列
X[,2]<-X[,2]+ third_col

X<-X[,-3]
#最小二乘估计方程beat
beta <- solve(t(X) %*% X) %*% t(X) %*% y
print(beta)

```

