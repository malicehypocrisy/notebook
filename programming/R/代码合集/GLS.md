# 一、广义二乘估计（GLS）

```R
# 定义矩阵 z
z <- matrix(c(1,1,0,0,1,0,0,0,0,0,1,0,0,1,1,1,0,0,0,1,0,0,0,0),
            nrow = 8)

# 计算 z 的转置与自身相乘
z_y <- z %*% t(z)

# 创建一个与 z_y 大小相同的矩阵，并初始化为 0
z_y_optimized <- matrix(0, nrow = ncol(z_y), ncol = ncol(z_y))

# 设置对角线元素为 625
diag(z_y_optimized) <- 625

# 设置非对角线元素为 49 如果原矩阵 z_y 中对应位置为 1
z_y_optimized[z_y == 1 & !diag(ncol(z_y))] <- 49

# 输出结果
print(z_y_optimized)
```

