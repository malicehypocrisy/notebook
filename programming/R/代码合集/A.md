```R
# 定义函数构建亲缘关系矩阵A
build_A_matrix <- function(ped) {
  # ped: 数据框，包含列"个体", "父", "母"，按出生顺序排列
  n <- nrow(ped)
  ids <- ped$个体
  A <- matrix(0, nrow = n, ncol = n)
  rownames(A) <- colnames(A) <- as.character(ids)
  
  # 标记基础群（双亲未知的个体）
  base <- ped$个体[is.na(ped$父) & is.na(ped$母)]
  
  # 初始化基础群个体
  for (i in base) {
    A[as.character(i), as.character(i)] <- 1
  }
  
  # 处理非基础群个体
  #setdiff(ids, base) 会返回那些在向量 ids 中存在但不在向量 base中出现的元素
  non_base <- setdiff(ids, base)
  for (i in non_base) {
    i_id <- as.character(i)
    p <- ped$父[ped$个体 == i]    # 父编号
    m <- ped$母[ped$个体 == i]    # 母编号
    
    # 获取父和母的行向量（若未知则设为0）
    row_p <- if (!is.na(p)) A[as.character(p), ] else rep(0, n)
    row_m <- if (!is.na(m)) A[as.character(m), ] else rep(0, n)
    
    # 计算当前个体的行：0.5*(父行 + 母行)
    current_row <- 0.5 * (row_p + row_m)
    A[i_id, ] <- current_row
    A[, i_id] <- current_row  # 对称赋值
    
    # 计算对角线元素：1 + 0.5*父母亲缘系数
    a_parents <- if (is.na(p) || is.na(m)) 0 else A[as.character(p), as.character(m)]
    A[i_id, i_id] <- 1 + 0.5 * a_parents
  }
  
  return(A)
}

# 构建矩阵A
A <- build_A_matrix(A)
print(A)

```

