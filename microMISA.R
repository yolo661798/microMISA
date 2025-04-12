source("utils.R")
library(logging)


## MISA
# 求fdp曲线
getFdp = function(W)
{
  W_abs = abs(W)
  
  ord = order(W_abs,decreasing = T)
  W = W[ord]
  W_abs = W_abs[ord]
  
  # fdp用于保存所有t\inW_abs时的FDP
  fdp = lapply(W_abs, function(x){(1 + length(which(W<=-x))) / max(length(which(W>=x)),1)})
  return(list(W_abs = as.numeric(W_abs), fdp = as.numeric(fdp)))
}

# 确定阈值
getThre = function(W, q){ # W should be vector or 1-col matrix
  result = getFdp(W)
  W_abs = getFdp(W)$W_abs
  fdp = getFdp(W)$fdp
  
  q_real = max(q, min(fdp))
  thre = W_abs[max(which(fdp <= q_real))]
  return(thre)
}

# 获得检验统计量
getWk = function(X,Z,M,Y,SIS_id,seed = 1)
{
  loginfo("计算检验统计量W开始", logger='microMISA')
  # 参数设定
  n = dim(M)[1] # 样本数
  p = dim(M)[2] # 中介变量维数
  MS = M[,SIS_id] # 经筛选的中介变量子集
  d = length(SIS_id)  # 经筛选的中介变量个数
  
  # 样本分割
  set.seed(seed)
  samp1_id = sort(sample(1:n,n/2))
  samp2_id = setdiff(1:n,samp1_id)
  loginfo("样本分割完成", logger='microMISA')
  
  # 参数估计
  loginfo("估计划分1的参数开始", logger='microMISA')
  est1 = abEst(X[samp1_id],Z[samp1_id,],MS[samp1_id,],Y[samp1_id])
  #est1 = abEstDLasso(X[samp1_id],Z[samp1_id,],MS[samp1_id,],Y[samp1_id])
  a_samp1 = est1$aEst
  b_samp1 = est1$bEst
  
  loginfo("估计划分2的参数开始", logger='microMISA')
  est2 = abEst(X[samp2_id],Z[samp2_id,],MS[samp2_id,],Y[samp2_id])
  #est2 = abEstDLasso(X[samp2_id],Z[samp2_id,],MS[samp2_id,],Y[samp2_id])
  a_samp2 = est2$aEst
  b_samp2 = est2$bEst
  
  #a1 = a_samp1[,1] / a_samp1[,2]
  #a2 = a_samp2[,1] / a_samp2[,2]
  #b1 = b_samp1[,1] / b_samp1[,2]
  #b2 = b_samp2[,1] / b_samp2[,2]
  
  ind_Mirror = a_samp1[,1]*a_samp2[,1]*b_samp1[,1]*b_samp2[,1] / (a_samp1[,2]*a_samp2[,2]*b_samp1[,2]*b_samp2[,2])
  Wk = rep(0,p)
  Wk[SIS_id] = ind_Mirror
  
  loginfo("计算检验统计量W完成", logger='microMISA')
  return(Wk)
}


# 拒绝的变量下标
getRejId = function(Wk,q = 0.1)
{
  p = length(Wk)
  L = getThre(Wk,q)
  loginfo("q = %f, 检验统计量W的筛选阈值为 %f", q, L, logger='microMISA')
  rej_id = (1:p)[which(Wk>=L)]
  return(rej_id)
}


# microMISA的一次检验
mmisa = function(X,Z,M,Y,qs = 0.1, SIS = TRUE, seed = 1)
{
  loginfo("----------------------------------------", logger='microMISA')
  loginfo("microMISA开始", logger='microMISA')

  if (SIS == TRUE)
  {
    est = abMarEst(X,Z,M,Y)
    a_est = est$aEst
    b_marEst = est$bMarEst
    SIS_id = screening_SIS(X,Z,M,Y,a_est,b_marEst)
  }
  else
    SIS_id = 1:(dim(M)[2])
  
  Wk = getWk(X,Z,M,Y,SIS_id,seed)
  
  rej_ids = list()
  for (q in qs)
  {
    rej_id = getRejId(Wk,q)
    rej_ids[[as.character(q)]] = rej_id
  }
  
  out_result = list(ID = rej_ids, Wk = Wk)
  loginfo("microMISA完成", logger='microMISA')
  return(out_result)
}


# bagging
getRejId_bagging = function(X,Z,M,Y,B,qs,SIS_id)
{
  # 参数设定
  n = dim(M)[1] # 样本量
  d = dim(M)[2] # 中介变量维数
  MS = M[,SIS_id] # 经筛选的中介变量M的子集
  dS = length(SIS_id)  # 经筛选的中介变量M的个数
  
  # 用于存储每一次检验结果的矩阵rej
  rejs = array(0, dim = c(length(qs), B, d), dimnames = list(as.character(qs), NULL, NULL))
  
  # 重复B次实验，记录每次实验哪些变量被拒绝
  for (i in 1:B)
  {
    loginfo("【第 %d 次重复实验】", i, logger='microMISA-bagging')
    # 构造对称样本并计算Wk
    Wk = getWk(X,Z,M,Y,SIS_id,i)
    # 拒绝原假设的下标
    for (q in qs)
    {
      rej_id = getRejId(Wk,q)
      rejs[as.character(q), i, rej_id] = 1
    }
  }
  
  loginfo("整合实验结果开始", logger='microMISA-bagging')
  rej_ids = list()
  for (q in qs)
  {
    rej = rejs[as.character(q),,]
    # 计算每个中介变量被拒绝的次数
    rej_times = apply(rej,2,sum)
    
    # 求被拒绝次数>=B/2的变量
    rej_idMost = which(rej_times>=(ceiling(B/2)))
    unRej_idMost = which(rej_times<(ceiling(B/2)))
    
    # 求最终被拒绝的变量
    overlap = rep(0,B)
    for (i in 1:B)
    {
      rej_id = which(rej[i,]>0)
      unRej_id = setdiff((1:d), rej_id)
      overlap[i] = length(intersect(rej_id,rej_idMost)) + length(intersect(unRej_id,unRej_idMost))
    }
    k = which(overlap == max(overlap))[1]
    
    rej_id = which(rejs[as.character(q), k,] > 0)
    
    rej_ids[[as.character(q)]] = rej_id
  }
  
  out_result = list(ID = rej_ids, Wk = Wk, all_IDs = rejs)
  loginfo("microMISA-bagging结束", logger='microMISA-bagging')
  return(out_result)
}

# 封装bagging
mmisa_bagging = function(X,Z,M,Y,B = 20,qs = 0.1,SIS = TRUE)
{
  loginfo("----------------------------------------", logger='microMISA-bagging')
  loginfo("microMISA-bagging开始", logger='microMISA-bagging')
  
  # 参数设定
  n = dim(M)[1] # 样本量
  d = dim(M)[2] # 中介变量维数
  
  if (SIS == TRUE)
  {
    a_est = aEst(X,Z,M)
    b_marEst = bMarEst(X,Z,M,Y)
    SIS_id = screening_SIS(X,Z,M,Y,a_est,b_marEst)
  } 
  else 
    SIS_id = 1:d
  
  out_result = getRejId_bagging(X,Z,M,Y,B,qs,SIS_id)
  return(out_result)
}

