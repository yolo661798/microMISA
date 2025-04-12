library(logging)


DLASSO_fun <- function(X,Y){
  library(glmnet)
  n <- dim(X)[1]
  p <- dim(X)[2]
  fit = glmnet::glmnet(X, Y, alpha = 1)
  cv.fit <- cv.glmnet(X, Y, alpha = 1)
  beta_0 <- coef(fit, s = cv.fit$lambda.min)[2:(p+1)]
  #
  fit <- glmnet::glmnet(X[,2:p], X[,1], alpha = 1)
  cv.fit <- cv.glmnet(X[,2:p], X[,1], alpha = 1)
  phi_hat <- coef(fit, s = cv.fit$lambda.min)[2:p]
  ##
  R <- X[,1] - X[,2:p]%*%t(t(phi_hat))
  E <- Y - X%*%t(t(beta_0))
  beta_1_hat <- beta_0[1] + sum(R*E)/sum(R*X[,1]) #  The de-biased lasso estimator
  ##
  sigma_e2 <- sum(E^2)/(n-length(which(beta_0!=0)))
  
  sigma_beta1_hat <- sqrt(sigma_e2)*sqrt(sum(R^2))/abs(sum(R*X[,1]))
  
  results <- c(beta_1_hat,sigma_beta1_hat)
  return(results)
  
}

## 计算FDP和TPP
getRealFDP = function(real_id,rej_id)
  if(length(rej_id)==0){
    return(0)
  }else{
    return(1 - (length(intersect(rej_id,real_id))/length(rej_id)))
  }

getRealTPP = function(real_id,rej_id)
  if(length(rej_id)==0){
    return(0)
  }else{
    return(length(intersect(rej_id,real_id))/length(real_id))
  }

## ilr变换
# 几何平均
geoMean = function(v)
{
  logv = log(v)
  logmean = mean(logv)
  return(exp(logmean))
}

ilrTrans = function(M)
{
  M = M/sum(M)
  d = length(M)
  lnM = log(M)
  cumsumLnM = cumsum(lnM)
  cumsumLnM_inv = sum(lnM) - cumsumLnM[1:(d-1)]
  alpha1 = (d-1):1
  alpha2 = sqrt(1/(((d-1):1)*(d:2)))
  M_ilr = alpha2*(alpha1*lnM[1:(d-1)]-cumsumLnM_inv)
  return(M_ilr)
}

# 重排序，并进行ilr变换
ilrTrans_l = function(M,l)
{
  M_l = c(M[l],M[-l])
  return(ilrTrans(M_l))
}

## 参数估计
# 线性回归参数估计（非截距的第一个系数）
getLm = function(x,y)
{
  mls = summary(lm(y~x))
  mc = mls$coefficients
  return(c(mc[2,1],mc[2,2]))  # 获取第一个自变量的系数估计值和标准误
}

# 计算线性回归第一个系数的标准化的值
getLmScale = function(x,y)
{
  est = getLm(x,y)
  return(est[1]/est[2])
}

# 估计参数a
aEst = function(X,Z,M)
{
  d = dim(M)[2]
  n = dim(M)[1]
  
  aEst = array(0,c(d,2))
  for (i in 1:d)
  {
    M_ilr = t(apply(M,1,function(M){ilrTrans_l(M,i)}))
    M_ilr = scale(M_ilr)
    a = getLm(cbind(X,Z),M_ilr[,1])
    aEst[i,] = a
  }
  return(aEst)
}

# 估计联合模型参数b
bEst = function(X,Z,M,Y)
{
  d = dim(M)[2]
  n = dim(M)[1]
  
  bEst = array(0,c(d,2))
  for (i in 1:d)
  {
    M_ilr = t(apply(M,1,function(M){ilrTrans_l(M,i)}))
    M_ilr = scale(M_ilr)
    b = getLm(cbind(M_ilr,X,Z),Y)
    bEst[i,] = b
  }
  return(bEst)
}

# 估计参数a和b
abEst = function(X, Z, M, Y)
{
  d = dim(M)[2]
  n = dim(M)[1]
  loginfo("联合模型参数估计开始，共有 %d 组参数需要估计", d, logger='参数估计')
  
  aEst = array(0,c(d,2))
  bEst = array(0,c(d,2))
  for (i in 1:d)
  {
    if (i %% 100 == 0)
      loginfo("正在估计第%d组参数", i, logger='参数估计')
    M_ilr = t(apply(M,1,function(M){ilrTrans_l(M,i)}))
    M_ilr = scale(M_ilr)
    a = getLm(cbind(X,Z),M_ilr[,1])
    aEst[i,] = a
    b = getLm(cbind(M_ilr,X,Z),Y)
    bEst[i,] = b
  }
  loginfo("联合模型参数估计完成", logger='参数估计')
  return(list(aEst = aEst, bEst = bEst))
}

abEstDLasso = function(X, Z, M, Y)
{
  d = dim(M)[2]
  n = dim(M)[1]
  loginfo("联合模型参数估计开始，共有 %d 组参数需要估计", d, logger='参数估计')
  
  aEst = array(0,c(d,2))
  bEst = array(0,c(d,2))
  for (i in 1:d)
  {
    if (i %% 100 == 0)
      loginfo("正在估计第%d组参数", i, logger='参数估计')
    M_ilr = t(apply(M,1,function(M){ilrTrans_l(M,i)}))
    M_ilr = scale(M_ilr)
    a = getLm(cbind(X,Z),M_ilr[,1])
    aEst[i,] = a
    b = DLASSO_fun(cbind(M_ilr,X,Z),Y)
    bEst[i,] = b
  }
  loginfo("联合模型参数估计完成", logger='参数估计')
  return(list(aEst = aEst, bEst = bEst))
}

# 估计边缘模型参数b
bMarEst = function(X,Z,M,Y)
{
  d = dim(M)[2]
  n = dim(M)[1]
  
  bMarEst = array(0,c(d,2))
  for (i in 1:d)
  {
    M_ilr = t(apply(M,1,function(M){ilrTrans_l(M,i)}))
    M_ilr = scale(M_ilr)
    b = getLm(cbind(M_ilr[,1],X,Z),Y)
    bMarEst[i,] = b
  }
  return(bMarEst)
}

# 估计参数a和边缘参数b
abMarEst = function(X, Z, M, Y)
{
  d = dim(M)[2]
  n = dim(M)[1]
  loginfo("边缘模型参数估计开始, 共有 %d 组参数需要估计", d, logger='SIS')
  
  aEst = array(0,c(d,2))
  bMarEst = array(0,c(d,2))
  for (i in 1:d)
  {
    if (i %% 100 == 0)
      loginfo("正在估计第%d组边缘参数", i, logger='SIS')
    M_ilr = t(apply(M,1,function(M){ilrTrans_l(M,i)}))
    M_ilr = scale(M_ilr)
    a = getLm(cbind(X,Z),M_ilr[,1])
    aEst[i,] = a
    b = getLm(cbind(M_ilr[,1],X,Z),Y)
    bMarEst[i,] = b
  }
  loginfo("边缘模型参数估计完成", logger='SIS')
  return(list(aEst = aEst, bMarEst = bMarEst))
}


## SIS
# get index after SIS
getSIS_id = function(iseq,fn,fk = 1){# iseq: indicator sequence; d = k*n/log(n)
  
  fp = length(iseq)
  fq = 1-fk*fn/(fp*log(fn))
  fL = quantile(iseq,fq)
  id = which(iseq >= fL)
  
  # 取 n - 1 大的
  # id = order(iseq, decreasing = TRUE)[1: (fn / 2 - 2)]
  # id = sort(id)
  return(id)
} 

# screening_SIS
screening_SIS = function(X,Z,M,Y,a_est,b_marEst,k = 1)
{
  loginfo("开始SIS变量初步筛选", logger='SIS')
  d = dim(M)[2]
  n = dim(M)[1]
  
  # 待筛选的中介变量指标
  ind_SIS = abs((a_est[,1]/a_est[,2])*(b_marEst[,1]/b_marEst[,2]))
  
  # 筛选
  SIS_id = getSIS_id(ind_SIS,n,k)
  loginfo("SIS变量初步筛选结果：%s", toString(SIS_id), logger='SIS')
  
  return(SIS_id)
}