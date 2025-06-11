setwd('D:/biowolf/guo')
library(bkmrhat)
library(coda)
library(ggplot2)
library(readxl)
future::plan(strategy = future::multisession,workers=20)
# 读取数据
dat <- read.csv('金属.csv')
# 定义变量
metals <- c("锌", "铁", "锰", "钙", "镁", "铅", "铜")   #暴露因素
covariates <- c("age", "BMI")                           #协变量 需为矩阵
y <- dat$optimal_embryo                                 # 二分类结局 (良胚: 1/0)
# 数据预处理
#Z <- scale(dat[, metals])         # 标准化金属浓度
X <- dat[, covariates]            # 协变量矩阵
set.seed(123)                     # 保证可复性
fitkm <- kmbayes_parallel(nchains=20,
                 y = y, Z = Z, X = X, 
                 family = "binomial", 
                 iter = 5000, 
                 verbose = FALSE, 
                 varsel = TRUE)
#多线程结果合并
fitkmcomb=kmbayes_combine(fitkm)
#验证模型收敛性 Investigate model convergence
TracePlot(fit = fitkmcomb, par = "beta")#β值
TracePlot(fit = fitkmcomb, par = "sigsq.eps")#
TracePlot(fit = fitkmcomb, par = "r", comp = 1)
#估计后验包含概率 Estimated posterior inclusion probabilities
ExtractPIPs(fitkmcomb)
#估计h Estimating h
#该模型较小，故采用approach2 
#取μh(θ)后验样本的均值，估计后验均值为E[μh(θ)]，取Vh(θ)后验样本的均
#值与μh(θ)后验样本的方差之和，估计后验方差为E[Vh(θ)]+Var[μh(θ)].
med_vals <- apply(Z, 2, median)
Znew <- matrix(med_vals, nrow = 1)
# 定义HFun函数
HFun <- function(x) {
  # 这里可以实现你的计算逻辑
  # 例如，假设我们计算输入向量的总和
  return(sum(x))  # 确保返回一个数值
}

# 直接调用HFun函数
h_true <- HFun(Znew)
h_true <- dat$HFun(Znew)#直接运行会提示 错误，不适用于非函数，不要运行此行

h_est1 <- ComputePostmeanHnew(fitkmcomb, Znew = Znew, method = "approx")
h_est2 <- ComputePostmeanHnew(fitkmcomb, Znew = Znew, method = "exact")
set.seed(666)
samps3 <- SamplePred(fitkmcomb, Znew = Znew, Xnew = cbind(0,0))


h_est_compare <- data.frame(
  method = c("truth", 1:3),
  post_mean = c(h_true, h_est1$postmean, h_est2$postmean, mean(samps3)),
  post_sd = c(NA, sqrt(h_est1$postvar), sqrt(h_est2$postvar), sd(samps3))
)
h_est_compare #取method3
############################
pred.resp.univar <- PredictorResponseUnivar(fit = fitkmcomb)
ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se,
                             ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z)")
###########################
pred.resp.bivar <- PredictorResponseBivar(fit = fitkmcomb, min.plot.dist = 1)
ggplot(pred.resp.bivar, aes(z1, z2, fill = est)) + 
  geom_raster() + 
  facet_grid(variable2 ~ variable1) +
  scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) +
  xlab("expos1") +
  ylab("expos2") +
  ggtitle("h(expos1, expos2)")
############################
pred.resp.bivar.levels <- PredictorResponseBivarLevels(
  pred.resp.df = pred.resp.bivar, 
  Z = Z, qs = c(0.1, 0.5, 0.9))
ggplot(pred.resp.bivar.levels, aes(z1, est)) + 
  geom_smooth(aes(col = quantile), stat = "identity") + 
  facet_grid(variable2 ~ variable1) +
  ggtitle("h(expos1 | quantiles of expos2)") +
  xlab("expos1")
############################################# 
# 检查和转换 X
if (!is.matrix(X)) {
  X <- as.matrix(X)
}

# 检查 beta
if (!is.vector(beta)) {
  beta <- as.vector(beta)  # 确保 beta 是一个向量 可不运行此条
}

# 运行 OverallRiskSummaries
risks.overall <- OverallRiskSummaries(fit = fitkmcomb, y = y, Z = Z, X = X, 
                                      qs = seq(0.25, 0.75, by = 0.05), 
                                      q.fixed = 0.5, method = "exact")
risks.overall
ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, 
                          ymax = est + 1.96*sd)) + 
  geom_pointrange()
################################################ 
risks.singvar <- SingVarRiskSummaries(fit = fitkmcomb, y = y, Z = Z, X = X, 
                                      qs.diff = c(0.25, 0.75), 
                                      q.fixed = c(0.25, 0.50, 0.75),
                                      method = "exact")
risks.singvar
ggplot(risks.singvar, aes(variable, est, ymin = est - 1.96*sd, 
                          ymax = est + 1.96*sd, col = q.fixed)) + 
  geom_pointrange(position = position_dodge(width = 0.75)) + 
  coord_flip()

##################################################
risks.int <- SingVarIntSummaries(fit = fitkmcomb, y = y, Z = Z, X = X, 
                                 qs.diff = c(0.25, 0.75), 
                                 qs.fixed = c(0.25, 0.75),
                                 method = "exact")
risks.int
###################################################
###################################################
set.seed(111)
d2 <- SimData(n = 100, M = 4, Zgen = "corr", sigsq.true = 2.2)
round(cor(d2$Z), 2)


set.seed(111)
fitkm_corr <- kmbayes(y = d2$y, Z = d2$Z, X = d2$X, iter = 10000, varsel = TRUE, 
                      verbose = FALSE)
fitkm_hier <- kmbayes(y = d2$y, Z = d2$Z, X = d2$X, iter = 10000, varsel = TRUE, 
                      groups = c(1,2,1,3), verbose = FALSE)
ExtractPIPs(fitkm_corr)
ExtractPIPs(fitkm_hier)
#####################################################未通
data.frame(fitkm$control.params)

priorfits <- InvestigatePrior(y = y, Z = Z, X = X, 
                              q.seq = c(2, 1/2, 1/4, 1/16))
PlotPriorFits(y = y, Z = Z, X = X, 
              fits = priorfits)



