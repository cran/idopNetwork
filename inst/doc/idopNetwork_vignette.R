## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  tidy = FALSE,
  cache = TRUE,
  collapse = TRUE,
  dev = "png",
  fig.width = 7, 
  fig.height = 3.5
)

## ----setup--------------------------------------------------------------------
library(idopNetwork)
backup_options <- options()
#load pre-computered results
test_result = idopNetwork:::test_result

## ----eval=FALSE---------------------------------------------------------------
#  data("gut_microbe")
#  View(gut_microbe)

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(gut_microbe[1:10,1:5])

## ---- eval=FALSE--------------------------------------------------------------
#  data("mustard_microbe")
#  View(mustard_microbe)

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(mustard_microbe[1:10,1:5])
knitr::kable(mustard_microbe[1:10,89:93])

## ---- echo=TRUE, eval=TRUE----------------------------------------------------
df = data_cleaning(gut_microbe)

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  result1 = power_equation_fit(df)

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
result1 = test_result$d1_power_fitting

## ---- echo=TRUE, eval=TRUE----------------------------------------------------
power_equation_plot(result1)

## ---- echo=TRUE---------------------------------------------------------------
matplot(t(power_equation(x = 1:30, matrix(c(2,1,3,0.2,0.5,-0.5),nrow = 3, ncol = 2))), 
        type = "l",
        xlab = "time", 
        ylab = "population")
legend("topright", 
       c("cluster 1", "cluster 2", "cluster 3"), 
       lty = c(1,2,3), 
       col = c(1,2,3), 
       box.lwd = 0)

## ---- echo=TRUE---------------------------------------------------------------
get_SAD1_covmatrix(c(2,0.5), n = 5)

## ---- echo=TRUE---------------------------------------------------------------
get_par_int(X = log10(df+1), k = 4, times = as.numeric(log10(colSums(df)+1)))

#use kmeans to get initial centers
tmp = kmeans(log10(df+1),4)$centers
tmp2 = power_equation_fit(tmp, trans = NULL)
power_equation_plot(tmp2, label = NULL,n = 4)


## ---- eval=TRUE---------------------------------------------------------------
options(max.print = 10)
fun_clu(result1$original_data, k = 3, iter.max = 5)

## ---- eval=FALSE--------------------------------------------------------------
#  result2 = fun_clu_parallel(result1$original_data, start = 2, end = 5)

## ---- eval=TRUE, echo=FALSE---------------------------------------------------
result2 = test_result$d1_cluster

## ---- eval=TRUE---------------------------------------------------------------
best.k = which.min(sapply( result2 , "[[" , 'BIC' )) + 1 #skipped k = 1
best.k

fun_clu_BIC(result = result2)

#we can direct give other k value
fun_clu_plot(result = result2, best.k = best.k)

## ---- eval=TRUE---------------------------------------------------------------
data("mustard_microbe")
df2 = data_cleaning(mustard_microbe, x = 160)

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  res_l = power_equation_fit(df2[,1:5]
#  res_r = power_equation_fit(df2[,89:95])
#  res1 = data_match(result1 = res_l, result2 = res_r)

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
res1 = test_result$d2_power_fitting

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  res2 = bifun_clu_parallel(data1 = res1$dataset1$original_data,
#                            data2 = res1$dataset2$original_data,
#                            Time1 = res1$dataset1$Time,
#                            Time2 = res1$dataset2$Time,
#                            start = 2,
#                            end = 10,
#                            thread = 9,
#                            iter.max = 10)

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
res2 = test_result$d2_cluster

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  res2 = bifun_clu_parallel(data1 = res1$dataset1$original_data,
#                            data2 = res1$dataset2$original_data,
#                            Time1 = res1$dataset1$Time,
#                            Time2 = res1$dataset2$Time,
#                            start = 2,
#                            end = 10,
#                            thread = 9,
#                            iter.max = 10)

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
res2 = test_result$d2_cluster

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
fun_clu_BIC(result = res2)

#we can set best.k directly
bifun_clu_plot(result = res2, best.k = 3, color1 = "#C060A1", color2 = "#59C1BD")

## ---- echo=TRUE, eval=FALSE---------------------------------------------------
#  res3 = bifun_clu_convert(res2, best.k = 3)
#  large.module = order(sapply(res3$a$Module.all,nrow))[5]
#  
#  res_suba = fun_clu_select(result_fit = res1$dataset1, result_funclu = res3$a, i = large.module)
#  res_subb = fun_clu_select(result_fit = res1$dataset2, result_funclu = res3$b, i = large.module)
#  dfsuba_l = power_equation_fit(res_suba$original_data)
#  dfsubb_r = power_equation_fit(res_subb$original_data)
#  ressub1 = data_match(result1 = dfsuba_l, result2 = dfsubb_r)
#  ressub2 = bifun_clu_parallel(data1 = ressub1$dataset1$original_data,
#                               data2 = ressub1$dataset2$original_data,
#                               Time1 = ressub1$dataset1$Time,
#                               Time2 = ressub1$dataset2$Time,
#                               start = 2,
#                               end = 5,
#                               iter.max = 3)

## ---- echo==FALSE-------------------------------------------------------------
ressub2 = test_result$d2_subcluster

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
fun_clu_BIC(result = ressub2)
bifun_clu_plot(result = ressub2, best.k = 2, color1 = "#C060A1", color2 = "#59C1BD", degree = 1)

## -----------------------------------------------------------------------------
result3 = fun_clu_convert(result2,best.k = best.k)
df.module = result3$original_data
get_interaction(df.module,1)


## -----------------------------------------------------------------------------
#we can the microbial relationship in Module1
df.M1 = result3$Module.all$`1`
get_interaction(df.M1,1)

## ---- eval=TRUE---------------------------------------------------------------
options(max.print = 10)
# first we test solving a qdODE
module.relationship = lapply(1:best.k, function(c)get_interaction(df.module,c))
ode.test = qdODE_all(result = result3, relationship = module.relationship, 1, maxit = 100)
# we can view the result
qdODE_plot_base(ode.test)

## ---- eval=FALSE--------------------------------------------------------------
#  # then we solve all qdODEs
#  ode.module = qdODE_parallel(result3)

## ---- echo=FALSE--------------------------------------------------------------
ode.module = test_result$d1_module

## ---- eval=TRUE---------------------------------------------------------------
qdODE_plot_all(ode.module)

## ---- eval=FALSE--------------------------------------------------------------
#  result_m1 = fun_clu_select(result_fit = result1, result_funclu = result3, i = 1)
#  ode.M1 = qdODE_parallel(result_m1)

## ---- echo=FALSE--------------------------------------------------------------
ode.M1 = test_result$d1_M1

## ---- eval=TRUE, fig.height=8, fig.width=10-----------------------------------
qdODE_plot_all(ode.M1)

## ---- eval=TRUE, fig.height=8, fig.width=10-----------------------------------
net_module = lapply(ode.module$ode_result, network_conversion)
network_plot(net_module, title = "Module Network")

## ---- eval=TRUE, fig.height=8, fig.width=10-----------------------------------
net_m1 = lapply(ode.M1$ode_result, network_conversion)
network_plot(net_m1, title = "M1 Network")

## ---- eval=FALSE--------------------------------------------------------------
#  mustard_module_a = qdODE_parallel(res3$a)
#  mustard_module_b = qdODE_parallel(res3$b)
#  
#  res_m1a = fun_clu_select(result_fit = res1$dataset1, result_funclu = res3$a, i = 3)
#  res_m1b = fun_clu_select(result_fit = res1$dataset2, result_funclu = res3$b, i = 3)
#  mustard_M1a = qdODE_parallel(res_m1a)
#  mustard_M1b = qdODE_parallel(res_m1b)

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
mustard_module_a = test_result$d2_module[[1]]
mustard_module_b = test_result$d2_module[[2]]

## ---- eval=TRUE---------------------------------------------------------------
mustard_m_a <- lapply(mustard_module_a$ode_result, network_conversion)
mustard_m_b <- lapply(mustard_module_b$ode_result, network_conversion)

#set seed to make same random layout
layout(matrix(c(1,2),1,2,byrow=TRUE))
set.seed(1)
network_plot(mustard_m_a, title = "Module Network a")
set.seed(1)
network_plot(mustard_m_b, title = "Module Network b")

## ---- eval=FALSE--------------------------------------------------------------
#  result_m1a = fun_clu_select(result_fit = res1$dataset1, result_funclu = res3$a, i = 1)
#  result_m1b = fun_clu_select(result_fit = res1$dataset2, result_funclu = res3$b, i = 1)
#  ode.m1a = qdODE_parallel(result_m1a, thread = 16)
#  ode.m1b = qdODE_parallel(result_m1b, thread = 16)

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
ode.m1a = test_result$d2_m1[[1]]
ode.m1b = test_result$d2_m1[[2]]

## ---- eval=TRUE---------------------------------------------------------------
net_m1a = lapply(ode.m1a$ode_result, network_conversion)
net_m1b = lapply(ode.m1b$ode_result, network_conversion)

#set seed to make same random layout
layout(matrix(c(1,2),1,2,byrow=TRUE))
set.seed(1)
network_plot(net_m1a, title = "Module1 a")
set.seed(1)
network_plot(net_m1b, title = "Module1 b")

## ---- eval=TRUE, echo=FALSE---------------------------------------------------
options(backup_options)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

