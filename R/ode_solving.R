#' @title quasi-dynamic lotka volterra model
#' @param Time vector of time point
#' @param State vector of ODE initial state
#' @param Pars vector for unknown ODE parameters
#' @param power_par matrix of power equation parameters for dependent effect
#' @return list used in ode function
qdODEmod <- function(Time, State, Pars, power_par) {
  nn = length(Pars)
  ind_effect = paste0("alpha","*",names(State)[1])
  dep_effect = sapply(2:nn, function(c) paste0(paste0("beta",c-1),"*",names(State)[c]))
  dep_effect = paste0(dep_effect, collapse = "+")
  all_effect = paste0(ind_effect, "+", dep_effect)
  expr = parse(text = all_effect)

  with(as.list(c(State, Pars)), {
    dx = eval(expr)
    dy <- power_par[,1]*power_par[,2]*Time^(power_par[,2]-1)
    dind = alpha*x
    for(i in c(1:(nn-1))){
      tmp = paste0(paste0("beta",i),"*",paste0("y",i))
      expr2 = parse(text = tmp)
      assign(paste0("ddep",i),eval(expr2))
    }
    return(list(c(dx, dy, dind, mget(paste0("ddep",1:(nn-1))))))
  })
}

#' @title least-square fit for qdODE model
#' @importFrom deSolve ode
#' @param pars vector for unknown ODE parameters
#' @param data data contain independent effect as first row and dependent effect
#' @param Time vector of time point
#' @param power_par matrix of power equation parameters for dependent effect
#' @return mean-square error
qdODE_ls <- function(pars, data, Time, power_par){
  n = length(pars)
  power_par = as.matrix(power_par)
  if (n==2) {
    Pars = c(alpha = pars[1], beta1 = pars[2:n])
    power_par = t(power_par)
    State = c(x=data[1,1],y1 = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  } else{
    Pars = c(alpha = pars[1], beta = pars[2:n])
    State = c(x=data[1,1],y = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  }
  out = as.data.frame(ode(func = qdODEmod, y = State, parms = Pars,
                          times = Time, power_par = power_par))
  X = as.numeric(data[1,])
  fit = as.numeric(out[,2])
  sse = crossprod(X-fit)
  return(sse)
}

#' @title legendre polynomials fit to qdODE model
#' @importFrom deSolve ode
#' @param pars vector of qdODE parameters
#' @param data dataframe of observed data
#' @param Time vector of time point
#' @param power_par matrix of power equation parameters for dependent effect
#' @param LOP_order scalar of LOP order
#' @param new_time vector produce new defined time point
#' @param n_expand scalar for how many interpolation needed
#' @return list contain legendre polynomials parameters, qdODE values and LOP fitted values
qdODE_fit <- function(pars, data, Time, power_par, LOP_order = 6, new_time = NULL, n_expand = 100){
  n = length(pars)
  if (n==2) {
    Pars = c(alpha = pars[1], beta1 = pars[2:n])
    power_par = t(power_par)
    State = c(x=data[1,1],y1 = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  } else{
    Pars = c(alpha = pars[1], beta = pars[2:n])
    State = c(x=data[1,1],y = matrix(data[-1,1], nrow = n-1, ncol=1),ind = data[1,1],dep = rep(0,n-1))
  }
  out = as.data.frame(ode(func = qdODEmod, y = State, parms = Pars,
                           times = Time, power_par = power_par))
  out2 = data.frame(x = out[,1], y = data[1,], y.fit = out[,2],
                    ind = out[,(n+2)], dep = out[,(n+3):(ncol(out))])
  colnames(out2)[4:ncol(out2)] = c(rownames(data)[1], rownames(data)[2:n])
  rownames(out2) = NULL

  all_LOP_par = sapply(2:ncol(out2),function(c)get_legendre_par(out2[,c], LOP_order, out2$x))

  if (is.null(new_time)) {
    time2 = seq(min(Time), max(Time), length = n_expand)
    out3 = apply(all_LOP_par, 2, legendre_fit, x = time2)
    out3 = cbind(time2, out3)
  } else{
    out3 = apply(all_LOP_par, 2, legendre_fit, x = new_time)
    out3 = cbind(new_time, out3)
  }
  colnames(out3) = colnames(out2)
  result = list(fit = out2,
                predict = data.frame(out3),
                LOP_par = all_LOP_par)
  return(result)
}



#' @title wrapper for qdODE model
#' @importFrom stats optim
#' @param result result from power_equation_fit
#' @param relationship list contain variable selection results
#' @param i scalar for which id used for qdODE solving, must <= nrow
#' @param init_pars scalar for initial parameters
#' @param LOP_order scalar of LOP order
#' @param method scalar of qdODE solving methodm, cuurent only support least square
#' @param new_time vector produce new defined time point
#' @param n_expand scalar for how many interpolation needed
#' @param maxit scalar of Optim iteration setting
#' @return list contain variable selection results and LOP parameters for every row
#' @export
qdODE_all <- function(result, relationship, i, init_pars = 1, LOP_order = 6, method = "ls",
                      new_time = NULL, n_expand = 100, maxit = 1e3){
  Time = as.numeric(colnames(result$power_fit))
  variable = c(relationship[[i]]$ind.name, relationship[[i]]$dep.name)
  data = result$power_fit[variable,]

  if (length(variable)<=1) {
    qdODE.est = NA
    result = NA
    return.obj <- append(result, list(ODE.value = NA,
                                      parameters = NA))
   } else{
      power_par = result$power_par[variable,][-1,]
      n = nrow(data)
      pars_int = c(init_pars,relationship[[i]]$coefficient)
      if (method == "ls") {
        qdODE.est <- optim(pars_int, qdODE_ls, data = data, Time = Time, power_par = power_par,
                           method = "BFGS",
                           control = list(trace = TRUE, maxit = maxit))

        result <- qdODE_fit(pars = qdODE.est$par,
                            data = data,
                            power_par = power_par,
                            Time = Time)
        return.obj <- append(result, list(ODE.value = qdODE.est$value,
                                          parameters = qdODE.est$par))
      } else{
        qdODE.est <- optim(pars_int, qdODE_ls, data = data, Time = Time, power_par = power_par,
                           method = "BFGS",
                           control = list(trace = TRUE, maxit = maxit))

        result <- qdODE_fit(pars = qdODE.est$par,
                            data = data,
                            power_par = power_par,
                            Time = Time)
        return.obj <- append(result, list(ODE.value = qdODE.est$value,
                                          parameters = qdODE.est$par))
      }
  }
  return(return.obj)
}


#' @title wrapper for qdODE_all in parallel version
#' @param result result from power_equation_fit
#' @param reduction use n/log(n) dimension reduction
#' @param thread scales for how many threads used
#' @param maxit scalar of Optim iteration setting
#' @return list contain variable selection results and LOP parameters for every row
#' @export
qdODE_parallel <- function(result, reduction = FALSE, thread = 2, maxit = 1e3){
  data = result$original_data
  relationship = lapply(1:nrow(data),function(c)get_interaction(data, c, reduction = reduction))
  core.number <- thread
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterEvalQ(cl, {require(orthopolynom)})
  clusterEvalQ(cl, {require(deSolve)})
  clusterExport(cl, c("qdODEmod", "qdODE_ls", "qdODE_fit", "qdODE_all","get_legendre_matrix",
                      "get_legendre_par","legendre_fit"), envir=environment())
  result = parLapply(1:nrow(data),function(c) qdODE_all(result = result,
                                                       relationship = relationship,
                                                       i = c,
                                                       maxit = maxit
  ), cl = cl)
  stopCluster(cl)
  names(result) = rownames(data)
  names(relationship) = rownames(data)
  return_obj <- list(ode_result = result,
                     relationship = relationship)
  return(return_obj)
}

#' @title plot single decompose plot
#' @import ggplot2
#' @importFrom reshape2 melt
#' @param result list of qdODE all
#' @return effect curve decompose plot
#' @export
qdODE_plot_base <- function(result){

  data = result$predict
  n = ncol(data)
  colnames(data)[4:n] = c(paste0("ind.",colnames(data)[4]),
                          paste0("dep.",colnames(data)[5:n]))

  plot.df = melt(data, id.vars = c("x"))

  name = levels(plot.df[,2])

  ind.name = name[grep("ind", name)]
  ind.name2 = strsplit(ind.name,split = "\\.")[[1]][2]
  ind.df <- subset(plot.df, plot.df[,2] == ind.name)
  ind.df$type = "ind"
  ind.df$variable = ind.name2

  depname = levels(plot.df[,2])[grep("dep",name )]
  dep.df <- subset(plot.df, plot.df[,2] %in% depname)
  dep.df$type = "dep"
  dep.df$variable = sapply(strsplit(as.character(dep.df$variable),"\\."),"[",2)


  original.df = subset(plot.df, plot.df[,2] == "y")
  original.df$type = "original"

  fit.df = subset(plot.df, plot.df[,2] == "y.fit")
  fit.df$type = "fit"

  plot.df2 = rbind(ind.df, dep.df,fit.df)

  name.df = subset(plot.df2, plot.df[,1] == max(plot.df2[,1]))
  name.df = name.df[-nrow(name.df),]
  name.df[,2][name.df[,2] == "y.fit"] = ind.name2

  name.df = name.df[-which(name.df[,4] == "fit"),]

  name.df[,1] = name.df[,1]*1.002

  p <- ggplot(plot.df2, mapping = aes_string(x = "x", y = "value")) +
    geom_line(mapping = aes_string(group = "variable", colour = "type"), size = 1.1) +
    geom_text(name.df, mapping = aes_string(label = "variable", colour = "type",
                                     x = "x", y = "value"), show.legend = FALSE) +
    geom_hline(yintercept = 0, size = 0.5,linetype = 'dashed') +
    scale_x_continuous(limits = c(min(plot.df2$x), max(plot.df2$x)*1.005))+
    scale_color_manual(
      name = "Effect Type",
      labels = c("Dep Effect", "All Effect", "Ind Effect"),
      values = c("green", "blue", "red")) +
    xlab("Habitat Index") + ylab("Microbial Abundance") +
    ggtitle(ind.name2) + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

  return(p)
}

#' @title plot all decompose plot
#' @import ggplot2 patchwork
#' @param result list of qdODE parallel
#' @return all effect curve decompose plot
#' @export
qdODE_plot_all <- function(result){
  #options(warn = -1)
  p = lapply(result$ode_result, qdODE_plot_base)
  p = lapply(p, "+", xlab(NULL))
  p = lapply(p, "+", ylab(NULL))

  y_lab <- ggplot() + annotate(geom = "text", size=7,x = 1, y = 1,
                               label = "Microbial Abundance", angle = 90) +
    coord_cartesian(clip = "off") + theme_void()
  pp = (y_lab | wrap_plots(p)) +
    plot_annotation(caption = "Habitat Index",
                    theme = theme(plot.caption = element_text(size = 20,hjust=.5))) +
    plot_layout(widths = c(.05, 1),guides = 'collect')

  return(pp)
}
