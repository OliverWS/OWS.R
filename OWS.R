#o.residuals
require(psych)
require(xtable)
require(ascii)
require(mvnormtest)
require(lawstat)
require(foreign)
require(MASS)
require(car)
require(lsr)
require(gvlma)
require(lattice)
require(nlme)


xyplot.lme <- function(formula,data) {
  xyplot(formula,data=data,panel=function(x,y){panel.xyplot(x,y); panel.lmline(x,y);})
}


ICClme <- function(out) {
  varests <- as.numeric(VarCorr(out)[1:2])
  return(paste("ICC =", varests[1]/sum(varests)))
}



o.lm.xtable <- function(model){
  tbl <- o.lm.table(model)
  return(xtable(tbl,digits=3,align="r|lllll"))
}



o.lm.table <- function(model) {
  cnames <- c("B","SE","Beta","t","p")
  rnames <- names(model$coefficients)
  result = matrix(nrow=length(rnames),ncol=length(cnames))
  model.sm <- summary(model)
  model.stdcoef <- standardCoefs(model)
  for(i in seq_len(length(rnames))){
    result[i,1] <- as.numeric(model.sm$coefficients[i,1])
    result[i,2] <- as.numeric(model.sm$coefficients[i,2])
    if(rnames[i] != "(Intercept)" && i > 1){
      result[i,3] <- as.numeric(model.stdcoef[i-1,2])
    }
    else {
      result[i,3] <- NA
    }
    result[i,4] <- as.numeric(model.sm$coefficients[i,3])
    result[i,5] <- as.numeric(model.sm$coefficients[i,4])
  }
  
  dimnames(result) <- list(rnames,cnames)
  return(result)
  
}





o.barplot <- function(freqs, main="Distribution",xlab="Group",ylab="Frequency",xpd=T) {
  range <- c(min(freqs)*0,max(freqs)*1.1)
  bp <- barplot(freqs, main=main,xlab=xlab, ylab=ylab, xpd=F,ylim=range)
  text(bp,freqs,labels=freqs,pos=3,offset=.5)

}

o.scatter <- function(x,y,main="Scatter Plot w/ LOWESS Fit") {
  plot(x,y,main=main, xlab=deparse(substitute(x)),ylab=deparse(substitute(y)))
  lines(lowess(x,y,f=0.75,iter=5),lwd=3,col="red")
}

o.rmcov <- function(model, variable) {
  mod.mat <- model.matrix(model)
  var.names <- colnames(mod.mat)
  var <- which(variable == var.names)
  if (0 == length(var))
    stop(paste(variable, "is not a column of the model matrix."))
  response <- model.response(model$model)
  responseName <- names(model$model)[[1]]
  if (is.null(weights(model)))
    wt <- rep(1, length(response))
  else wt <- weights(model)
  res <- lsfit(mod.mat[, -var], cbind(mod.mat[, var], response),
               wt = wt, intercept = FALSE)$residuals
  colnames(res) <- c(var.names[var], responseName)
  rownames(res) <- rownames(mod.mat)
  
  return(res)
}


std <- function(x, na.rm=T) {
  mu <- mean(x,na.rm=na.rm)
  s <- sd(x,na.rm=na.rm)
  return((x-mu)/s)
  
}

o.ppplot <- function(dat, main = "P-P Plot", xlab = "Observed Probability", ylab = "Expected Probability", ...) {
  probDist <- pnorm(dat)
  plot(ppoints(length(dat)), sort(probDist), main=main, xlab=xlab, ylab=ylab, ...)
  abline(0,1)
}

o.outliers <- function(model, cutoff=3.0) {
  return(which(abs(stdres(model)) >= cutoff))
}

o.residuals <- function(model,  cd.threshold = 1.0) {
  askSetting = par("ask")
  #Now we do our stuff
  data  = model$model
  yname = strsplit(names(data)[[1]],"\\$")[[1]][2]
  xname = strsplit(names(data)[[2]],"\\$")[[1]][2]
  y = data[[1]]
  x = data[[2]]
  res = residuals(model)
  res.std = stdres(model)
  outs = o.outliers(model)
  b0 = coef(model)[[1]]
  b1 = coef(model)[[2]]
  
  #First, lets plot a histogram
  par(mfrow=c(1,2))
  o.hist(res.std)
  
  #Now we'll do a pp plot
  o.ppplot(res.std)
  par(ask=T)

  #Now lets do a plot showing ourliers
  par(mfrow=c(1,1))
  o.plot.residuals(model)
  print(outs)
  outlier_table <- as.data.frame(cbind(outs,y[outs],x[outs],res.std[outs]))
  View(outlier_table)
  colnames(outlier_table) <- c("index",yname,xname,"Standardized Residual")
  print(outlier_table)
  print(xtable(outlier_table, caption="Outliers"))
  
  #Now compare to a model fit without the outliers, adding an abline
  alt_x = x
  alt_x[outs] = NA
  alt_y = y
  alt_y[outs] = NA
  alt_model = lm(alt_y~alt_x)
  print(summary(alt_model))
  abline(alt_model,lwd=2,col="darkblue")
  legend("topleft",c("Model","Model Fit w/out Outliers","Residuals","Outliers"),col=c("black","darkblue","darkgreen","red"), inset=0.01,pch=c(3,3,16,16))
  #Now lets run some diagnostics on influence
  
  res.cd = cooks.distance(model)
  print(describe(res.cd))
  if(length(which(res.cd > cd.threshold)) > 0) {
    res.cd.x = x[which(res.cd > cd.threshold)]
    res.cd.y = y[which(res.cd > cd.threshold)]
    res_table = as.data.frame(cbind(res.cd.y,res.cd.x,res.std[which(res.cd > cd.threshold)], res.cd[which(res.cd > cd.threshold)]) )
    names(res_table) <- list(yname,xname,"Standardized Residual","Cook's Distance")
    print(res_table)
    print(xtable(res_table))
  }
  else {
    print(sprintf("There were no residuals with Cook's Distance greater than %f",cd.threshold))
  }

  #Restore settings
  par(ask=askSetting)
  par(mfrow=c(1,1))
}
skew.se <- function(x){
  n = length(x)
  se = 6/sqrt(n)
  return(se)
}

skew.ratio <- function(x){return(skew(x)/skew.se(x))}


o.plot.lm <- function(model, col="darkgreen", ...) {
  res = residuals(model)
  res.cd = cooks.distance(model)
  data  = model$model
  yname = names(data)[[1]]
  ncoef = length(model$coefficients)
  nrows = ceiling(sqrt(ncoef-1))
  par(mfrow=c(nrows,nrows))
  for(fact in seq(2,ncoef)){
    xname = names(data)[[fact]]
    y = data[[1]]
    x = data[[fact]]
    plot(x,y,pch=16,col=col,type="p",xlab=xname,ylab=yname, ...)
    outs = o.outliers(model)
    abline(coef(model)[[1]],coef(model)[[fact]])
    points(x[outs],y[outs],col="red",pch=16,...)
    points(x[which(res.cd > 1)],y[which(res.cd > 1)],col="yellow",type="p", pch=0,... )
    points(x[which(res.cd > 4)],y[which(res.cd > 4)],col="red",type="p",pch=2,... )
  }
  par(mfrow=c(1,1))
  
}

o.plot.residuals <- function(model, col="darkgreen", ...) {
  res = residuals(model)
  res.cd = cooks.distance(model)
  data  = model$model
  yname = strsplit(names(data)[[1]],"\\$")[[1]][2]
  xname = strsplit(names(data)[[2]],"\\$")[[1]][2]
  y = data[[1]]
  x = data[[2]]
  plot(x,y,pch=16,col=col,type="p",xlab=xname,ylab=yname, ...)
  outs = o.outliers(model)
  abline(model,lwd=2)
  points(x[outs],y[outs],col="red",pch=16,...)
  points(x[which(res.cd > 1)],y[which(res.cd > 1)],col="yellow",type="p", pch=0,... )
  points(x[which(res.cd > 4)],y[which(res.cd > 4)],col="red",type="p",pch=2,... )
}



stem <- function(x,y,mdl,pch=16,linecol=1,clinecol=1,showAll=T,...){
  x = as.numeric(as.character(x))
  xf = as.factor(x)
  if (missing(y)){
    print("No Y value supplied")
    y = x
    x = 1:length(x) }
  c = coefficients(mdl)[[1]]
  m = coefficients(mdl)[[2]]
  if(showAll == F){
    px = vector(length=2*length(levels(xf)))
    py = vector(length=2*length(levels(xf)))
    nl = 1
    for (l in levels(xf)){
      for(yi in c(min(y[x==l]),max(y[x==l]))){
        xl = as.numeric(l)
        px[nl] = xl
        py[nl] = yi
        nl = nl + 1
      }
    }
  }
  else {
    px = x
    py = y
  }
  
  plot(px,py,pch=pch,type="p",...)
  for (l in levels(xf)){
    for(yi in c(min(y[x==l]),max(y[x==l]))){
      xl = as.numeric(l)
      lines(c(xl,xl), c(c + xl*m, yi),col=linecol,lwd=2)
    }
  }
  abline(mdl,col=clinecol,lw=2)
  #lines(c(x[1]-2,x[length(x)]+2), c(0,0),col=clinecol)
}


omega.squared <- function (model, sseffects=NULL) {
  if (class(model)[1]=="aov") {
    df.factors=0
    for (i in 2:dim(model$model)[2]) {
      df.factors=df.factors+length(levels(model$model[,i]))-1
    }
  }
  if (class(model)[1]=="lm") {
    df.factors=summary(model)$fstatistic[2]
  }
  
  nvariables = dim(model$model)[2]
  nobservations = dim(model$model)[1]
  ncontrasts = 0
  for (i in 2:nvariables) {
    ncontrasts = ncontrasts + length(as.factor(levels(as.factor(model$model[,i]))))-1
  }
  endeffects = ncontrasts+1
  startresiduals = ncontrasts+2
  sseffects = sum(model$effects[2:endeffects]^2)
  
  ssresiduals = sum(model$effects[startresiduals:nobservations]^2)
  
  omega.sq = (sseffects-(df.factors*(ssresiduals/model$df.residual)))/(sseffects+ssresiduals+(ssresiduals/model$df.residual))
  attr(omega.sq,"names") <- NULL
  return(omega.sq)
}

pcor.test <- function(x,y,z,use="mat",method="p",na.rm=T){
  # The partial correlation coefficient between x and y given z
  #
  # pcor.test is free and comes with ABSOLUTELY NO WARRANTY.
  #
  # x and y should be vectors
  #
  # z can be either a vector or a matrix
  #
  # use: There are two methods to calculate the partial correlation coefficient.
  #	 One is by using variance-covariance matrix ("mat") and the other is by using recursive formula ("rec").
  #	 Default is "mat".
  #
  # method: There are three ways to calculate the correlation coefficient, 
  #	    which are Pearson's ("p"), Spearman's ("s"), and Kendall's ("k") methods.
  # 	    The last two methods which are Spearman's and Kendall's coefficient are based on the non-parametric analysis.
  #	    Default is "p".
  #
  # na.rm: If na.rm is T, then all the missing samples are deleted from the whole dataset, which is (x,y,z).
  #        If not, the missing samples will be removed just when the correlation coefficient is calculated.
  #	   However, the number of samples for the p-value is the number of samples after removing 
  #	   all the missing samples from the whole dataset.
  #	   Default is "T".
  
  x <- c(x)
  y <- c(y)
  z <- as.data.frame(z)
  
  if(use == "mat"){
    p.use <- "Var-Cov matrix"
    pcor = pcor.mat(x,y,z,method=method,na.rm=na.rm)
  }else if(use == "rec"){
    p.use <- "Recursive formula"
    pcor = pcor.rec(x,y,z,method=method,na.rm=na.rm)
  }else{
    stop("\'use\' should be either \"rec\" or \"mat\"!\n")
  }
  
  # print the method
  if(gregexpr("p",method)[[1]][1] == 1){
    p.method <- "Pearson"
  }else if(gregexpr("s",method)[[1]][1] == 1){
    p.method <- "Spearman"
  }else if(gregexpr("k",method)[[1]][1] == 1){
    p.method <- "Kendall"
  }else{
    stop("\'method\' should be \"pearson\" or \"spearman\" or \"kendall\"!\n")
  }
  
  # sample number
  n <- dim(na.omit(data.frame(x,y,z)))[1]
  
  # given variables' number
  gn <- dim(z)[2]
  
  # p-value
  if(p.method == "Kendall"){
    statistic <- pcor/sqrt(2*(2*(n-gn)+5)/(9*(n-gn)*(n-1-gn)))
    p.value <- 2*pnorm(-abs(statistic))
    
  }else{
    statistic <- pcor*sqrt((n-2-gn)/(1-pcor^2))
    p.value <- 2*pnorm(-abs(statistic))
  }
  
  data.frame(estimate=pcor,p.value=p.value,statistic=statistic,n=n,gn=gn,Method=p.method,Use=p.use)
}			

# By using var-cov matrix
pcor.mat <- function(x,y,z,method="p",na.rm=T){
  
  x <- c(x)
  y <- c(y)
  z <- as.data.frame(z)
  
  if(dim(z)[2] == 0){
    stop("There should be given data\n")
  }
  
  data <- data.frame(x,y,z)
  
  if(na.rm == T){
    data = na.omit(data)
  }
  
  xdata <- na.omit(data.frame(data[,c(1,2)]))
  Sxx <- cov(xdata,xdata,m=method)
  
  xzdata <- na.omit(data)
  xdata <- data.frame(xzdata[,c(1,2)])
  zdata <- data.frame(xzdata[,-c(1,2)])
  Sxz <- cov(xdata,zdata,m=method)
  
  zdata <- na.omit(data.frame(data[,-c(1,2)]))
  Szz <- cov(zdata,zdata,m=method)
  
  # is Szz positive definite?
  zz.ev <- eigen(Szz)$values
  if(min(zz.ev)[1]<0){
    stop("\'Szz\' is not positive definite!\n")
  }
  
  # partial correlation
  Sxx.z <- Sxx - Sxz %*% solve(Szz) %*% t(Sxz)
  
  rxx.z <- cov2cor(Sxx.z)[1,2]
  
  rxx.z
}

# By using recursive formula
pcor.rec <- function(x,y,z,method="p",na.rm=T){
  # 
  
  x <- c(x)
  y <- c(y)
  z <- as.data.frame(z)
  
  if(dim(z)[2] == 0){
    stop("There should be given data\n")
  }
  
  data <- data.frame(x,y,z)
  
  if(na.rm == T){
    data = na.omit(data)
  }
  
  # recursive formula
  if(dim(z)[2] == 1){
    tdata <- na.omit(data.frame(data[,1],data[,2]))
    rxy <- cor(tdata[,1],tdata[,2],m=method)
    
    tdata <- na.omit(data.frame(data[,1],data[,-c(1,2)]))
    rxz <- cor(tdata[,1],tdata[,2],m=method)
    
    tdata <- na.omit(data.frame(data[,2],data[,-c(1,2)]))
    ryz <- cor(tdata[,1],tdata[,2],m=method)
    
    rxy.z <- (rxy - rxz*ryz)/( sqrt(1-rxz^2)*sqrt(1-ryz^2) )
    
    return(rxy.z)
  }else{
    x <- c(data[,1])
    y <- c(data[,2])
    z0 <- c(data[,3])
    zc <- as.data.frame(data[,-c(1,2,3)])
    
    rxy.zc <- pcor.rec(x,y,zc,method=method,na.rm=na.rm)
    rxz0.zc <- pcor.rec(x,z0,zc,method=method,na.rm=na.rm)
    ryz0.zc <- pcor.rec(y,z0,zc,method=method,na.rm=na.rm)
    
    rxy.z <- (rxy.zc - rxz0.zc*ryz0.zc)/( sqrt(1-rxz0.zc^2)*sqrt(1-ryz0.zc^2) )
    return(rxy.z)
  }			
}	

o.pcor <- function(x, method="p", na.rm=T,latex=T) {
  #Loop through all the possible co-variate combinations to create matrix
  x <- as.data.frame(x)
  variables <- colnames(x)
  ncols <- length(variables)
  zdims <- c("r","p.value")
  output <- array(dim=c(ncols,ncols,length(zdims)),dimnames=list(variables,variables,zdims))
  for(row in variables) {
    for(col in variables) {
      covariates <- variables[which(variables != row)]
      covariates <- covariates[which(covariates != col)]
      pc <- pcor.test(x[row],x[col],as.data.frame(x[covariates]),method=method,na.rm=na.rm)
      output[row,col,"r"] <- pc[["estimate"]]
      output[row,col,"p.value"] <- pc[["p.value"]]
    }
  }
  print("Partial Correlations")
  print(output[,,"r"],digits=3)
  print("p-values")
  print(output[,,"p.value"],digits=3)
  if(latex){
    print(xtable(output[,,"r"],digits=3,caption="Partial Correlation Coefficients"))
    
    print(xtable(output[,,"p.value"],digits=3,caption="P-Values"))
    
  }
}


o.describe <- function (x, digits = 2,na.rm=TRUE, type = 3, ...)   #basic stats after dropping non-numeric data
{                                   #first, define a local function
  valid <- function(x) {      
    return(sum(!is.na(x)))
  }
  if (is.vector(x) )          #do it for vectors or 
  {
    stats = matrix(rep(NA,7),ncol=7)    #create a temporary array
    stats[1, 1] <-  valid(x )
    stats[1, 2] <-  mean(x, na.rm=na.rm )
    stats[1, 3] <-  sd(x, na.rm=na.rm )
    stats[1, 4] <-  min(x, na.rm=na.rm )
    stats[1, 5] <-  max(x, na.rm=na.rm )
    stats[1, 6] <-  skew(x,na.rm=na.rm)
    stats[1, 7] <-  kurtosi(x,na.rm=na.rm)
    len <- 1;
  }
  
  
  else  {
    len = dim(x)[2]      #do it for matrices or data.frames
    if(is.null(len)){len = 1}
    stats = matrix(rep(NA,len*7),ncol=7)    #create a temporary array
    for (i in 1:len) {
      if ((len==1 && is.numeric(x)) || is.numeric(x[,i])) {   #just do this for numeric data
        stats[i, 1] <-  valid(x[,i] )
        stats[i, 2] <-  mean(x[,i], na.rm=na.rm )
        stats[i, 3] <-  sd(x[,i], na.rm=na.rm )
        stats[i, 4] <-  min(x[,i], na.rm=na.rm )
        stats[i, 5] <-  max(x[,i], na.rm=na.rm )
        stats[i, 6] <-  skew(x[,i], na.rm=na.rm)
        stats[i, 7] <-  kurtosi(x[,i], na.rm=na.rm)
      }
    }
  }
  temp <-  data.frame(n = stats[,1],mean = stats[,2], sd = stats[,3], min= stats[,4],max=stats[,5],skew = stats[, 6],kurtosis=stats[,7])
  answer <-  round(data.frame(temp),  digits) #, se = temp$sd/sqrt(temp$n)
  rownames(answer) <- colnames(x)
  return(answer)
}

o.describe.by <- function (x, group = NULL, mat = FALSE, type = 3) 
{
  if (is.null(group)) {
    answer <- o.describe(x, type = type)
    warning("no grouping variable requested")
  }
  else {
    answer <- by(x, group, o.describe, type = type)
    c(answer,Total=o.describe(x))
  }
  if (mat) {
    ncol <- length(answer[[1]])
    n.var <- nrow(answer[[1]])
    n.col <- ncol(answer[[1]])
    n.grouping <- length(dim(answer))
    n.groups <- prod(dim(answer))
    names <- names(answer[[1]])
    row.names <- attr(answer[[1]], "row.names")
    dim.names <- attr(answer, "dimnames")
    mat.ans <- matrix(NA, ncol = ncol, nrow = n.var * n.groups)
    labels.ans <- matrix(NA, ncol = n.grouping + 1, nrow = n.var * 
                           n.groups)
    colnames(labels.ans) <- c("item", paste("group", 1:n.grouping, 
                                            sep = " "))
    colnames(mat.ans) <- colnames(answer[[1]])
    rn <- 1:(n.var * n.groups)
    k <- 1
    labels.ans[, 1] <- seq(1, (n.var * n.groups))
    group.scale <- cumprod(c(1, dim(answer)))
    for (var in 1:(n.var * n.groups)) {
      for (group in 1:n.grouping) {
        groupi <- ((trunc((var - 1)/group.scale[group]))%%dim(answer)[group]) + 
          1
        labels.ans[var, group + 1] <- dim.names[[group]][[groupi]]
      }
    }
    k <- 1
    for (var in 1:n.var) {
      for (group in 1:n.groups) {
        rn[k] <- paste(row.names[var], group, sep = "")
        for (stat in 1:n.col) {
          if (!is.null(answer[[group]][[stat]][var])) {
            mat.ans[k, stat] <- answer[[group]][[stat]][var]
          }
          else {
            mat.ans[k, stat] <- NA
          }
        }
        k <- k + 1
      }
    }
    answer <- data.frame(labels.ans, mat.ans)
    rownames(answer) <- rn
  }
  return(answer)
}



o.describeBy <- function(x,groups,tex=TRUE,digits=3) {
  desc <- describeBy(x,groups,mat=TRUE)
  if(tex){
    print(xtable(desc[,2:ncol(desc)]),digits=digits,include.rownames=F)
  }
  else {
    print(desc)
  }
  
}

o.hist <- function(x, group=NULL, xlab=deparse(substitute(x)), main=paste("Histogram with Normal Curve: ",deparse(substitute(x))), save=FALSE, ...){
  if(!is.null(group)){
    nplots = length(levels(group))
    rows = nplots/2
    if(nplots <= 2){rows=1}
    cols = 2
    par(mfrow=c(rows,cols))
    groups = split(x,group)
    namelist = names(group)
    for(n in seq_len(length(groups)) ) {
      name <- names(groups[n])[1]
      lapply(groups[n],o.hist,main=name,xlab=xlab)
      print(name)
    }
  }
  else {
    # Add a Normal Curve (Thanks to Peter Dalgaard)
    h<-hist(x, col="light gray", xlab=xlab, 
            main=main, ...) 
    x = na.omit(x)
    xfit<-seq(min(x,na.rm=na.omit),max(x,na.rm=na.omit),length=length(x)) 
    yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
    yfit <- yfit*diff(h$mids[1:2])*length(x) 
    lines(xfit, yfit, col="black", lwd=3)
    if(save) {
      
    }
  }
}

o.assumptions <- function(y,group=NULL,type="ANOVA",sig=0.05,n_threshold=30) {
  result = TRUE
  #========================Test Normality==========================
  if(is.null(group)){
    isnorm <- o.isnorm(y,sig=sig)
    if(!isnorm){result = FALSE}
  }
  else {
    by(y,group,function(x){  
      isnorm <- o.isnorm(x,sig=sig)
      if(!isnorm){result = FALSE}
      
    })
  }

  #========================Test Homogeneity of Variance==========================
  if(is.null(group)){return}
  hov <- o.hov(y,group,sig=sig)
  if(!hov){result = FALSE}
}

o.isnorm <- function(y,sig=0.05,n_threshold=30) {
  print("Testing Assumption of Normality Using Shapiro–Wilk",quote=FALSE)
  if(length(y) < 3){
    print("Shaprio test requires N > 2",quote=FALSE)
    return(FALSE)
  }
  isnorm <- shapiro.test(y)
  print(isnorm)
  if(isnorm["p.value"] < sig){
    print(sprintf("Since p=%f is less than threshold of %0.2f, we can reject H0 and conclude that distribution IS NOT normal.",isnorm["p.value"],sig), quote=FALSE)
    if(length(y) > n_threshold) {
      print(sprintf("However, since N=%d is larger than %d we can assume that the distribution of means will be normal via the Central Limit Theorem.", length(y), n_threshold), quote=FALSE)
      return(TRUE)
    }
    else {
      return(TRUE)
    }
  }
  else {
    print(sprintf("Since p=%0.3f is greater than threshold of %0.3f, we fail to reject H0 and conclude that distribution IS normal.",isnorm["p.value"],sig), quote=FALSE)
    return(TRUE)
  }
}

o.hov <- function(y,group,sig=0.05) {
  print("Testing Homogeneity of Variance Using Levene's Test")
  homogeny <- leveneTest(y, group)
  print(homogeny)
  if(homogeny["p.value"] < sig){
    print(sprintf("Since p=%f is less than threshold of %0.2f, we can reject H0 and conclude that the variances are NOT EQUAL.",homogeny["p.value"],sig), quote=FALSE)
    return(FALSE)
  }
  else {
    print(sprintf("Since p=%0.3f is greater than threshold of %0.3f, we fail to reject H0 and conclude that the variances ARE EQUAL",homogeny["p.value"],sig), quote=FALSE)
    return(TRUE)
  }
}

o.aov <- function(formula) {
  fit <- aov(formula)
  drop1(fit,~.,test="F") #Correct for the fact that R does Type I ANOVA by default; we want Type III
  print(summary(fit))
  return(fit)
}

o.anova <- function(formula,data=NULL) {
  vars <- get_all_vars(formula,data=data)
  y = vars[,1]
  groups = vars[,2]:vars[,3]
  assumptions_met <- o.assumptions(y,groups)
  #if(!assumptions_met){return(FALSE)}
  print.table(o.anova.table(formula,TRUE),quote=FALSE,digits=2)
  
}

o.ancova <- function(model,sig_cutoff=FALSE)
{
  fit <- model
  sig_levels <- c(0.001,0.01,0.05)
  anova_info <- Anova(fit)
  sum_sqs <- anova_info["Sum Sq"]
  totalSS = sum(sum_sqs)
  cnames <- c("df","Type III Sum of Squares","eta^2","F Statistic","Significance")
  rnames <- row.names(anova_info)
  result = matrix(nrow=length(rnames),ncol=length(cnames))
  for(i in seq_len(nrow(anova_info))){
    ss <- anova_info[i,"Sum Sq"]
    result[i,1] <- anova_info[i,"Df"]
    result[i,2] <- ss
    result[i,3] <- ss/totalSS
    result[i,4] <- anova_info[i,"F value"]
    p <- anova_info[i,"Pr(>F)"]
    if(!is.na(p) && !sig_cutoff){
      if(p < 0.001) {
        p = 0.000
      }
      result[i,5] <- signif(p,3)
    }
    else if(!is.na(p)){
      result[i,5] <- sprintf("p = %0.3f",p)
      for(sig in sort(sig_levels)) {
        if(p<sig){
          result[i,5] <- sprintf("p < %s",sig)
          break
        }
      }
    }
    else {
      result[i,5] <- NA
    }
  }
  
  dimnames(result) <- list(rnames,cnames)
  return(result)
}

o.anova.table <- function(formula,sig_cutoff=FALSE)
{
  fit <- o.aov(formula)
  sig_levels <- c(0.001,0.01,0.05)
  anova_info <- anova(fit)
  sum_sqs <- anova_info["Sum Sq"]
  totalSS = sum(sum_sqs)
  cnames <- c("df","Type III Sum of Squares","eta^2","F Statistic","Significance")
  rnames <- row.names(anova_info)
  result = matrix(nrow=length(rnames),ncol=length(cnames))
  for(i in seq_len(nrow(anova_info))){
    ss <- anova_info[i,"Sum Sq"]
    result[i,1] <- anova_info[i,"Df"]
    result[i,2] <- ss
    result[i,3] <- ss/totalSS
    result[i,4] <- anova_info[i,"F value"]
    p <- anova_info[i,"Pr(>F)"]
    if(!is.na(p) && !sig_cutoff){
      if(p < 0.001) {
        p = 0.000
      }
      result[i,5] <- signif(p,3)
    }
    else if(!is.na(p)){
      result[i,5] <- sprintf("p = %0.3f",p)
      for(sig in sort(sig_levels)) {
        if(p<sig){
          result[i,5] <- sprintf("p < %s",sig)
          break
        }
      }
    }
    else {
      result[i,5] <- ""
    }
  }
  
  dimnames(result) <- list(rnames,cnames)
  return(result)
}