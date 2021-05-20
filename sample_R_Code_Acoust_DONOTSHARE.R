####  acoust.GMM.regression.r


require(R.matlab)
require(stringr)
require(dplyr)
require(lubridate)
require(glmnet)
require(ggplot2)


dir.mat <- 'C:/Users/mshapiro/Documents/MATLAB/Acoustic Sensors Matlab/Acoustics Matlab Sandbox/Data'
dir.dat <-  './Data'

# d.mpm <- readRDS('d.33.smooth.LOWESS_1min.rds')
d.mpm <- readRDS('d.33.LOESS_1min.rds')

sfx <- 'GMM_Moments_2mx.reg'
files.mat <- dir(path = dir.mat, pattern = sfx, all.files = FALSE,
                 full.names = TRUE)

names.tmp <- character(length = 0)

ix.names <- 1:21
for(ix in ix.names){
  names.tmp <- c(names.tmp,paste0('GMM.MU.1.',ix),paste0('GMM.MU.2.',ix),
                 paste0('GMM.SIG.1.',ix),paste0('GMM.SIG.2.',ix))
}
gmm.names <- c(names.tmp,'GMM.RHO.1', 'GMM.RHO.2')
mu.names <- sapply(ix.names, function(x) paste0('MOMENT.MU.',x))
sig.names <- sapply(ix.names, function(x) paste0('MOMENT.SIG.',x))
skew.names <- sapply(ix.names, function(x) paste0('MOMENT.SKEW.',x))
kurt.names <- sapply(ix.names, function(x) paste0('MOMENT.KURT.',x))

coeff.names <- c(gmm.names,mu.names,sig.names,skew.names,kurt.names,'DATE.TIME')

d.g <- readMat(files.mat[1]) %>% as.data.frame(stringsAsFactors=FALSE)
for(file in files.mat[-1]){
  d.temp <- readMat(file) %>% as.data.frame(stringsAsFactors=FALSE)
  d.g <- rbind(d.g, d.temp)
}
names(d.g) <- coeff.names
d.g$DATE.TIME <- dmy_hms(d.g$DATE.TIME) %>% round_date('minute')
d.g.1 <- d.g[, c('DATE.TIME', head(names(d.g),-1))]

## join gmm and meter.data
d.reg <- dplyr:: left_join(d.g.1, d.mpm, by='DATE.TIME')

saveRDS(d.reg, 'd.reg.moments.rds')

############## elastic net

d.reg <- readRDS('d.reg.moments.rds')
## 70% of the sample size
n.train <- floor(0.7 * nrow(d.reg))
## set the seed to make the partition reproductible
# set.seed(42)
ix.train <- sample(seq_len(nrow(d.reg)), size = n.train)

##### Choose predictor variables
### 'MOMENT.MU.','MOMENT.SIG.', 'MOMENT.SKEW.','MOMENT.KURT.'

sfx.x <- '(MOMENT.MU|MOMENT.SIG)' 
sfx.save <- 'MOMENT.MU.SIG'

names.x <- c(names(d.reg)[grep(sfx.x,names(d.reg))], "WL.PROD.FLOWLINE.PRES", "WL.PROD.FLOWLINE.TEMP")

#### Choose response variables
names.mpm.y <- c("MPM.ACTG.STD.OIL.FLOW.RATE",
                 "MPM.ACTG.STD.H2O.FLOW.RATE",               
                 "MPM.ACTG.STD.GAS.FLOW.RATE", 
                 "MPM.GAS.VOLUMETRIC.FLOW.STD") 
names.tsep.y <- c( "TSEP.STD.OIL.FLOW.RATE" ,
                   "TSEP.STD.WATER.FLOW.RATE",      
                   "TSEP.STD.GAS.FLOW.RATE",
                   "TSEP.GOR")

#### meter  <- either 'mpm' or 'tsep'
meter <- 'mpm'
if(meter=='mpm'){
  names.y <- names.mpm.y ##[1]  # if [i] will run single-variate regression
}else if(meter=='tsep'){
  names.y <- names.tsep.y ##[1]
}else{
  stop(paste0('meter should be either \'mpm\' or \'tsep\' '))
}
time.train <- d.reg$DATE.TIME[ix.train]
time.test <- d.reg$DATE.TIME[-ix.train]

x.train <- d.reg[ix.train, names.x] %>% as.matrix()
x.test <- d.reg[-ix.train, names.x] %>% as.matrix()

y.train <- d.reg[ix.train, names.y] %>% as.matrix()
y.test <- d.reg[-ix.train, names.y] %>% as.matrix()

if(length(names.y)>1){
  mfit <- glmnet(x.train, y.train , family='mgaussian')  ### multivariate
  plot(mfit, xvar = "lambda", label = TRUE, type.coef = "2norm")
  plot(mfit, xvar='dev', label = TRUE, type.coef = "2norm")
}else{
  mfit <- glmnet(x.train, y.train)   ### single
  plot(mfit, xvar = "lambda", label = TRUE)
  plot(mfit, xvar='dev', label = TRUE)
}

mfit

if(length(names.y)>1){
  cv.mfit <- cv.glmnet(x.train, y.train , family='mgaussian')  ### multivariate
}else{
  cv.mfit <- cv.glmnet(x.train, y.train )  ### single
}

plot(cv.mfit)

if(length(names.y)>1){
  yhat <- predict(mfit, x.test , family='mgaussian')  ### multivariate
  
  yerr <- data.frame(matrix(nrow=100, ncol=ncol(y.test)))
  names(yerr) <- dimnames(y.test)[[2]]
  for( ix.name in 1:ncol(y.test)){
    yerr[,ix.name] <- as.vector(apply((yhat[,dimnames(yhat)[[2]][ix.name],] 
                                       - y.test[,ix.name])^2,2,mean))
  }
  # points(log(mfit$lambda), yerr[,2], col='blue', pch='*')
}else{
  yhat <- predict(mfit, x.test)   ### single
}


yhat.cv <- predict.cv.glmnet(cv.mfit, x.test, s="lambda.1se")

if(length(names.y)>1){
  for(i in 1:ncol(y.test)){
    plot(y.test[,i], type='l', main=dimnames(y.test)[[2]][i],  col='green', lwd=2)
    lines(yhat.cv[,i,], col='blue', lwd=2)
  }
}else{
  plot(y.test[,i], type='l', main=dimnames(y.test)[[2]][i],  col='green', lwd=2)
  lines(yhat.cv[,i], col='blue', lwd=2)
}  

##### Finding the Rsquared
Rsquared<-list()
MAPE <- list()
RMSE <- list()

if(length(names.y)>1){
  for(iname in dimnames(y.test)[[2]]){
    ## Rsquared
    error<-sum((y.test[ , iname] - yhat.cv [ ,iname,1])^2)
    tss <- sum((y.test[ ,iname] - mean(y.test[,iname]))^2)
    Rsquared[[iname]] <- 1-(error/tss)
    RMSE[[iname]] <- sqrt(mean(error))
    ## MAPE
    err.mape <- sum(abs((y.test[ , iname] - yhat.cv [ ,iname,1])))
    MAPE[[iname]] <- err.mape / sum(y.test[ , iname]) *100
    
  }
}else{
  error<-sum((y.test - yhat.cv)^2)
  tss <- sum((y.test - mean(y.test))^2)
  Rsquared[[names.y]] <-1-(error/tss)
  RMSE <- sqrt(mean(error))
  MAPE[[names.y]] <- sum(abs((y.test - yhat.cv ))) / sum(y.test) *100
}


# ######PLOTTING CODE
d.test<-data.frame(y.test)

if(length(names.y)>1){
  d.yhat.cv<-data.frame(yhat.cv[,,1])
  
  for(iname in dimnames(d.test)[[2]] ){
    # plot(d.test[,i], type='l', main=paste(dimnames(d.test)[[2]][i],round(Rsquared[[i]],3),sep=" "),  col='blue', lwd=2)
    # lines(d.yhat.cv[,i,],col='red',lwd = 2)
    
    # p <- ggplot(d.test,aes_string(x = 1:nrow(d.test),y = iname)) + geom_line(color='green',size=1) +
    #   geom_line(data = d.yhat.cv,aes_string(x = 1:nrow(d.yhat.cv),y = iname),color = 'blue', size=1) +
    #   ggtitle(paste0(iname, ',  R2 = ' , round(Rsquared[[iname]], digits=2), ' (no WL)'))+ xlab('')
    
    p <- ggplot(d.test,aes(x = time.test,y = d.test[iname])) + geom_line(color='green',size=1) +
      geom_line(data = d.yhat.cv, aes(x = time.test,y = d.yhat.cv[iname]),color = 'blue', size=1) +
      ggtitle(paste0(iname, ',  R2 = ' , round(Rsquared[[iname]], digits=2), 
                     ',  MAPE= ',round(MAPE[[iname]], digits=1))) + 
      xlab('') + ylab(iname)
    ggsave(p, filename=paste0('Plot.Yhat.',sfx.save,'_',iname,'.png'), width = 6, height = 5)
  }
}else{
  d.yhat.cv<-data.frame(names.y=yhat.cv)
  
  # plot(d.test[,i], type='l', main=paste(dimnames(d.test)[[2]][i],round(Rsquared[[i]],3),sep=" "),  col='blue', lwd=2)
  # lines(d.yhat.cv[,i,],col='red',lwd = 2)
  
  # p <- ggplot(d.test,aes_string(x = 1:nrow(d.test),y = iname)) + geom_line(color='green',size=1) +
  #   geom_line(data = d.yhat.cv,aes_string(x = 1:nrow(d.yhat.cv),y = iname),color = 'blue', size=1) +
  #   ggtitle(paste0(iname, ',  R2 = ' , round(Rsquared[[iname]], digits=2), ' (no WL)'))+ xlab('')
  
  p <- ggplot(d.test,aes(x = time.test, y = y.test)) + geom_line(color='green',size=1) +
    geom_line(data = d.yhat.cv, aes(x = time.test, y = yhat.cv),color = 'blue', size=1) +
    ggtitle(paste0(names.y, ',  R2 = ' , round(Rsquared[[names.y]], digits=2),
                   ',  MAPE= ',round(MAPE[[names.y]], digits=1))) + 
    xlab('') + ylab(names.y)
  ggsave(p, filename=paste0('Plot.Yhat.',sfx.save,'_',iname,'_single.png'), width = 6, height = 5)
}

# This is the code to extract the coefficients
coef.cv.mfit<-coef(cv.mfit) 
d.coef.cv.mfit<- data.frame(as.matrix(coef.cv.mfit[[1]]))
# d.coef.cv.mfit$COEF <- rownames(d.coef.cv.mfit)

for( i in 2:length(coef.cv.mfit)){
  d.coef.cv.mfit<- data.frame(d.coef.cv.mfit , data.frame(as.matrix(coef.cv.mfit[[i]])))
}
#resetting colnames
colnames(d.coef.cv.mfit)<-dimnames(d.test)[[2]]
d.coef.cv.mfit <- rbind(as.data.frame(Rsquared),as.data.frame(MAPE),as.data.frame(RMSE),
                        d.coef.cv.mfit)
rownames(d.coef.cv.mfit)[1:3] <- c('R^2', 'MAPE', 'RMSE')
# writing the coef file
write.csv(d.coef.cv.mfit, paste0('regression.',toupper(meter),'.', sfx.save,'.csv'))