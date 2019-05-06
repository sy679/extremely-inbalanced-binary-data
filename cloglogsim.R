install.packages("dplyr")
library(dplyr)

set.seed(1)
install.packages("evd")
library(evd)

install.packages("VGAM")
library(VGAM)

install.packages("mvtnorm")
library(mvtnorm)

install.packages("emdbook")
library(emdbook)

install.packages("dae")
library(dae)

library(MASS)

install.packages('truncnorm')
library(truncnorm)

install.packages("coda")
library(coda)

install.packages('splitstackshape')
library(splitstackshape)


set.seed(1)
###simulation

N=1000

#x1 <- sample(LETTERS[1:4], N, replace = T, prob = c(0.23, 0.54, 0.12, 0.11))
x2 <- sample(c("upper", "lower"), N, replace = T, prob = c(0.25, 0.75))
x3 <- rnorm(N, mean=0, sd= 0.4)

###simulation 
#sim.beta.true <- as.vector(c(-3, 0.1, 0.2, 0.5, 0.35, 0.25, 0.3))
sim.beta.true <- as.vector(c(-2, 0.6, 0.8))
sim.shape.true <- -0.7
#A1 <- ifelse(x1=="B", 1, 0)
#A2 <- ifelse(x1=="C", 1, 0)
#A3 <- ifelse(x1=="D", 1, 0)


#U1 <- ifelse(x2=="median", 1, 0)
U2 <- ifelse(x2 =="upper", 1, 0)

###test
#xpred <- as.matrix(cbind( A1, A2, A3, U1, U2, x3))
#sim.beta.true <- as.vector(c( 0.1, 0.2, 0.5, 0.35, 0.25, 0.3))
#xbeta <- tcrossprod(xpred, t(sim.beta.true))
#min(xbeta)


#xpred <- as.matrix(cbind(rep(1, N), A1, A2, A3, U1, U2, x3))
xpred <- as.matrix(cbind(rep(1, N),  U2, x3))

prob <- pgev(-tcrossprod(xpred, t(sim.beta.true)), shape = sim.shape.true, scale=1, lower.tail = F)
y <- rbinom(N, 1, prob = prob)
length(which(y==1))

model_log <- glm(y~xpred-1, family = binomial(link="logit"))
summary(model_log)

save(y, file="sim_y_gev.RDATA")
save(xpred, file="sim_x_gev.RDATA")

load("sim_y_gev.RDATA")
load("sim_x_gev.RDATA")


log.lik <- function(y, x, beta, xi){
  p<- pgev(-tcrossprod(x, t(beta)), shape = xi, scale = 1, lower.tail = F)
  p.up <- ifelse(p < 0.001, 0.001, p)
  p.up <- ifelse(p.up ==1, 0.999, p.up)
  loglik <- sum(dbinom(y, 1, prob = p.up, log = T))
  return (loglik)
}


iter <- 50000

est.simbeta <- matrix(NA, nrow =iter, ncol=3)
est.simshape <- matrix(NA, nrow = iter, ncol=1)
est.simdev <- matrix(NA, nrow=iter, ncol=1)

model_log <- glm(y~x2+x3, family = binomial(link="cloglog"))
summary(model_log)
sigma_beta <- diag(summary(model_log)$coefficients[, 2], nrow=3)

est.simbeta[1, ]<- c(-3, 0.5, 0.8)
#as.vector(model_log$coefficients)
est.simshape[1] <- -0.5
est.simdev[1] <- -2*log.lik(y,xpred, est.simbeta[1, ], est.simshape[1])


U.beta <- runif(iter)
U.shape <- runif(iter)

accept.simbeta <- 0
accept.simshape <- 0


pb <- utils::txtProgressBar(min=1, max=iter, style=3)
for (r in 2:iter){
  #print(r)
  shape.can <- rnorm(1, mean=est.simshape[r-1], sd=0.3)
  #a <- pnorm(0, mean = est.simshape[r-1], sd = 0.5)
  #b <- runif(1, 0, 1-a)
  #shape.can <- qnorm(a+b, mean = est.simshape[r-1], sd = 0.5)
  #shape.can <- rtnorm(1, mean = est.simshape[r-1], sd = 0.5, lower = 0, upper = Inf)
  #shape.can <- rnorm(1, mean = est.simshape[r-1], sd = 0.01)
  #shape.can <- rgamma(1, shape = est.simshape[r-1],  scale =2)
  
  g.fraction <- exp(log.lik(y, xpred, est.simbeta[r-1, ], shape.can) +
                      log(dnorm(shape.can, mean = 0, sd=0.3))
                    #+log(dtruncnorm(shape.can, a =0, b=Inf, mean = est.simshape[r-1], sd = 0.1))-
                    -log.lik(y, xpred, est.simbeta[r-1, ], est.simshape[r-1])-
                      log(dnorm(est.simshape[r-1], mean=0, sd=0.3))
                    #-log(dtruncnorm(est.simshape[r-1], a =0, b=Inf, mean =shape.can, sd = 0.1))
  )
  
  mh.g.fraction <- min(1, g.fraction)
  
  if(U.shape[r] < mh.g.fraction){
    est.simshape[r] <- shape.can
    accept.simshape <- accept.simshape+1
  }else{
    est.simshape[r] <- est.simshape[r-1]
  }
  
  #em.mbeta.best #rep(0, 7)
  beta.update <- function(index, beta.current, ...) {
    beta.star <- rnorm(n = 1, mean =beta.current[index], 
                       sd = 0.1)
    #beta.star1 <- runif(1, -3, 3) summary(model_log)$coefficients[, 1]
    if (index == 1)
      beta.temp <- c(beta.star, beta.current[2:3])
    else beta.temp <- beta.current ; beta.temp[index] <- beta.star
    r <- exp(log.lik(y, xpred, beta.temp, est.simshape[r]) +
               log(dmvnorm(beta.temp, mean=rep(0, 3),
                           sigma = diag(rep(1, 3),nrow=3)))-
               log.lik(y, xpred, est.simbeta[r-1, ], est.simshape[r])-
               log(dmvnorm(beta.current, mean =rep(0, 3), 
                           sigma =  diag(rep(1,3),nrow=3))))
    r <- min(r, 1, na.rm = T)
    if (runif(1) < r){
      beta.t <- beta.star
      #accept.simbeta[index] <- accept.simbeta[index]+1
    }else {beta.t<- beta.current[index]}
    return (beta.t)
  }
  
  est.simbeta[r, 1] <-beta.update(1, beta.current=est.simbeta[r-1, ])
  est.simbeta[r, 2] <-beta.update(2, beta.current=est.simbeta[r-1, ])
  est.simbeta[r, 3] <-beta.update(3, beta.current=est.simbeta[r-1, ])
  #est.simbeta[r, 4] <-beta.update(4, beta.current=est.simbeta[r-1, ])
  #est.simbeta[r, 5] <-beta.update(5, beta.current=est.simbeta[r-1, ])
  #est.simbeta[r, 6] <-beta.update(6, beta.current=est.simbeta[r-1, ])
  #est.simbeta[r, 7] <-beta.update(7, beta.current=est.simbeta[r-1, ])
  
  
  est.simdev[r] <- -2*log.lik(y, xpred, est.simbeta[r, ], est.simshape[r])
  utils::setTxtProgressBar(pb, r)
}
close(pb)

accept.simshape/iter
acfsample <- seq(500, iter, 80)

windows()
par(mfrow=c(2,2))
a1 <- quantile(est.simshape[acfsample ], 0.025)
a2 <- quantile(est.simshape[acfsample ], 0.975)
plot(density(est.simshape[acfsample]), main="shape parameter",xlab="",type="l")
abline(v=sim.shape.true, col="red")
abline(v=a1, col="blue")
abline(v=a2, col="blue")
for (i in 1:3){
  a1 <- quantile(est.simbeta[acfsample , i], 0.025)
  a2 <- quantile(est.simbeta[acfsample , i], 0.975)
  plot(density(est.simbeta[acfsample , i]),main=paste("Density for beta", i-1),xlab="",type="l")
  abline(v=sim.beta.true[i], col="red")
  abline(v=a1, col="blue")
  abline(v=a2, col="blue")
}

windows()
par(mfrow=c(2,2))
acf(est.simshape[acfsample], lag.max = 200, main="ACF for shape")
for(i in 1:3){
  acf(est.simbeta[acfsample, i], lag.max = 100,main=paste("ACF plot for beta", i-1))
}

length(unique(est.simshape[]))
length(unique(est.simbeta[, 2]))/iter

windows()
par(mfrow=c(2,2))
plot(est.simshape[acfsample], type="l", main="Trace plot of shape")
abline(h=sim.shape.true, col="red")

for(i in 1:3){
  plot(est.simbeta[acfsample, i], main=paste("trace plot for beta", i-1), type="l")
  abline(h=sim.beta.true[i], col="red")
}


effectiveSize(est.simshape[acfsample])
effectiveSize(est.simbeta[acfsample, 2])

burn = c(1:1000)

est.simbeta.mean <-colMeans(est.simbeta[acfsample, ])

devtheta.mean <- -2*log.lik(y, xpred, est.simbeta.mean, mean(est.simshape[acfsample]))

devhat <- mean(est.simdev[acfsample])

DIC.mean <- 2*devhat-devtheta.mean

DIC.mean

summary(model_log)

HPDinterval.mcmc <- function(obj, prob = 0.95, ...)
{
  obj <- as.matrix(obj)
  vals <- apply(obj, 2, sort)
  if (!is.matrix(vals)) stop("obj must have nsamp > 1")
  nsamp <- nrow(vals)
  npar <- ncol(vals)
  gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
  init <- 1:(nsamp - gap)
  inds <- apply(vals[init + gap, ,drop=FALSE] - vals[init, ,drop=FALSE],
                2, which.min)
  ans <- cbind(vals[cbind(inds, 1:npar)],
               vals[cbind(inds + gap, 1:npar)])
  dimnames(ans) <- list(colnames(obj), c("lower", "upper"))
  attr(ans, "Probability") <- gap/nsamp
  ans
}

HPDsimbeta <- HPDinterval.mcmc(est.simshape[acfsample ], prob=0.95)
write.csv(HPDsimbeta, file="hpdbeta.csv")

save(est.simshape, file="/home/syin/R/estsimshapegev0205.RDATA")
save(est.simbeta, file="/home/syin/R/estsimbetagev02052.RDATA")
save(est.simdev, file="/home/syin/R/estsimdevgev0205.RDATA")


load("estsimshapegev0205.RDATA")
load("estsimbetagev02052.RDATA")
load("estsimdevgev0205.RDATA")

#prob <- pweibull(-tcrossprod(xpred, t(sim.beta.true)), shape = sim.shape.true, scale=1, lower.tail = F)

prob_gev <- matrix(NA, nrow=nrow(xpred), ncol=10000)
for (i in 1:10000){
prob_pre <- pgev(-tcrossprod(xpred, t(est.simbeta.mean)), shape = mean(est.simshape[acfsample]),
                 scale=1, lower.tail = F)
prob_gev[, i] <- prob_pre
}
y_gev <- rbinom(N, 1, prob = prob_pre)
#colMeans(prob_gev)
#prob_w <- pweibull(-tcrossprod(xpred, t(est.simbeta.mean)), 
#                   shape = mean(est.simshape[acfsample]), scale=1, lower.tail = F)
#y_wfit <- rbinom(N, 1, prob = prob)
length(which(y_gev==1))

l1=0;l2=0;l3=0;l4=0
for (i in 1:N){
  if (y[i]==1 && y_gev[i]==1) l1 =l1+1
  if (y[i]==1 && y_gev[i]==0) l2 =l2+1
  if (y[i]==0 && y_gev[i]==1) l3 =l3+1
  if (y[i]==0 && y_gev[i]==0) l4 =l4+1
}
c(l1,l2, l3,l4)

