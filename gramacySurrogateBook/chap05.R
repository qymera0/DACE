
# 5.0 LOAD PACKAGE --------------------------------------------------------

library(plgp)
library(mvtnorm)
library(lhs)
eps <- sqrt(.Machine$double.eps)

# 5.1 GP Prior ------------------------------------------------------------

n <- 100

X <- matrix(seq(0, 10, length = n), ncol = 1)

# Euclidian distance between all points

D <- distance(X)

# Covariance Matrix

eps <- sqrt(.Machine$double.eps)

sigma <- exp(-D) + diag(eps, n)

Y <- rmvnorm(1, sigma = sigma)

plot(X, Y, type = "l")

Y <- rmvnorm(3, sigma = sigma)

matplot(X, t(Y), type = "l", ylab = "Y")

# Simple 1d GP prediction example

rm(list = ls())

n <- 8

X <- matrix(seq(0, 2*pi, length = n), ncol = 1)

y <- sin(X)

D <- distance(X)

eps <- sqrt(.Machine$double.eps)

sigma <- exp(-D) + diag(eps, ncol(D))

XX <- matrix(seq(-0.5, 2*pi + 0.5, length = 100), ncol = 1)

DXX <- distance(XX)

SXX <- exp(-DXX) + diag(eps, ncol(DXX))

DX <- distance(XX, X)

SX <- exp(-DX)

Si <- solve(sigma)

mup <- SX %*% Si %*% y

sigmap <- SXX - SX %*% Si %*% t(SX)

YY <- rmvnorm(100, mup, sigmap)

q1 <- mup + qnorm(0.05, 0, sqrt(diag(sigmap)))

q2 <- mup + qnorm(0.95, 0, sqrt(diag(sigmap)))

matplot(XX, t(YY), type="l", col="gray", lty=1, xlab="x", ylab="y")
points(X, y, pch=20, cex=2)
lines(XX, sin(XX), col="blue")
lines(XX, mup, lwd=2)
lines(XX, q1, lwd=2, lty=2, col=2)
lines(XX, q2, lwd=2, lty=2, col=2)

# Higher dimensions

nx <- 20

x <- seq(0, 2, length = nx)

X <- expand.grid(x, x)

D <- distance(X)

eps <- sqrt(.Machine$double.eps)

Sigma <- exp(-D) + diag(eps, nrow(X))

Y <- rmvnorm(2, sigma = Sigma)

par(mfrow=c(1,2)) 

persp(x, x, matrix(Y[1,], ncol=nx), theta=-30, phi=30, xlab="x1", 
      ylab="x2", zlab="y")

persp(x, x, matrix(Y[2,], ncol=nx), theta=-30, phi=30, xlab="x1", 
      ylab="x2", zlab="y")

# toy example

X <- randomLHS(40, 2)

X[,1] <- (X[,1] - 0.5)*6 + 1

X[,2] <- (X[,2] - 0.5)*6 + 1

y <- X[,1]*exp(-X[,1]^2 - X[,2]^2)

xx <- seq(-2, 4, length = 40)

XX <- expand.grid(xx, xx)

D <- distance(X)

Sigma <- exp(-D)

DXX <- distance(XX)

SXX <- exp(-DXX) + diag(eps, ncol(DXX))

DX <- distance(XX, X)

SX <- exp(-DX)

Si <- solve(Sigma)

mup <- SX %*% Si %*% y

Sigmap <- SXX - SX %*% Si %*% t(SX)

sdp <- sqrt(diag(Sigmap))

par(mfrow=c(1,2))

cols <- heat.colors(128)

image(xx, xx, matrix(mup, ncol=length(xx)), xlab="x1", ylab="x2", col=cols)

points(X[,1], X[,2])

image(xx, xx, matrix(sdp, ncol=length(xx)), xlab="x1", ylab="x2", col=cols)

points(X[,1], X[,2])

persp(xx, xx, matrix(mup, ncol=40), theta=-30, phi=30, xlab="x1", 
      ylab="x2", zlab="y")

# 5.2 GP hyperparameters --------------------------------------------------

# 5.2.1 Scale

par(mfrow = c(1, 1))

rm(list = ls())

eps <- sqrt(.Machine$double.eps)

n <- 100

X <- matrix(seq(0, 10, length = n), ncol = 1)

D <- distance(X)

C <- exp(-D) + diag(eps, n)

tau2 <- 25

Y <- rmvnorm(10, sigma = tau2*C)

matplot(X, t(Y), type="l")

n <- 8

X <- matrix(seq(0, 2*pi, length = n), ncol = 1)

y <- 5*sin(X)

D <- distance(X)

Sigma <- exp(-D)

XX <- matrix(seq(-0.5, 2*pi + 0.5, length=100), ncol=1)

DXX <- distance(XX)

SXX <- exp(-DXX) + diag(eps, ncol(DXX))

DX <- distance(XX, X)

SX <- exp(-DX)

Si <- solve(Sigma); 

mup <- SX %*% Si %*% y

Sigmap <- SXX - SX %*% Si %*% t(SX)

YY <- rmvnorm(100, mup, Sigmap)

q1 <- mup + qnorm(0.05, 0, sqrt(diag(Sigmap)))

q2 <- mup + qnorm(0.95, 0, sqrt(diag(Sigmap)))

matplot(XX, t(YY), type="l", col="gray", lty=1, xlab="x", ylab="y")

points(X, y, pch=20, cex=2)

lines(XX, mup, lwd=2)

lines(XX, 5*sin(XX), col="blue")

lines(XX, q1, lwd=2, lty=2, col=2)

lines(XX, q2, lwd=2, lty=2, col=2)

CX <- SX

Ci <- Si

CXX <- SXX

tau2hat <- drop(t(y) %*% Ci %*% y / length(y))

mup2 <- CX %*% Ci %*% y

Sigmap2 <- tau2hat*(CXX - CX %*% Ci %*% t(CX))

YY <- rmvnorm(100, mup2, Sigmap2)

q1 <- mup + qnorm(0.05, 0, sqrt(diag(Sigmap2)))

q2 <- mup + qnorm(0.95, 0, sqrt(diag(Sigmap2)))

matplot(XX, t(YY), type="l", col="gray", lty=1, xlab="x", ylab="y")
points(X, y, pch=20, cex=2)
lines(XX, mup, lwd=2)
lines(XX, 5*sin(XX), col="blue")
lines(XX, q1, lwd=2, lty=2, col=2); lines(XX, q2, lwd=2, lty=2, col=2)

score <- function(Y, mu, Sigma, mah=FALSE)
{
        Ymmu <- Y - mu
        Sigmai <- solve(Sigma)
        mahdist <- t(Ymmu) %*% Sigmai %*% Ymmu
        if(mah) return(sqrt(mahdist))
        return (- determinant(Sigma, logarithm=TRUE)$modulus - mahdist)
}

Ytrue <- 5*sin(XX)

df <- data.frame(score(Ytrue, mup, Sigmap, mah=TRUE), 
                 score(Ytrue, mup2, Sigmap2, mah=TRUE))

colnames(df) <- c("tau2=1", "tau2hat")

df

# 5.2.2 Noise and nuggets

nlg <- function(g, D, Y) 
{
        n <- length(Y)
        K <- exp(-D) + diag(g, n)
        Ki <- solve(K)
        ldetK <- determinant(K, logarithm=TRUE)$modulus
        ll <- - (n/2)*log(t(Y) %*% Ki %*% Y) - (1/2)*ldetK
        counter <<- counter + 1
        return(-ll)
}

eps <- sqrt(.Machine$double.eps)

n <- 100

X <- matrix(seq(0, 10, length=n), ncol=1)

X <- rbind(X, X)

n <- nrow(X)

y <- 5*sin(X) + rnorm(n, sd=1)

D <- distance(X)

counter <- 0

g <- optimize(nlg, interval=c(eps, var(y)), D=D, Y=y)$minimum

g

K <- exp(-D) + diag(g, n)

Ki <- solve(K)

tau2hat <- drop(t(y) %*% Ki %*% y / n)

c(tau=sqrt(tau2hat), sigma=sqrt(tau2hat*g))

DX <- distance(XX, X)

KX <- exp(-DX)

DXX <- distance(XX)

KXX <- exp(-DXX) + diag(g, nrow(DXX))

mup <- KX %*% Ki %*% y
Sigmap <- tau2hat*(KXX - KX %*% Ki %*% t(KX))
q1 <- mup + qnorm(0.05, 0, sqrt(diag(Sigmap)))
q2 <- mup + qnorm(0.95, 0, sqrt(diag(Sigmap)))

Sigma.int <- tau2hat*(exp(-DXX) + diag(eps, nrow(DXX)) 
                      - KX %*% Ki %*% t(KX))
YY <- rmvnorm(100, mup, Sigma.int)

matplot(XX, t(YY), type="l", lty=1, col="gray", xlab="x", ylab="y")
points(X, y, pch=20, cex=2)
lines(XX, mup, lwd=2)
lines(XX, 5*sin(XX), col="blue")
lines(XX, q1, lwd=2, lty=2, col=2)
lines(XX, q2, lwd=2, lty=2, col=2)

# 5.2.3 Derivative-based hyperparameter optimization ----------------------

gnlg <- function(g, D, Y)
{
        n <- length(Y)
        K <- exp(-D) + diag(g, n)
        Ki <- solve(K)
        KiY <- Ki %*% Y
        dll <- (n/2) * t(KiY) %*% KiY / (t(Y) %*% KiY) - (1/2)*sum(diag(Ki))
        return(-dll)
}

counter <- 0

out <- optim(0.1*var(y), nlg, gnlg, method="L-BFGS-B", lower=eps, 
             upper=var(y), D=D, Y=y)

c(g, out$par)


# # 5.2.4 Lengthscale: rate of decay of correlation -----------------------

nl <- function(par, D, Y) 
{
  theta <- par[1]                                       ## change 1
  g <- par[2]
  n <- length(Y)
  K <- exp(-D/theta) + diag(g, n)                       ## change 2
  Ki <- solve(K)
  ldetK <- determinant(K, logarithm=TRUE)$modulus
  ll <- - (n/2)*log(t(Y) %*% Ki %*% Y) - (1/2)*ldetK
  counter <<- counter + 1
  return(-ll)
}

X2 <- randomLHS(40, 2)

X2 <- rbind(X2, X2)

X2[,1] <- (X2[,1] - 0.5)*6 + 1

X2[,2] <- (X2[,2] - 0.5)*6 + 1

y2 <- X2[,1]*exp(-X2[,1]^2 - X2[,2]^2) + rnorm(nrow(X2), sd=0.01)

D <- distance(X2)

counter <- 0

out <- optim(c(0.1, 0.1*var(y2)), nl, method="L-BFGS-B", lower=eps, 
             upper=c(10, var(y2)), D=D, Y=y2) 

out$par

brute <- c(out$counts, actual=counter)

brute

gradnl <- function(par, D, Y)
{
  ## extract parameters
  theta <- par[1]
  g <- par[2]
  
  ## calculate covariance quantities from data and parameters
  n <- length(Y)
  K <- exp(-D/theta) + diag(g, n)
  Ki <- solve(K)
  dotK <- K*D/theta^2
  KiY <- Ki %*% Y
  
  ## theta component
  dlltheta <- (n/2) * t(KiY) %*% dotK %*% KiY / (t(Y) %*% KiY) - 
    (1/2)*sum(diag(Ki %*% dotK))
  
  ## g component
  dllg <- (n/2) * t(KiY) %*% KiY / (t(Y) %*% KiY) - (1/2)*sum(diag(Ki))
  
  ## combine the components into a gradient vector
  return(-c(dlltheta, dllg))
}

counter <- 0

outg <- optim(c(0.1, 0.1*var(y2)), nl, gradnl, method="L-BFGS-B", 
              lower=eps, upper=c(10, var(y2)), D=D, Y=y2) 

rbind(grad=outg$par, brute=out$par)

rbind(grad=c(outg$counts, actual=counter), brute)

K <- exp(- D/outg$par[1]) + diag(outg$par[2], nrow(X2))

Ki <- solve(K)

tau2hat <- drop(t(y2) %*% Ki %*% y2 / nrow(X2))

gn <- 40

xx <- seq(-2, 4, length=gn)

XX <- expand.grid(xx, xx)

DXX <- distance(XX)

KXX <- exp(-DXX/outg$par[1]) + diag(outg$par[2], ncol(DXX))

DX <- distance(XX, X2)

KX <- exp(-DX/outg$par[1])

mup <- KX %*% Ki %*% y2

Sigmap <- tau2hat*(KXX - KX %*% Ki %*% t(KX))

sdp <- sqrt(diag(Sigmap))

par(mfrow=c(1,2))

image(xx, xx, matrix(mup, ncol=gn), main="mean", xlab="x1",
      ylab="x2", col=cols)

points(X2)

image(xx, xx, matrix(sdp, ncol=gn), main="sd", xlab="x1",
      ylab="x2", col=cols)

points(X2)

# 5.2.5 Anisotropic modeling

fried <- function(n=50, m=6)
{
  if(m < 5) stop("must have at least 5 cols")
  X <- randomLHS(n, m)
  Ytrue <- 10*sin(pi*X[,1]*X[,2]) + 20*(X[,3] - 0.5)^2 + 10*X[,4] + 5*X[,5]
  Y <- Ytrue + rnorm(n, 0, 1)
  return(data.frame(X, Y, Ytrue))
}

m <- 7

n <- 200

nprime <- 1000

data <- fried(n + nprime, m)

X <- as.matrix(data[1:n,1:m])

y <- drop(data$Y[1:n])

XX <- as.matrix(data[(n + 1):(n + nprime),1:m])

yy <- drop(data$Y[(n + 1):(n + nprime)])

yytrue <- drop(data$Ytrue[(n + 1):(n + nprime)])

D <- distance(X)

out <- optim(c(0.1, 0.1*var(y)), nl, gradnl, method="L-BFGS-B", lower=eps, 
             upper=c(10, var(y)), D=D, Y=y)

out

K <- exp(- D/out$par[1]) + diag(out$par[2], nrow(D))

Ki <- solve(K)

tau2hat <- drop(t(y) %*% Ki %*% y / nrow(D))

DXX <- distance(XX)

KXX <- exp(-DXX/out$par[1]) + diag(out$par[2], ncol(DXX))

DX <- distance(XX, X)

KX <- exp(-DX/out$par[1])

mup <- KX %*% Ki %*% y

Sigmap <- tau2hat*(KXX - KX %*% Ki %*% t(KX))

rmse <- c(gpiso=sqrt(mean((yytrue - mup)^2)))

rmse

library(mda)

fit.mars <- mars(X, y)

p.mars <- predict(fit.mars, XX)

rmse <- c(rmse, mars=sqrt(mean((yytrue - p.mars)^2)))

rmse

nlsep <- function(par, X, Y) 
{
  theta <- par[1:ncol(X)]  
  g <- par[ncol(X)+1]
  n <- length(Y)
  K <- covar.sep(X, d=theta, g=g)
  Ki <- solve(K)
  ldetK <- determinant(K, logarithm=TRUE)$modulus
  ll <- - (n/2)*log(t(Y) %*% Ki %*% Y) - (1/2)*ldetK
  counter <<- counter + 1
  return(-ll)
}

tic <- proc.time()[3]

counter <- 0 

out <- optim(c(rep(0.1, ncol(X)), 0.1*var(y)), nlsep, method="L-BFGS-B", 
             X=X, Y=y, lower=eps, upper=c(rep(10, ncol(X)), var(y)))

toc <- proc.time()[3]

out$par

brute <- c(out$counts, actual=counter)

brute

gradnlsep <- function(par, X, Y)
{
  theta <- par[1:ncol(X)]
  g <- par[ncol(X)+1]
  n <- length(Y)
  K <- covar.sep(X, d=theta, g=g) 
  Ki <- solve(K)
  KiY <- Ki %*% Y
  
  ## loop over theta components
  dlltheta <- rep(NA, length(theta))
  for(k in 1:length(dlltheta)) {
    dotK <- K * distance(X[,k])/(theta[k]^2)
    dlltheta[k] <- (n/2) * t(KiY) %*% dotK %*% KiY / (t(Y) %*% KiY) - 
      (1/2)*sum(diag(Ki %*% dotK))
  }
  
  ## for g   
  dllg <- (n/2) * t(KiY) %*% KiY / (t(Y) %*% KiY) - (1/2)*sum(diag(Ki))
  
  return(-c(dlltheta, dllg))
}

tic <- proc.time()[3]

counter <- 0

outg <- optim(c(rep(0.1, ncol(X)), 0.1*var(y)), nlsep, gradnlsep, 
              method="L-BFGS-B", lower=eps, upper=c(rep(10, ncol(X)), var(y)), X=X, Y=y) 
toc <- proc.time()[3]


thetahat <- rbind(grad=outg$par, brute=out$par)

colnames(thetahat) <- c(paste0("d", 1:ncol(X)), "g")

thetahat

rbind(grad=c(outg$counts, actual=counter), brute)

toc - tic

K <- covar.sep(X, d=outg$par[1:ncol(X)], g=outg$par[ncol(X)+1])
Ki <- solve(K)
tau2hat <- drop(t(y) %*% Ki %*% y / nrow(X))
KXX <- covar.sep(XX, d=outg$par[1:ncol(X)], g=outg$par[ncol(X)+1]) 
KX <- covar.sep(XX, X, d=outg$par[1:ncol(X)], g=0)
mup2 <- KX %*% Ki %*% y
Sigmap2 <- tau2hat*(KXX - KX %*% Ki %*% t(KX))

rmse <- c(rmse, gpsep=sqrt(mean((yytrue - mup2)^2)))
rmse
