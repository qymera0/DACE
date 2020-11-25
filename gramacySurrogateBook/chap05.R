
# 5.0 LOAD PACKAGE --------------------------------------------------------

library(plgp)
library(mvtnorm)
library(lhs)

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


