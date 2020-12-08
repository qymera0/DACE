# 7.1 SURROGATE-ASSISTED OPTIMIZATION -------------------------------------


# Real function

f <- function(X)
{
        if(is.null(nrow(X))) X <- matrix(X, nrow=1)
        m <- 8.6928
        s <- 2.4269
        x1 <- 4*X[,1] - 2
        x2 <- 4*X[,2] - 2
        a <- 1 + (x1 + x2 + 1)^2 * 
                (19 - 14*x1 + 3*x1^2 - 14*x2 + 6*x1*x2 + 3*x2^2)
        b <- 30 + (2*x1 - 3*x2)^2 * 
                (18 - 32*x1 + 12*x1^2 + 48*x2 - 36*x1*x2 + 27*x2^2)
        f <- log(a*b)
        f <- (f - m)/s
        return(f)
}

# Initial design

library(lhs)

ninit <- 12

X <- randomLHS(ninit, 2)

y <- f(X)

# Initial GP

library(laGP)

da <- darg(list(mle = TRUE, max = 0.5), randomLHS(1000, 2))

gpi <- newGPsep(X, y, d = da$start, g = 1e-6, dK = T)

mle <- mleGPsep(gpi, param="d", tmin=da$min, tmax=da$max, ab=da$ab)$msg

# Optimization objective

obj.mean <-
        function(x, gpi){
                predGPsep(gpi, matrix(x, nrow = 1), lite = T)$mean
        }

m <- which.min(y)

opt <- optim(X[m, ], obj.mean, lower = 0, upper = 1, method="L-BFGS-B", gpi=gpi)

opt$par

plot(X[1:ninit,], xlab="x1", ylab="x2", xlim=c(0,1), ylim=c(0,1))

arrows(X[m,1], X[m,2], opt$par[1], opt$par[2], length=0.1)

# New evaluation

ynew <- f(opt$par)

updateGPsep(gpi, matrix(opt$par, nrow=1), ynew)

mle <- mleGPsep(gpi, param="d", tmin=da$min, tmax=da$max, ab=da$ab)

X <- rbind(X, opt$par)

y <- c(y, ynew)

# Solve for the new point

m <- which.min(y)

opt <- optim(X[m,], obj.mean, lower=0, upper=1, method="L-BFGS-B", gpi=gpi)

opt$par

plot(X, xlab="x1", ylab="x2", xlim=c(0,1), ylim=c(0,1))

arrows(X[m,1], X[m,2], opt$par[1], opt$par[2], length=0.1)

# Wrap the code on While loop

while(1) {
        
        m <- which.min(y)
        
        opt <- optim(X[m,], obj.mean, lower=0, upper=1, 
                     method="L-BFGS-B", gpi=gpi)
        
        ynew <- f(opt$par)
        
        if(abs(ynew - y[length(y)]) < 1e-4) break
        
        updateGPsep(gpi, matrix(opt$par, nrow=1), ynew)
        
        mle <- mleGPsep(gpi, param="d", tmin=da$min, tmax=da$max, ab=da$ab)
        
        X <- rbind(X, opt$par)
        
        y <- c(y, ynew)
}

deleteGPsep(gpi)

bov <- function(y, end=length(y))
{
        prog <- rep(min(y), end)
        prog[1:min(end, length(y))] <- y[1:min(end, length(y))]
        for(i in 2:end) 
                if(is.na(prog[i]) || prog[i] > prog[i-1]) prog[i] <- prog[i-1]
        return(prog)
}

prog <- bov(y)

plot(prog, type="l", col="gray", xlab="n: blackbox evaluations", 
     ylab="best objective value")

abline(v=ninit, lty=2)

legend("topright", "seed LHS", lty=2, bty="n")


# 7.2 EXPECTED IMPROVEMENT ------------------------------------------------

# Classic illustration

x <- c(1, 2, 3, 4, 12)

y <- c(0, -1.75, -2, -0.5, 5)

gpi <- newGP(matrix(x, ncol=1), y, d=10, g=1e-8)

xx <- seq(0, 13, length=1000)

p <- predGP(gpi, matrix(xx, ncol=1), lite=TRUE)

# EI calculation

m <- which.min(y)

fmin <- y[m]

d <- fmin - p$mean

s <- sqrt(p$s2)

dn <- d/s

ei <- d*pnorm(dn) + s*dnorm(dn)

par(mfrow=c(1,2))

plot(x, y, pch=19, xlim=c(0,13), ylim=c(-4,9), main="predictive surface")

lines(xx, p$mean)

lines(xx, p$mean + 2*sqrt(p$s2), col=2, lty=2)

lines(xx, p$mean - 2*sqrt(p$s2), col=2, lty=2)

abline(h=fmin, col=3, lty=3)

legend("topleft", c("mean", "95% PI", "fmin"), lty=1:3, 
       col=1:3, bty="n")

plot(xx, ei, type="l", col="blue", main="EI", xlab="x", ylim=c(0,0.15))

mm <- which.max(ei)

x <- c(x, xx[mm])

y <- c(y, p$mean[mm])

# Re evaluation

updateGP(gpi, matrix(xx[mm], ncol=1), p$mean[mm])

p <- predGP(gpi, matrix(xx, ncol=1), lite=TRUE)

deleteGP(gpi)

m <- which.min(y)

fmin <- y[m]

d <- fmin - p$mean

s <- sqrt(p$s2)

dn <- d/s

ei <- d*pnorm(dn) + s*dnorm(dn)

par(mfrow=c(1,2))

plot(x, y, pch=19, xlim=c(0,13), ylim=c(-4,9), main="predictive surface")

lines(xx, p$mean)

lines(xx, p$mean + 2*sqrt(p$s2), col=2, lty=2)

lines(xx, p$mean - 2*sqrt(p$s2), col=2, lty=2)

abline(h=fmin, col=3, lty=3)

legend("topleft", c("mean", "95% PI", "fmin"), lty=1:3, 
       col=1:3, bty="n")

plot(xx, ei, type="l", col="blue", main="EI", xlab="x", ylim=c(0,0.15))

# EI on our running example

EI <- function(gpi, x, fmin, pred=predGPsep)
{
        if(is.null(nrow(x))) x <- matrix(x, nrow=1)
        p <- pred(gpi, x, lite=TRUE)
        d <- fmin - p$mean
        sigma <- sqrt(p$s2)
        dn <- d/sigma
        ei <- d*pnorm(dn) + sigma*dnorm(dn)
        return(ei)
}

obj.EI <- function(x, fmin, gpi, pred=predGPsep){
        - EI(gpi, x, fmin, pred)
} 

eps <- sqrt(.Machine$double.eps) ## used lots below

EI.search <- function(X, y, gpi, pred=predGPsep, multi.start=5, tol=eps)
{
        m <- which.min(y)
        fmin <- y[m]
        start <- matrix(X[m,], nrow=1)
        if(multi.start > 1) 
                start <- rbind(start, randomLHS(multi.start - 1, ncol(X)))
        xnew <- matrix(NA, nrow=nrow(start), ncol=ncol(X)+1)
        for(i in 1:nrow(start)) {
                if(EI(gpi, start[i,], fmin) <= tol) { out <- list(value=-Inf); next }
                out <- optim(start[i,], obj.EI, method="L-BFGS-B", 
                             lower=0, upper=1, gpi=gpi, pred=pred, fmin=fmin)
                xnew[i,] <- c(out$par, -out$value)
        }
        solns <- data.frame(cbind(start, xnew))
        names(solns) <- c("s1", "s2", "x1", "x2", "val")
        solns <- solns[solns$val > tol,]
        return(solns)
}      

X <- randomLHS(ninit, 2)

y <- f(X)

gpi <- newGPsep(X, y, d=0.1, g=1e-6, dK=TRUE)

da <- darg(list(mle=TRUE, max=0.5), randomLHS(1000, 2))

solns <- EI.search(X, y, gpi)

m <- which.max(solns$val)

maxei <- solns$val[m]

plot(X, xlab="x1", ylab="x2", xlim=c(0,1), ylim=c(0,1))

arrows(solns$s1, solns$s2, solns$x1, solns$x2, length=0.1)

points(solns$x1[m], solns$x2[m], col=2, pch=20)

xnew <- as.matrix(solns[m,3:4])

X <- rbind(X, xnew)

y <- c(y, f(xnew))

updateGPsep(gpi, xnew, y[length(y)])

mle <- mleGPsep(gpi, param="d", tmin=da$min, tmax=da$max, ab=da$ab)

solns <- EI.search(X, y, gpi)

m <- which.max(solns$val)

maxei <- c(maxei, solns$val[m])

xnew <- as.matrix(solns[m,3:4])

X <- rbind(X, xnew)

y <- c(y, f(xnew))

updateGPsep(gpi, xnew, y[length(y)])

mle <- mleGPsep(gpi, param="d", tmin=da$min, tmax=da$max, ab=da$ab)

plot(X, xlab="x1", ylab="x2", xlim=c(0,1), ylim=c(0,1))

arrows(solns$s1, solns$s2, solns$x1, solns$x2, length=0.1)

points(solns$x1[m], solns$x2[m], col=2, pch=20)

end <- 50

for(i in nrow(X):end) {
        solns <- EI.search(X, y, gpi)
        m <- which.max(solns$val)
        maxei <- c(maxei, solns$val[m])
        xnew <- as.matrix(solns[m,3:4])
        ynew <- f(xnew)
        X <- rbind(X, xnew)
        y <- c(y, ynew)
        updateGPsep(gpi, xnew, y[length(y)])
        mle <- mleGPsep(gpi, param="d", tmin=da$min, tmax=da$max, ab=da$ab)
}

deleteGPsep(gpi)

optim.EI <- function(f, ninit, end)
{
        ## initialization
        X <- randomLHS(ninit, 2)
        y <- f(X)
        gpi <- newGPsep(X, y, d=0.1, g=1e-6, dK=TRUE)
        da <- darg(list(mle=TRUE, max=0.5), randomLHS(1000, 2))
        mleGPsep(gpi, param="d", tmin=da$min, tmax=da$max, ab=da$ab)
        
        ## optimization loop of sequential acquisitions
        maxei <- c()
        for(i in (ninit+1):end) {
                solns <- EI.search(X, y, gpi)
                m <- which.max(solns$val)
                maxei <- c(maxei, solns$val[m])
                xnew <- as.matrix(solns[m,3:4])
                ynew <- f(xnew)
                updateGPsep(gpi, xnew, ynew)
                mleGPsep(gpi, param="d", tmin=da$min, tmax=da$max, ab=da$ab)
                X <- rbind(X, xnew)
                y <- c(y, ynew)
        }
        
        ## clean up and return
        deleteGPsep(gpi)
        return(list(X=X, y=y, maxei=maxei))
}

reps <- 100

prog.ei <- matrix(NA, nrow=reps, ncol=end)

for(r in 1:reps) {
        os <- optim.EI(f, ninit, end)
        prog.ei[r,] <- bov(os$y)
}

# Conditional improvement