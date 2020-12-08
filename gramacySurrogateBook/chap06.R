library(lhs)
library(laGP)

# 6.2 SEQUENTIAL DESIGN ---------------------------------------------------

# Generate design

ninit <- 12

X <- randomLHS(ninit, 2)

f <- function(X, sd=0.01) 
{
        X[,1] <- (X[,1] - 0.5)*6 + 1
        X[,2] <- (X[,2] - 0.5)*6 + 1
        y <- X[,1] * exp(-X[,1]^2 - X[,2]^2) + rnorm(nrow(X), sd=sd)
}

y <- f(X)

# Initial model fitting

g <- garg(list(mle = TRUE, max = 1), y) # Priors GP correlation

d <- darg(list(mle = TRUE, max = 0.25), X) # Priors GP correlation

gpi <- newGP(X, y, d = d$start, g = g$start, dK = T)

mle <- jmleGP(gpi, c(d$min, d$max), c(g$min, g$max), d$ab, g$ab)

# Test grid

x1 <- x2 <- seq(0, 1, length = 100)

XX <- expand.grid(x1, x2)

yytrue <- f(XX, sd = 0)

initialPred <- predGP(gpi, XX, lite=T)

rmse <- sqrt(mean((yytrue - predGP(gpi, XX, lite=TRUE)$mean)^2))

# Determine where is the maximum variance

obj.alm <- function(x, gpi){
        # Return the diagonal predictive covariance matrix
        - sqrt(predGP(gpi, 
                      matrix(x, nrow = 1), 
                      lite = TRUE)$s2) # - because the answer will be minimized
}

mymaximin <- function(n, m, T=100000, Xorig=NULL) 
{   
        X <- matrix(runif(n*m), ncol=m)     ## initial design
        d <- distance(X)
        d <- d[upper.tri(d)]
        md <- min(d)
        if(!is.null(Xorig)) {               ## new code
                md2 <- min(distance(X, Xorig))
                if(md2 < md) md <- md2
        }
        
        for(t in 1:T) {
                row <- sample(1:n, 1)
                xold <- X[row,] ## random row selection
                X[row,] <- runif(m) ## random new row
                d <- distance(X)
                d <- d[upper.tri(d)]
                mdprime <- min(d)
                if(!is.null(Xorig)) {             ## new code
                        mdprime2 <- min(distance(X, Xorig))
                        if(mdprime2 < mdprime) mdprime <- mdprime2
                }
                if(mdprime > md) { md <- mdprime  ## accept
                } else { X[row,] <- xold }        ## reject
        }
        
        return(X)
}

xnp1.search <- function(X, gpi, obj=obj.alm, ...)
{
        start <- mymaximin(nrow(X), 2, T=100*nrow(X), Xorig=X)
        xnew <- matrix(NA, nrow=nrow(start), ncol=ncol(X) + 1)
        for(i in 1:nrow(start)) {
                out <- optim(start[i,], obj, method="L-BFGS-B", lower=0, 
                             upper=1, gpi=gpi, ...)
                xnew[i,] <- c(out$par, -out$value)
        }
        solns <- data.frame(cbind(start, xnew))
        names(solns) <- c("s1", "s2", "x1", "x2", "val")
        return(solns)
}

solns <- xnp1.search(X, gpi)

plot(X, xlab="x1", ylab="x2", xlim=c(0,1), ylim=c(0,1))

arrows(solns$s1, solns$s2, solns$x1, solns$x2, length=0.1)

m <- which.max(solns$val)

prog <- solns$val[m]

points(solns$x1[m], solns$x2[m], col=2, pch=20)

xnew <- as.matrix(solns[m, 3:4])

X <- rbind(X, xnew)

y <- c(y, f(xnew))

updateGP(gpi, xnew, y[length(y)])

mle <- rbind(mle, jmleGP(gpi, c(d$min, d$max), c(g$min, g$max), 
                         d$ab, g$ab))

rmse <- c(rmse, sqrt(mean((yytrue - predGP(gpi, XX, lite=TRUE)$mean)^2)))

solns <- xnp1.search(X, gpi)

m <- which.max(solns$val)

prog <- c(prog, solns$val[m])

xnew <- as.matrix(solns[m, 3:4])

X <- rbind(X, xnew)

y <- c(y, f(xnew))

updateGP(gpi, xnew, y[length(y)])

mle <- rbind(mle, jmleGP(gpi, c(d$min, d$max), c(g$min, g$max), 
                         d$ab, g$ab))
p <- predGP(gpi, XX, lite=TRUE)

rmse <- c(rmse, sqrt(mean((yytrue - p$mean)^2)))

plot(X, xlab="x1", ylab="x2", xlim=c(0,1), ylim=c(0,1))

arrows(solns$s1, solns$s2, solns$x1, solns$x2, length=0.1)

m <- which.max(solns$val)

points(solns$x1[m], solns$x2[m], col=2, pch=20)

# For Loop

for(i in nrow(X):24) {
        solns <- xnp1.search(X, gpi)
        m <- which.max(solns$val)
        prog <- c(prog, solns$val[m])
        xnew <- as.matrix(solns[m, 3:4])
        X <- rbind(X, xnew)
        y <- c(y, f(xnew))
        updateGP(gpi, xnew, y[length(y)])
        mle <- rbind(mle, jmleGP(gpi, c(d$min, d$max), c(g$min, g$max), 
                                 d$ab, g$ab))
        p <- predGP(gpi, XX, lite=TRUE)
        rmse <- c(rmse, sqrt(mean((yytrue - p$mean)^2)))
}

mle[seq(1,nrow(mle),by=2),]

plot((ninit+1):nrow(X), prog, xlab="n")

# Active learning Cohn

obj.alc <-
        function (x, gpi, Xref){
                - sqrt(alcGP(gpi, matrix(x, nrow = 1), Xref))
        }

deleteGP(gpi)

X <- X[1:ninit,]

y <- y[1:ninit]

g <- garg(list(mle=TRUE, max=1), y)

d <- darg(list(mle=TRUE, max=0.25), X)

gpi <- newGP(X, y, d=d$start, g=g$start, dK=TRUE)

mle <- jmleGP(gpi, c(d$min, d$max), c(g$min, g$max), d$ab, g$ab)

p <- predGP(gpi, XX, lite = TRUE)

rmse.alc <- sqrt(mean((yytrue - p$mean)^2))

Xref <- randomLHS(100, 2)

solns <- xnp1.search(X, gpi, obj=obj.alc, Xref=Xref)

m <- which.max(solns$val)

xnew <- as.matrix(solns[m, 3:4])

prog.alc <- solns$val[m]

plot(X, xlab="x1", ylab="x2", xlim=c(0,1), ylim=c(0,1))

arrows(solns$s1, solns$s2, solns$x1, solns$x2, length=0.1)

points(solns$x1[m], solns$x2[m], col=2, pch=20)

points(Xref, cex=0.25, pch=20, col="gray")

X <- rbind(X, xnew)

y <- c(y, f(xnew))

updateGP(gpi, xnew, y[length(y)])

mle <- rbind(mle, jmleGP(gpi, c(d$min, d$max), c(g$min, g$max), 
                         d$ab, g$ab))

p <- predGP(gpi, XX, lite=TRUE)

rmse.alc <- c(rmse.alc, sqrt(mean((yytrue - p$mean)^2)))

d <- darg(list(mle=TRUE), X)

for(i in nrow(X):99) {
        Xref <- randomLHS(100, 2)
        solns <- xnp1.search(X, gpi, obj=obj.alc, Xref=Xref)
        m <- which.max(solns$val)
        prog.alc <- c(prog.alc, solns$val[m])
        xnew <- as.matrix(solns[m, 3:4])
        X <- rbind(X, xnew)
        y <- c(y, f(xnew))
        updateGP(gpi, xnew, y[length(y)])
        mle <- rbind(mle, jmleGP(gpi, c(d$min, d$max), c(g$min, g$max), 
                                 d$ab, g$ab))
        p <- predGP(gpi, XX, lite=TRUE)
        rmse.alc <- c(rmse.alc, sqrt(mean((yytrue - p$mean)^2)))
}

par(mfrow=c(1,2))

plot((ninit+1):nrow(X), prog.alc, xlab="n: design size", 
     ylab="ALC progress")

plot(ninit:nrow(X), rmse, xlab="n: design size", ylab="OOS RMSE")

points(ninit:nrow(X), rmse.alc, col=2, pch=20)

legend("topright", c("alm", "alc"), pch=c(21,20), col=1:2)