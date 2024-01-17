#### OLD: Current priors (Cole did this) ####

r = rlnorm(500, log(2.7),0.5)
h = rlnorm(500, log(10),0.5)
b = rlnorm(500, log(18),0.5)
c = rlnorm(500, log(0.02),0.5)
p = seq(0,25, length = 500)

plot(density(r), type = 'l',xlab="x",ylab="f(x)", main = 'parameters, sd = 0.5', col = 1, xlim = c(0,50))
#lines(density(r),col="red")
lines(density(h),col=2)
lines(density(b),col=3)
lines(density(c),col=4)
legend("topright",c("R", "O", 'B', 'C'),lty=1,col=1:4)

#lognormals with the 0 value
r = rlnorm(500, 2.7,0.5)
h = rlnorm(500, log(10),0)
b = rlnorm(500, log(18),0)
c = rlnorm(500, log(0.02),0)
grid = seq(0,25,.1)

plot(density(r), type = 'l',xlab="x",ylab="f(x)", main = 'parameters, sd = 0', col = 1)
lines(density(h),col=2)
lines(density(b),col=3)
lines(density(c),col=4)
legend("topright",c("R", "H", 'B', 'C'),lty=1,col=1:4)


legend("topright",c("R", "O", 'B', 'C'),lty=1,col=1:2)

p = seq(0,1, length = 100)
plot(p, dbeta(p, 2,2), ylab = 'density', type = 'l', col = 4)
lines(p, dbeta(p, 2,2), col = 'red')

#the two betas have the same dist

plot()

#### OLD: Fucking around ####

# Normal dist.
r = rnorm(1000,2.5,1) # 2.5
mean(r)
h = rnorm(1000,0.06,0.5) # 0.06
mean(h)
b = rnorm(1000,35,1) # 35
mean(b)
c = rnorm(1000,0.2,1) # 0.1
mean(c)
u = rbeta(1000,2,2) # 0.2
mean(u)
O = rbeta(1000,2,2) # 0.01
mean(O)
plot(density(r), type = 'l',xlab="x",ylab="f(x)", main = 'Prior Probability Density', col = 1, xlim = c(0,100), ylim = c(0,1))
lines(density(h),col=2)
lines(density(b),col=3)
lines(density(c),col=4)
lines(density(u),col=5)
lines(density(O),col=6)
legend("topright",c("r = 2.5", "h = 0.06", 'b = 35', 'c = 0.2',"u = 0.2","O = 0.008"),lty=1,col=1:6)

O = rbeta(100000,0.5,50) # 0.008
mean(O)
plot(density(O), type = 'l',xlab="x",ylab="Density", main = 'Prior Probability Density', col = 1, xlim = seq(0,1), ylim = c(0,1))

O = rbinom(1000,size=1,prob=0.2)
mean(O)
plot(density(O), type = 'l',xlab="x",ylab="Density", main = 'Prior Probability Density', col = 1, xlim = seq(0,1), ylim = c(0,1))


# Boop
r = rnorm(1000,2.5,1) # 2.5
mean(r)
h = rnorm(1000,0.075,0.5) # 0.06
mean(h)
b = rnorm(1000,35,1) # 35
mean(b)
c = rnorm(1000,0.3,1) # 0.1
mean(c)
u = rnorm(1000,0,0.5) # 0.2
mean(u)
O = rnorm(1000,0,0.5) # 0.01
mean(O)
plot(density(r), type = 'l',xlab="x",ylab="f(x)", main = 'Prior Probability Density', col = 1, xlim = c(0,100), ylim = c(0,1))
lines(density(h),col=2)
lines(density(b),col=3)
lines(density(c),col=4)
lines(density(u),col=5)
lines(density(O),col=6)
legend("topright",c("r = 2.5", "h = 0.06", 'b = 35', 'c = 0.2',"u = 0.2","O = 0.008"),lty=1,col=1:6)








h = rnorm(1000,0.075,0.5) # 0.06
mean(h)
c = rnorm(1000,0.3,0.5) # 0.1
mean(c)
u = rnorm(1000,0.4,0.5) # 0.2
mean(u)
O = rnorm(1000,0.015,0.5) # 0.01
mean(O)
plot(density(h), type = 'l',xlab="x",ylab="f(x)", main = 'Prior Probability Density', col = 1, xlim = c(0,1), ylim = c(0,1))
lines(density(c),col=2)
lines(density(u),col=3)
lines(density(O),col=4)
legend("topright",c("h = 0.075", 'c = 0.3',"u = 0.41","O = 0.012"),lty=1,col=1:4)

h = rnorm(1000,0.075,0.5) # 0.06
mean(h)
c = rnorm(1000,0.3,0.5) # 0.1
mean(c)
u = rnorm(1000,0.4,0.5) # 0.2
mean(u)
O = rnorm(1000,0.015,0.5) # 0.01
mean(O)
r = rnorm(1000,2.5,1) # 2.5
mean(r)
b = rnorm(1000,35,1) # 35
mean(b)
Z_init = rlnorm(1000,log(3),1)
plot(density(r), type = 'l',xlab="x",ylab="f(x)", main = 'Prior Probability Density', col = 1, xlim = c(0,50), ylim = c(0,1))
lines(density(O),col=2)
lines(density(h),col=3)
lines(density(b),col=4)
lines(density(c),col=5)
lines(density(u),col=6)
legend("topright",c("r = 2.5","O = 0.012","h = 0.075","b = 35","c = 0.3","u = 0.41"),lty=1,col=1:6)

Z_init = rlnorm(1000,log(3),1)
plot(density(Z_init), type = 'l',xlab="x",ylab="f(x)", main = 'Prior Probability Density', col = 1, xlim = c(0,400), ylim = c(0,1))
mean(Z_init)
log(80)
log(200)

sigma = rlnorm(1000,-1,1)
plot(density(sigma), type = 'l',xlab="x",ylab="f(x)", main = 'Prior Probability Density', col = 1, xlim = c(0,1), ylim = c(0,1))
mean(sigma)


h = rlnorm(1000,log(0.075),0.5) # 0.075
mean(h)
c = rlnorm(1000,log(0.3),0.5) # 0.3
mean(c)
u = rlnorm(1000,log(0.4),0.5) # 0.41
mean(u)
O = rlnorm(1000,log(0.015),0.5) # 0.01
mean(O)
r = rlnorm(1000,log(2.5),0.5) # 2.5
mean(r)
b = rlnorm(1000,log(35),0.5) # 35
mean(b)
plot(density(r), type = 'l',xlab="x",ylab="f(x)", main = 'Prior Probability Density', col = 1, xlim = c(0,50), ylim = c(0,1))
lines(density(O),col=2)
lines(density(h),col=3)
lines(density(b),col=4)
lines(density(c),col=5)
lines(density(u),col=6)
legend("topright",c("r = 2.5","O = 0.012","h = 0.075","b = 35","c = 0.3","u = 0.41"),lty=1,col=1:6)


x_log = runif(1000, log(10^(-6)), log(1)) # 10^(-6) = 0.000001
plot(density(x_log))
x = exp(x_log)
plot(density(x))
# > min(x)
# [1] 1.063668e-06
# > max(x)
# [1] 0.9885564

x <- rlnorm(100,3)
plot(density(x))

x <- rnorm(1000, 0.5, 1)
x <- sort(x); x <- rev(x); fn <- which(x < 0)[1]; x <- x[x[1]:fn]
plot(density(x))
mean(x); range(x)

x <- rnorm(1000, 1, 3)
plot(density(x))
mean(x); range(x)

x <- rlnorm(1000, 5, 1)
plot(density(x))
mean(x); range(x)

#### Antia et al. model =========================
## True values
r_true <- 0.2; k_true <- 0.01; p_true <- 1; o_true <- 1000

## r, greater than 0
r <- rnorm(1000, 1, 3)
r <- sort(r); r <- rev(r); fn <- which(r < 0)[1]; r <- r[r[1]:fn]
plot(density(r))
abline(v = r_true)
mean(r); range(r)

## k, between 0 and 1
k <- rlnorm(1000, log(0.1), 1)
k <- sort(k); g1 <- which(k >= 1)[1]; k <- k[1:g1-1]
#k <- rev(k); fn <- which(k < 0)[1]; k <- k[k[1]:fn]
plot(density(k))
abline(v = k_true)
mean(k); range(k)

## p, greater than 0
p <- rnorm(1000, 1, 1)
p <- sort(p); p <- rev(p); fn <- which(p < 0)[1]; p <- p[p[1]:fn]
plot(density(p))
abline(v = p_true)
mean(p); range(p)

## o, greater than 0
o <- rnorm(1000, 1000, 10)
#o <- sort(o); o <- rev(o); fn <- which(o < 0)[1]; o <- o[o[1]:fn]
plot(density(o))
abline(v = o_true)
mean(o); range(o)

#### Nowak and May model ========================
#### ... Damping oscillations ===================
## True values
r_true <-2.5; O_true <- 0.008; h_true <- 0.06; b_true <- 35; c_true <- 0.2; u_true <- 0.2

## r, greater than 0
r <- rnorm(1000, 1, 4)
r <- sort(r); r <- rev(r); fn <- which(r < 0)[1]; r <- r[r[1]:fn]
plot(density(r))
abline(v = r_true)
mean(r); range(r)

## O, between 0 and 1
O <- rnorm(1000, 0, 1)
O <- sort(O); g1 <- which(O >= 1)[1]; O <- O[1:g1-1]
O <- rev(O); fn <- which(O < 0)[1]; O <- O[O[1]:fn]
plot(density(O))
abline(v = O_true)
mean(O); range(O)

## h, greater than 0
h <- rnorm(1000, 0, 1)
h <- sort(h); h <- rev(h); fn <- which(h < 0)[1]; h <- h[h[1]:fn]
plot(density(h))
abline(v = h_true)
mean(h); range(h)

## b, greater than 0
b <- rnorm(1000, 100, 100)
b <- sort(b); b <- rev(b); fn <- which(b < 0)[1]; b <- b[b[1]:fn]
plot(density(b))
abline(v = b_true)
mean(b); range(b)

## c, greater than 0
c <- rnorm(1000, 0, 1)
c <- sort(c); c <- rev(c); fn <- which(c < 0)[1]; c <- c[c[1]:fn]
plot(density(c))
abline(v = c_true)
mean(c); range(c)

## u, between 0 and 1
u <- rnorm(1000, 0, 1)
u <- sort(u); g1 <- which(u >= 1)[1]; u <- u[1:g1-1]
u <- rev(u); fn <- which(u < 0)[1]; u <- u[u[1]:fn]
plot(density(u))
abline(v = u_true)
mean(u); range(u)

#### ... Stable limit cycle =====================
## True values
r_true <- 2.5; O_true <- 0.0012; h_true <- 0.075; b_true <- 35; c_true <- 0.3; u_true <- 0.41



