
# generate data
n <- 30
x1 <- c(1:n)
# time dependency
x1 <- x1/10 + 2
a <- 4
b <- 2

# regression model
y1 <- a + b *x1 + 0.1*rnorm(n)

# state: n*2 (intercept, slope)
x0 <- rep(1,n) # for intercept
xx <- cbind(x0,x1)
xx

# observations: n*1
y <- matrix(y1, nrow = n, ncol = 1)
y

# observation matrix: n*2
F <- matrix(xx, nrow=n, ncol=2)
F
# transition matrix: 2*2
G <- matrix(c(1,0,0,1), nrow = 2, ncol = 2)
G

# parameter variance: 2*2
W <- matrix(c(1,0,0,1), nrow = 2, ncol = 2)
W
# observation variance: 1*1
V <- matrix(1)
V

m0 <- matrix(c(5,1.5), nrow = 2, ncol = 1)
C0 <- matrix(c(.1,0,0,.1), nrow = 2, ncol = 2)

a<-0; R<-0; f<-0; Q<-0; e<-0; A <-0; m <-0; C<-0; tt<-0; 
m<-m0; C<-C0;
kfilt.m <- cbind(rep(0,n), rep(0,n))

#for (tt in 1:1) {
for (tt in 1:n) {  
  Fmat <- matrix(c(F[tt,1],F[tt,2]), nrow = 2, ncol = 1)
  print(Fmat)
  a <- G %*% m
  print(a)
  R <- G %*% C %*% t(G) + W
  print(R)
  f <- t(Fmat) %*% a
  print(f)
  Q <- t(Fmat) %*% R %*% Fmat + V
  print(Q)
  e <- y[tt]-f
  print(e)
  A <-R %*% Fmat %*% solve(Q)
  print(A)
  m <- a + A %*% e
  print(m)
  C <- R - A %*% Q %*% t(A)
  print(C)
  kfilt.m[tt,1] <- m[1,1]
  kfilt.m[tt,2] <- m[2,1]
  print(kfilt.m)
}
plot(kfilt.m[1:n,1])
plot(kfilt.m[1:n,2])
