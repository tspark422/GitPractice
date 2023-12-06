### Q2
### (a)
x <- -2:2
X <- cbind(rep(1, 5), x)
W <- cbind(x^2, x^3)
Z <- cbind(X, W)

beta_hat <- solve(t(X)%*%X)%*%t(X)
t(c(1, 0))%*%solve(t(X)%*%X)%*%t(X)%*%W
t(c(0, 1))%*%beta_hat%*%W
t(c(1, -2))%*%solve(t(X)%*%X)%*%t(X)%*%W # theta1
t(c(1, -1))%*%solve(t(X)%*%X)%*%t(X)%*%W # theta2
t(c(1, 1))%*%solve(t(X)%*%X)%*%t(X)%*%W  # theta4
t(c(1, 2))%*%solve(t(X)%*%X)%*%t(X)%*%W  # theta5
t(W)%*%(diag(5)-X%*%solve(t(X)%*%X)%*%t(X))%*%W

### (b)
t(c(1, -2))%*%solve(t(X)%*%X)%*%c(1, -2)
d <- t(t(c(1, -2))%*%solve(t(X)%*%X)%*%t(X))
d0 <- Z%*%solve(t(Z)%*%Z)%*%c(1, -2, 4, -8)
d0perp <- d0-d
t(d0perp)%*%W

t(c(1, -1))%*%solve(t(X)%*%X)%*%c(1, -1) # MSE(theta2)
d <- t(t(c(1, -1))%*%solve(t(X)%*%X)%*%t(X))
d0 <- Z%*%solve(t(Z)%*%Z)%*%c(1, -1, 1, -1)
d0perp <- d0-d
t(d0perp)%*%W

t(c(1, 0))%*%solve(t(X)%*%X)%*%c(1, 0) # MSE(theta2)
d <- t(t(c(1, 0))%*%solve(t(X)%*%X)%*%t(X))
d0 <- Z%*%solve(t(Z)%*%Z)%*%c(1, 0, 0, 0)
d0perp <- d0-d
t(d0perp)%*%W

t(c(1, 1))%*%solve(t(X)%*%X)%*%c(1, 1) # MSE(theta2)
d <- t(t(c(1, 1))%*%solve(t(X)%*%X)%*%t(X))
d0 <- Z%*%solve(t(Z)%*%Z)%*%c(1, 1, 1, 1)
d0perp <- d0-d
t(d0perp)%*%W

t(c(1, 2))%*%solve(t(X)%*%X)%*%c(1, 2) # MSE(theta2)
d <- t(t(c(1, 2))%*%solve(t(X)%*%X)%*%t(X))
d0 <- Z%*%solve(t(Z)%*%Z)%*%c(1, 2, 4, 8)
d0perp <- d0-d
t(d0perp)%*%W
