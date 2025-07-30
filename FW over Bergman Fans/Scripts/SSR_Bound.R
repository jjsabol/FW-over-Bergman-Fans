
q <- 8
Z <- matrix(rnorm(1000000*q, mean=0, sd=1), ncol=q)
Z_tr <- as.numeric(diff(rowMinsMaxs(Z)))
Vq <- var(Z_tr)
Eq <- mean(Z_tr)

w_min <- 1
n <- 25
sigma <- seq(from=0, to=w_min/(2*n*Eq), length.out=10000)
p_lb <- 1-(Vq/n)/(w_min/(2*n*sigma)-Eq)^2
kdata <- data.table(Sigma=sigma, Pr_LB = p_lb)
ggplot(kdata, aes(x=Sigma, y=Pr_LB))+
  geom_line() + geom_vline(xintercept=w_min/(2*n*Eq), lty=2, lwd=1, col='red') +
  ylim(0,1)