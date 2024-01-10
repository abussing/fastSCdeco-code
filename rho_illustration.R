

library(ggplot2)

n <- 1292

xvec <- rnbinom(n=n, mu=1.576, size=1/0.95)

factr <- rbinom(n=n, size=1, prob=0.59)


rho <- tanh(-0.42 - .196*xvec - 0.008*factr + 0.194*xvec*factr)

factr <- ifelse(factr==1, "level 1", "level 2")


ourdf <- data.frame(x1=xvec, rho=rho, factr=factr)


ggplot(data=ourdf, aes(x=x1, y=rho, group=factr)) + geom_point(aes(shape=factr)) +scale_shape_manual(values=c(16, 0)) + xlab(expression(x[1])) + ylab(expression(rho))


