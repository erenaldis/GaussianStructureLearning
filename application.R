library(pcalg)
library(glasso)
library(mice)

data <- read.csv(file="hw3_data/data.txt", header=TRUE, sep=" ")
n <- nrow(data)
V <- colnames(data)

alphas <- seq(0.01, 0.1, 0.01)
scores <- c()
num_edges <- c()
for (alpha in alphas) {
  ## estimate CPDAG
  pc.fit <- pc(suffStat = list(C = cor(data), n = n),
               indepTest = gaussCItest, ## indep.test: partial correlations
               alpha=alpha, labels = V)
  sample_dag <- pdag2dag(pc.fit@graph)$graph
  score <- new("GaussL0penObsScore", data)
  scores <- c(scores, score$global.score(as(sample_dag, "GaussParDAG")))
}

optimal.alpha <- alphas[which.max(scores)]
pc.fit <- pc(suffStat = list(C = cor(data), n = n),
             indepTest = gaussCItest, ## indep.test: partial correlations
             alpha=optimal.alpha, labels = V)
if (require(Rgraphviz)) {
  ## show estimated CPDAG
  par(mfrow=c(1,2))
  plot(pc.fit, main = "Estimated CPDAG")
}
jpeg('BIC_PC.jpg')
plot(alphas, scores, type='l', main='BIC Score vs. Alpha')
dev.off()

nr <- 100
rhos <- seq(0.1, 1, length=nr)
penalized_log_lik <- c()
for (rho in rhos) {
  g <- glasso(var(data), rho=rho)
  p_off_d <- sum(g$wi!=0 & col(var(data))<row(var(data)))
  penalized_log_lik <- c(penalized_log_lik, -2*(g$loglik) + p_off_d*log(n))
}
jpeg('loglik_glasso.jpg')
plot(rhos, penalized_log_lik, type='l', main='BIC vs. Penalty Term for Glasso')
dev.off()


data_mis <- read.csv(file="hw3_data/data_missing.txt", header=TRUE, sep=" ")
imp <- mice(data_mis, m=10, printFlag = FALSE)

A <- matrix( 0, nrow = 20, ncol = 20)

for (m in 1:10) {
  data_complete <- complete(imp, m)
  n <- nrow(data_complete)
  V <- colnames(data_complete)
  
  alphas <- seq(0.01, 0.1, 0.01)
  scores <- c()
  num_edges <- c()
  
  
  ## estimate CPDAG
  pc.fit <- pc(suffStat = list(C = cor(data_complete), n = n),
               indepTest = gaussCItest, ## indep.test: partial correlations
               alpha=0.05, labels = V)
  A <- A + as(pc.fit, "amat")
}

A_mean <- round(A/10)

d<-graph_from_adjacency_matrix(A_mean)
jpeg('Averaged_CPDAG_PC.jpg')
plot(d, label=V)
dev.off()


P <- matrix( 0, nrow = 20, ncol = 20)

for (m in 1:10) {
  data_complete <- complete(imp, m)
  n <- nrow(data_complete)
  V <- colnames(data_complete)
  
  
  g <- glasso(var(data_complete), rho=0.38)
  P <- P + ifelse(g$wi!=0 & row(g$wi)!=col(g$wi),1,0)

}

P <- P/10

A <- ifelse(P > 0.5 & row(P)!=col(P),1,0)

g <- network(A)
jpeg('Lasso_Averaged_CPDAG.jpg')
ggnet2(g,label=V,main="Averaged CPDAG using Lasso")
dev.off()
