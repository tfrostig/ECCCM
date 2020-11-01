library(MASS)


corr2 <- function(fam.size , rho, v.d = 1) {
  corr2 <- matrix(0 , ncol = fam.size , nrow = fam.size)
  for (i in 1:(fam.size - 1)) {
    corr2[i,((i + 1):fam.size)] <- rho ^ (1:(fam.size - i))
  }
  diag(corr2) <- v.d
  corr2[lower.tri(corr2)] <- t(corr2)[lower.tri(corr2)]
  return(corr2)
}

#### Create SNPs transforms normal vector into 0,1,2 while keeping MAF
snpMaker <- function(x.vec, maf) {
  g.vec <- x.vec
  g.vec[x.vec <= qnorm(1 - maf)] <- 0
  g.vec[x.vec >  qnorm(1 - maf)]  <- 1
  g.vec[x.vec >  qnorm(1 - (1 / 3) * maf)]  <- 2
  return(g.vec)
}

## Creating Block Covariance Matrix
## Input: Family size - dimension required, correlation = correlation inside each block, Var - var on diagonal
## output: required matrix

corrBlock <- function(famSize , blockSize , correlation = 0.8 , var = 1) {
  if ((famSize < blockSize) | (((famSize / blockSize) %% 1) != 0)) {
    print('Block Size Must Be Smaller Than Familly Size Or The Devision Of Them Isnt a Whole Number')
  }
  corrBlock <-  matrix(0,ncol = famSize , nrow = famSize)
  for (i in 1:(famSize/blockSize)) {
    startPos <- 1 + (i - 1) * blockSize
    corrBlock[startPos:(startPos + blockSize -1),startPos:(startPos + blockSize - 1)] <- correlation
  }
  diag(corrBlock) <- var
  return(corrBlock)
}



#### Functions for simulation
summarisePval <- function(pval.vec, qu, false.null.ind) {
  hypo.num          <- length(pval.vec)
  hypo.rejected.ind <- which(pval.vec < qu)
  fdr.power   <- c('FDR'   = calcFDR(hypo.rejected.ind, false.null.ind, hypo.num),
                   'power' = calcPower(hypo.rejected.ind, false.null.ind))
  return(fdr.power)
}


calcFDR <- function(hyp.rej, true.rej, hypo.num){
  if (length(hyp.rej) == 0) {
    return(0)
  } else {
    return(length(intersect(hyp.rej, setdiff(1:hypo.num, true.rej))) / length(hyp.rej))
  }
}

calcPower <- function(hyp.rej, true.rej) {
  if (length(true.rej) == 0) {
    return(0)
  } else {
    return(length(intersect(true.rej, hyp.rej)) / length(true.rej))
  }
}

findT <- function(beta.coef, ste, method = 'BH') {
  t.stat    <- beta.coef / ste
  p.vec     <- 2 * (1 - pnorm(abs(t.stat)))
  return(cbind('stat' = t.stat,
               'pval' = p.vec,
               'p.adjust' = p.adjust(p.vec, method)))
}



TEST <- function() {
  eff.dim <<- c(1,3)
  n.org <<- 1000
  n.ref <<- 1000
  h     <<- 0.05
  qu    <<- 0.05
  cov.mat  <<- corr2(15, 0.92)
  maf.vec  <<- runif(15)
  is.scale <<- TRUE
  k.tag    <<- 8
}


### Create Genes from a multivariate normal distirubtion
### Ensures estimated MAF > 0
createGene <- function(n, cov.mat, maf.vec, center = TRUE, scale = FALSE) {
  print(paste('Simulating, scaled SNPs, center =', center, 'scale =', scale))
  ### Estiamted MAF
  est.maf.vec      <- 0
  ### Simulating
  p       <- ncol(cov.mat)
  X <- mvrnorm(n, rep(0, p), cov.mat)
  while(any(est.maf.vec == 0) | any(is.nan(est.maf.vec))) {
    #### Sampling population
    X   <- mvrnorm(n, rep(0, p), cov.mat)
    #### Transforming them into Genotypes
    G   <- sweep(X, 2, maf.vec, snpMaker)
    #### Estimating MAF
    est.maf.vec      <- colMeans(G)
    #### Seperating into refrerence data and real data
    G.scale          <- scale(G, center = center, scale = scale)
  }
  return(G.scale)
}

TEST()

cas.dim <- length(eff.dim)
p       <- ncol(cov.mat)
#### Creating model
beta.vec    <- rep(0, p)
effect.size <- rep(1, p) ## Constant effect size
beta.vec[eff.dim] <- effect.size[eff.dim]
x.org <- createGene(n.org, cov.mat, maf.vec, scale = TRUE)
x.ref <- createGene(n.ref, cov.mat, maf.vec, scale = FALSE)
#### Creating the model
if (length(eff.dim) > 0) {
  true.y           <- drop(rbind(x.org, scale(x.ref)) %*% beta.vec)
  #unexplain.var   <- ((1 - h) / h) * beta.vec %*% cov.mat %*% beta.vec
  unexplain.var    <- (1 - h) / h
  sigma            <- sqrt(unexplain.var)
  y                <- true.y + rnorm(n.org + n.ref, mean = 0, sd = sqrt(unexplain.var))
}
if (length(eff.dim) == 0) {
  sigma           <- 1
  y <- 0 + rnorm(n.org + n.ref, mean = 0, sd = sigma)
}
### Scaling y
sd.y.org <- sd(y)
y.org    <- scale(y[1:n.org])
sigma    <- sigma / sd.y.org
#### The full estimator
inv.o        <- solve((t(x.org) %*% x.org) / (n.org - 1))
beta.org     <- solve(t(x.org) %*% x.org) %*% t(x.org) %*% y.org
var.beta.org <- sigma^2 * diag(inv.o) / n.org
test.oracle  <- testCoef(drop(beta.org),
                         drop(var.beta.org),
                         method = 'BH')
beta.report  <- drop(solve(diag(apply(x.org, 2, var))) %*% t(x.org) %*% y.org / n.org)
##### Naive estimation
p   <- ncol(x.org)
n.r <- nrow(x.org)
cov.mat        <- cov(x.ref)
cov.mat.inv    <- solve(cov.mat)
cov.list.r     <- list('cov'   = cov.mat,
                       'omega' = cov.mat.inv)
marg.to.joint <- ECCCM::marginalToJoint(marg.beta.hat = beta.report,
                                        n.o = n.org,
                                        cor.r = cov.list.r$cov,
                                        sigma = 1)

test.naive.df <- ECCCM::testCoef(est.beta = marg.to.joint$est.beta.hat,
                                 var.beta = marg.to.joint$naive.var.beta.hat,
                                 method   = 'BH')
analyzeRef(beta.report, x.r = x.ref, n.o = n.org)

