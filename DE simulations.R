## This file implements simulations to find the distribution of the (scaled) Maximally selected CUSUM for data X_i, which are L^2[0,1] valued random functions. 

## This code was written by me, and subsequently optimized by Claude to significantly improve run-time. 



## =====================================================================

##  Functional CUSUM  ->  Gumbel limit:  Monte-Carlo study  (fast version)

##

##  Same statistic, same RNG stream, ~2 orders of magnitude faster.

##  Key changes are flagged with  ##<<  in the comments.

## =====================================================================

## ---------------------------------------------------------------------

## 1.  Model ingredients

## ---------------------------------------------------------------------

get.eigenvalues <- function(m = 1, n = 100, alpha = 2) { # This generates eigenvalues which decay less slowly.
  
  if (m > n) {
    
    warning("Multiplicity higher than n")
    
    rep(1, n)
    
  } else if (m == n) {
    
    rep(1, m)
    
  } else {
    
    c(rep(1, m), 1 / (2:(n-m+1))^alpha)
    
  }
  
}

##<< The basis is deterministic, so evaluate it on the grid ONCE.

##<< Rows = basis index j, cols = grid point t.   e_j(t) = sqrt(2) sin((j-1/2) pi t)


library("pracma")

# get.basis.matrix <- function(N.order, t.grid) {
#   
#   sqrt(2) * sin(outer((seq_len(N.order) - 0.5) * pi, t.grid))
#   
# }

get.basis.matrix.orthonormal <- function(N.order, t.grid) {
  # requires that the length of the t-grid is larger than N.order
    
  BM.raw <- sqrt(2) * sin(outer((seq_len(N.order) - 0.5) * pi, t.grid))
  
  BM.raw.GS <- gramSchmidt(t(BM.raw))
  
  t(sqrt(length(t.grid)) * BM.raw.GS$Q)
}

##<< Fold sqrt(lambda_j) into the basis: W[j, ] = sqrt(lambda_j) e_j(t).

##<< Then a whole sample is one matrix product:  X = Z %*% W.

get.weight.matrix <- function(sqrt.eigenvalues, basis.matrix) {
  
  sqrt.eigenvalues * basis.matrix        # recycles down rows (column-major)
  
}

a.fn <- function(u) sqrt(2 * log(u))

# get.k.sigma <- function(alpha = 2, n.order = 100) {
#     if(n.order < 2){
#       return(1)
#     }
#   
#     lambdas <- (2:n.order)^(-alpha)
#     
#     sqrt(prod(1 / (1 - lambdas)))
#     
# }

get.k.sigma.fixed <- function(alpha = 2, n.order = 100, m = 1) {
  if(n.order <= m){
    return(1)
  }
  
  lambdas <- (2:(n.order-m+1))^(-alpha)
  
  sqrt(prod(1 / (1 - lambdas)))
  
}
  

make.b.fn <- function(m, k.sigma) {
  
  function(u) {
    
    2 * log(u) + m/2 * log(log(u)) -
      
      lgamma(m/2) + log(k.sigma)
    
  }
  
}

## ---------------------------------------------------------------------

## 2.  Core simulation

## ---------------------------------------------------------------------

##<< Column-wise cumsum: 50 C-level calls instead of apply()'s split/relist.

col.cumsum <- function(M) {
  
  S <- M
  
  for (k in seq_len(ncol(M))) S[, k] <- cumsum(M[, k])
  
  S
  
}

## One Monte-Carlo replicate: returns T_N and H_N for every N in N.seq.

## The samples are nested -- the first N rows are reused for each N,

## exactly as in the original code.

get.H.N.seq.fast <- function(N.seq, W, random.fn, a.seq, b.seq) {
  
  N.final  <- max(N.seq)
  
  N.order  <- nrow(W)
  
  t.length <- ncol(W)
  
  ##<< One RNG call + one BLAS call replace the N.final x N.order R loop.
  
  ##<< byrow = TRUE keeps the draw order identical to the original code.
  
  Z <- matrix(random.fn(N.final * N.order),
              
              nrow = N.final, ncol = N.order, byrow = TRUE)
  
  X <- Z %*% W                                   # N.final x t.length
  
  S <- col.cumsum(X)                             # partial sums S_k
  
  ##<< ||S_k||^2 does not depend on N, so compute it once for all of N.seq.
  
  ss <- rowSums(S * S)
  
  T.N <- numeric(length(N.seq))
  
  for (u in seq_along(N.seq)) {
    
    N  <- N.seq[u]
    
    k  <- seq_len(N - 1L)                        ##<< integer index, no floor()
    
    tk <- k / N
    
    ## B_N(t_k) = (S_k - t_k S_N) / N ; expand the norm to avoid building B.
    
    ##   ||S_k - t_k S_N||^2 = ||S_k||^2 - 2 t_k <S_k, S_N> + t_k^2 ||S_N||^2
    
    cross <- as.vector(S[k, , drop = FALSE] %*% S[N, ])
    
    num   <- ss[k] - 2 * tk * cross + tk * tk * ss[N]
    
    ## T_N(t)^2 = N ||B_N(t)||^2 / (t (1-t)),  ||.||^2 = (1/t.length) sum over grid
    
    T.N.sq <- num / (N * t.length * tk * (1 - tk))
    
    T.N[u] <- sqrt(max(T.N.sq))                  ##<< one sqrt, not N-1 of them
    
  }
  
  rbind(T = T.N, H = a.seq * T.N - b.seq)
  
}

## ---------------------------------------------------------------------

## 3.  Driver

## ---------------------------------------------------------------------


## ---------------------------------------------------------------------

## 4.  Comparison with the limit law

## ---------------------------------------------------------------------

## Limit: Gumbel with location log 2, i.e. F(x) = exp(-2 exp(-x))

d.limiting <- function(x, loc = 0, scale = 1) {
  
  ex <- exp(-(x - loc) / scale)
  
  2 * exp(-2 * ex) * ex / scale
  
}

p.limiting <- function(x) exp(-2 * exp(-x))       ##<< drops the 'ordinal' dependency



get.random.fn <- function(name){
  if(name == "normal"){
    function(n) rnorm(n)
  } else if (name == "exponential") {
    function(n) rexp(n) - 1
  } else if (name == "t"){
    function(n, df = 5) rt(n, df = df)*sqrt((df-2)/df)
  } else if (name == "uniform"){
    function(n, min = -1, max = 1) runif(n, min = -1, max = 1)*sqrt(3)
  }
}



## Simulation settings

#getwd()
#setwd("C:\\Users\\rgajendr\\AppData\\Local\\OneDrive - University of Waterloo\\Research\\DE asymptotics simulations")
#setwd("/Users/Rohit/Library/CloudStorage/OneDrive-UniversityofWaterloo/Research/DE asymptotics simulations")

output <- c()

p.limiting <- function(x) exp(-2 * exp(-x))

N.runs  <- 1e4
N.seq <- c(100, 500, 1000, 5000, 10000, 100000)
t.grid  <- seq(0, 1, length.out = 101)


random.func <- c("normal","exponential","t","uniform")
alpha.values <- c(1.5,2)
d.values <- c(1,3,50)
multiplicity.values <- c(1,3)


hyper.param.df <- droplevels(expand.grid(d.values,random.func, alpha.values, multiplicity.values))
names(hyper.param.df) <- c("d","RV","alpha","m")

rows.to.exclude <- ((hyper.param.df$d == 1)&((hyper.param.df$m >1)|(hyper.param.df$alpha < 2)))|((hyper.param.df$d == 3)&(hyper.param.df$m == 3)&(hyper.param.df$alpha < 2))
hyper.param.df <- hyper.param.df[!rows.to.exclude,]


seed <- sample(1e6, size = 1)
set.seed(seed)

start.time <- Sys.time()

folder.name <- paste("Outputs for (seed = ",seed,") (",substring(format(start.time, "%Y%m%d %H%M"),1,16),")")
# 2. Create the folder if it does not already exist
if (!dir.exists(folder.name)) {
  dir.create(folder.name)
}

for(i in 1:nrow(hyper.param.df)){
  N.order <- hyper.param.df[i,"d"]
  alpha        <- hyper.param.df[i,"alpha"]
  random.fn <- get.random.fn(hyper.param.df[i,"RV"])
  m <- hyper.param.df[i,"m"]
  
  sqrt.eigenvalues <- sqrt(get.eigenvalues(m = m, n = N.order, alpha = alpha))
  W        <- get.weight.matrix(sqrt.eigenvalues, get.basis.matrix.orthonormal(N.order, t.grid))
  k.sigma  <- get.k.sigma.fixed(alpha, N.order, m = m)
  b.sigma  <- make.b.fn(m, k.sigma)
  a.seq <- a.fn(log(N.seq))
  b.seq <- b.sigma(log(N.seq))
  
  
  ##<< vapply pre-allocates; no growing objects.
  
  res <- vapply(seq_len(N.runs),
                
                function(r) get.H.N.seq.fast(N.seq, W, random.fn, a.seq, b.seq),
                
                matrix(0, nrow = 2, ncol = length(N.seq)))
  
  
  
  T.N.vals.seq  <- res["T", , ]                     # length(N.seq) x N.runs
  H.N.vals.seq  <- res["H", , ]
  
  ks.dist <- vapply(seq_along(N.seq),
                    
                    function(x) unname(ks.test(H.N.vals.seq[x, ], p.limiting)$statistic),
                    
                    numeric(1))
  
  output <- rbind(output, cbind(m,N.order,alpha,random.func[hyper.param.df[i,"RV"]], N.seq, ks.dist))
  output.df <- data.frame(output)
  colnames(output.df) <- c("m","d","alpha","random var","N", "KS dist")
  write.csv(output.df,paste("Output table (",substring(format(start.time, "%Y%m%d %H%M"),1,16),") (seed ",seed,") (Runs = ",N.runs,").csv",sep=""))
  
  saveRDS(res, file =  paste(folder.name,"//HN values (d = ",N.order,")(alpha = ",alpha,") (rv = ",random.func[hyper.param.df[i,"RV"]],") (m = ",m,").rds",sep=""))
}

print(Sys.time() - start.time)



## ---------------------------------------------------------------------

## 5.  Optional: run the replicates in parallel  (Linux / macOS)

## ---------------------------------------------------------------------

if (FALSE) {
  
  library(parallel)
  
  RNGkind("L'Ecuyer-CMRG"); set.seed(2024)
  
  res <- simplify2array(mclapply(seq_len(N.runs),
                                 
                                 function(r) get.H.N.seq.fast(N.seq, W, random.fn, a.seq, b.seq),
                                 
                                 mc.cores = detectCores() - 1))
  
}



# limiting quantiles

limiting.quantiles <- function(q) -log(-log(q)/2)

limiting.quantiles(c(0.9,0.95,0.99))

# check how many in a vector are larger than vector of values

count.greater <- function(x, v){
  sapply(v, FUN = function(val) sum(x >= val))
}


count.greater(curr.H[1,] , quantile.values[2])
