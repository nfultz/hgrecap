s <- replicate(4, sample(LETTERS,3), simplify = F)




hgcapmcmc <- function(s, mcmc.control=c()) {

  r <- length(s)
  marked <- Reduce(union, s, c(), accumulate = TRUE)

  N_min <- length(marked[[r + 1]])

  n_drawn   <-  sapply(s, length)
  K_marked   <-  sapply(marked[-r], length)
  d_unmarked <-  sapply(Map(setdiff, s, marked[-r]), length)


  # n_drawn <- n_drawn[2]
  # K_marked <- K_marked[2]
  # d_unmarked <- d_unmarked[2]


  dhyper(d_unmarked, N_min-K_marked, K_marked, n_drawn )

  prior = 0.0
  Q = local({
    q <- c()
    function(i) {
      if(i <= length(q)) return(q[i])
      N <- N_min + i
      q[i] = sum(log(N - n_drawn) - log(N) + log(N-K_marked) - log(N-K_marked - d_unmarked)) - prior
      q[i]
    }
  })

  i = 1



  # q[i] <- log(N_min+i - n_drawn) - log(N_min+i) + log(N_min+i-K_marked) - log(N_min+i-K_marked - d_unmarked)
  out = c()

  for(j in 1:20000){

    a <- runif(1);
    a <- a + a - 1;

    w = sign(a)

    if(i == 0 && w == -1) {
      out <- c(out, N_min)
      next;
    }

    a = abs(a)

    P = if(w == -1) -Q(i) else Q(i+1)

    if(log(a) < P) {
      i = i + w

    }
    out <- c(out, N_min + i)
  }

rcppmm <- hgrecap::rcpp_hgrecap_metropolis(n_drawn, K_marked, d_unmarked, N_min, list(lambda=prior), list(mcmc.iterations=5000000, mcmc.thin=10000, mcmc.burnin=1000, verbose=0))
}
