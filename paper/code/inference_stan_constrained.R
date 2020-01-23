require(rstan)
require(phytools)

rstan_options(mc.cores = parallel::detectCores(),auto_write = TRUE)

IIr.trees <- read.nexus('IIr.trees')


states <- read.csv('states.csv',header=F,sep=' ',row.names=1)
levels(states$V2) <- c(levels(states$V2),'(-clf,+obl.pl)','(-clf,-obl.pl)','(+clf,+obl.pl)','(+clf,-obl.pl)')
states[states$V2=='(-clf,morph.pl)',] <- '(-clf,+obl.pl)'
states[states$V2=='(-clf,obl.pl)',] <- '(-clf,+obl.pl)'
states[states$V2=='(-clf,opt.pl)',] <- '(-clf,-obl.pl)'
states[states$V2=='(+clf,obl.pl)',] <- '(+clf,+obl.pl)'
states[states$V2=='(+clf,opt.pl)',] <- '(+clf,-obl.pl)'
states <- droplevels(states)
states$V2 <- factor(states$V2,levels=sort(levels(states$V2)))


model.code <- "data {
  int<lower=1> N; //number of tips+internal nodes+root
  int<lower=1> T; //number of tips
  int<lower=1> B; //number of branches
  int<lower=1> F;                       //number of states
  int<lower=1> child[B];                //child of each branch
  int<lower=1> parent[B];               //parent of each branch
  real<lower=0> brlen[B];                //length of each branch
  int<lower=0,upper=1> tiplik[T,F];     //likelihoods for data at tips in tree
  }
parameters {
  real<lower=0> R[F*(F-2)];             //rates to be put into matrix
  }
transformed parameters {
  matrix[F,F] Q;                       //rate matrix
  //fill relevant cells of rate matrix
  //Q[1,1] = 0;
  Q[1,2] = R[1];
  Q[1,3] = R[2];
  Q[1,4] = 0;
  Q[2,1] = R[3];
  //Q[2,2] = 0;
  Q[2,3] = 0;
  Q[2,4] = R[4];
  Q[3,1] = R[5];
  Q[3,2] = 0;
  //Q[3,3] = 0;
  Q[3,4] = R[6];
  Q[4,1] = 0;
  Q[4,2] = R[7];
  Q[4,3] = R[8];
  //Q[4,4] = 0;
  //fill diagonal rows of rate matrix
  for (i in 1:F) {
    real z;
    z = 0;
    for (j in 1:F) {
      if (i != j) {
        z = z - Q[i,j];
      }
    }
    Q[i,i] = z;
  }
}
model {
  matrix[N,F] lambda; //matrix of likelihoods for root, internal nodes, and tips
  matrix[F,F] pi_matrix;                //holds approximation of stationary probabilities in each row
  vector[F] pi;                         //stationary probability
  for (i in 1:F*(F-2)) {
    R[i] ~ gamma(1,5);                //shape/rate parameterization
    //R[i] ~ uniform(0,1);
  }
  for (t in 1:T) {
    for (f in 1:F) {
      lambda[t,f] = log(tiplik[t,f]);
    }
  }
  for (n in (T+1):N) {
    for (f in 1:F) {
      lambda[n,f] = 0;
    }
  }
  for (b in 1:B) {
    matrix[F,F] P;
    P = matrix_exp(brlen[b]*Q); //via matrix exponentiation
    for (f in 1:F) {
      lambda[parent[b],f] = lambda[parent[b],f] + log(dot_product(P[f],exp(lambda[child[b]])));
    }
  }
  pi_matrix = matrix_exp(100000*Q);     //
  for (f in 1:F) {                      //dirty way of approximating
    pi[f] = pi_matrix[1,f];             //stationary probabilities of
  }                                     //each character state;
  target += log(dot_product(pi,exp(lambda[parent[B]])));   //increment log posterior by tree likelihood
}"

fit.list <- list()
for (i in floor(seq(10,length(IIr.trees),length.out=100))) {
    tree <- IIr.trees[[i]]
    #tree <- IIr.trees[[i]]
    tree <- reorder.phylo(tree,'pruningwise')
    bin.states <- states[tree$tip.label,]
    bin.states <- to.matrix(bin.states,seq=levels(bin.states))
    parent <- tree$edge[,1]
    child <- tree$edge[,2]
    b.lens <- tree$edge.length
    N <- length(unique(c(parent,child)))
    T <- length(child[which(!child %in% parent)])
    data.list <- list(N=N,
              T=T,
              B=length(parent),
              brlen=b.lens/1000,
              child=child,
              parent=parent,
              tiplik=bin.states,
              F=ncol(bin.states))
    fit <- stan(model_code=model.code,data=data.list,thin=5,chains=3)
    fit.list[[length(fit.list)+1]] <- fit
}


fit.full <- sflist2stanfit(fit.list)

#save.image(file = 'all_rates_constrained.Rdata')
save.image(file = 'all_rates.Rdata')