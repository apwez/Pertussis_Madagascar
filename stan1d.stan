data {
  int<lower=0> N;  		//number of individuals
  int<lower=0> T;  		//number of years+1 for initial year
  int<lower=0> nloc;  		//number of sites
  matrix[N,T] x;  		//fraction each year individual present
  int<lower=0,upper=1> y[N];    //if seropositive or seronegative IgG
  int<lower=0,upper=nloc> locs[N];    //if seropositive or seronegative IgG

  real lambdaStart;		// Prior on lambda
  real lambdaStartSigma;	// Prior on lambda

  real vStart[nloc];			// Prior on v
  real vStartSigma[nloc];	        // Prior on v


}
parameters {
  matrix<lower=0.0001,upper=5>[nloc,T] lambda;       		// coeff of year with location adjustment
  vector<lower=0,upper=0.8>[nloc] v; 		//prob was vaccinated individual level
}

model {

   // get notFoi/notVacc for IgG seropos/neg data
   vector[N] estSeroNeg;

    // Define the priors
    //lambda ~ normal(lambdaStart,lambdaStartSigma);	// Mean transmission, by biweek
    v ~ normal(vStart,vStartSigma);			// Mean rate vacc, by biweek
    
   for (n in 1:N)
        estSeroNeg[n] = exp(-sum(x[n].*lambda[locs[n]]))*(1-v[locs[n]]);
	
  // likleihoods
  y ~ bernoulli(1-estSeroNeg);

 
}


generated quantities {

   vector[N] logLik;
  vector[N] estSeroNeg;

    
   for (n in 1:N)
        estSeroNeg[n] = exp(-sum(x[n].*lambda[locs[n]]))*(1-v[locs[n]]);

   for (n in 1:N)
        logLik[n] = bernoulli_lpmf(y[n]|1-estSeroNeg[n]);

 
}

