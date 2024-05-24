data{
    int<lower=1> L;
    int E[L];
    vector[L] months_to_event;
    real A[L];
    int F[L];
    int d[L];
}
parameters{
    matrix[2,2] z;
    real a_bar;
    real beta;
    real<lower=0> sigma;
}
model{
    vector[L] lambda;
    vector[L] mu;
    sigma ~ exponential( 1 );
    a_bar ~ normal( 0 , 1 );
    beta ~ normal( 0 , 1 );
    to_vector( z ) ~ normal( 0 , 1 );
    for ( i in 1:L ) {
        mu[i] = a_bar + z[d[i], F[i]] * sigma + beta * A[i];
        mu[i] = exp(mu[i]);
        lambda[i] = 1/mu[i];
    }
    for ( i in 1:L ) 
        if ( E[i] == 0 ) target += exponential_lccdf(months_to_event[i] | lambda[i]);
    for ( i in 1:L ) 
        if ( E[i] == 1 ) months_to_event[i] ~ exponential( lambda[i] );
}
generated quantities{
    matrix[2,2] a;
    a = a_bar + z * sigma;
}
