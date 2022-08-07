PhiHeston <- function(u,Tt,S0,r,q, v0, kappa, nu, theta, rho){
  #The characteristic function for the Heston model
  #   V = PHIHESTON(u,T,S0,r,q, v0, kappa, nu, theta, rho) returns the value
  #   of the characteristic function where the entries can both be numeric
  #   or arrays.
  sigma = sqrt(v0)
  d = sqrt( (rho*theta*u*1i-kappa)^2-theta^2*(-1i*u-u^2) )
  g = (kappa-rho*theta*1i*u-d)/(kappa-rho*theta*1i*u+d)
  value1 = exp(1i*u*(log(S0)+(r-q)*Tt))
  value2 = exp(nu*kappa*theta^(-2)*( (kappa-rho*theta*1i*u-d)*Tt - 2*log( (1-g*exp(-d*Tt))/(1-g))))
  value3 = exp(sigma^2*theta^(-2)* (kappa-rho*theta*1i*u-d) * (1-exp(-d*Tt))/(1-g*exp(-d*Tt)))
  value = value1*value2*value3
}


HestonPricer <- function(S0,Tt, q, r, v0, K, kappa, nu_H, theta, rho_H, P_C_flag){
  # Calculate the price of a put or call option according to the Heston model
  # Uses the FFT and Carr-Madan formula with Simpson's rule
  nu = 0.25
  N = 4096 # Take power of 2
  alpha = 1.5
  
  lambda = 2*pi/(N*nu)
  b = N*lambda/2
  
  j = 1:N
  k = -b + lambda*(j-1)
  v = (j-1)*nu
  
  # If there are duplicate greek letters, e.g., both C-M and Heston use "rho",
  # then the variables used in the characteristic function have an additional
  # H (for Heston)
  rho = (exp(-r*Tt)*PhiHeston(v-(alpha+1)*1i, Tt, S0, r,q, v0, kappa, nu_H, theta, rho_H))/(alpha^2+alpha-v^2+1i*(2*alpha+1)*v)
  delta <- c(1, rep(0,N-1))
  A = exp(1i*v*b)*rho*nu*( (3+(-1)^j-delta)/3)
  FFTr = fft(A)
  a = Re(FFTr)
  calls = 1/pi * exp(-alpha*k)*a #Price of a call
  prices = calls
  KK = exp(k)
  f = splinefun(KK, prices)
  price = f(K) -P_C_flag*(exp(-q*Tt)*S0 - exp(-r*Tt)*K); #Use Put-Call Parity to obtain the right price depending on de P_C_flag
  return(price)
}
