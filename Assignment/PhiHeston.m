function [value] = PhiHeston(u,T,S0,r,q, v0, kappa, nu, theta, rho)
% PHIHESTON  The characteristic function for the Heston model
%   V = PHIHESTON(u,T,S0,r,q, v0, kappa, nu, theta, rho) returns the value
%   of the characteristic function where the entries can both be numeric
%   or arrays.

sigma = sqrt(v0);
d = sqrt( (rho*theta*u*1i-kappa).^2-theta^2*(-1i*u-u.^2) );
g = (kappa-rho*theta*1i*u-d)./(kappa-rho*theta*1i*u+d);
value1 = exp(1i*u*(log(S0)+(r-q)*T));
value2 = exp(nu.*kappa.*theta^(-2).*( (kappa-rho.*theta.*1i.*u-d).*T - 2*log( (1-g.*exp(-d*T))./(1-g))));
value3 = exp(sigma^2*theta^(-2)* (kappa-rho*theta*1i*u-d) .* (1-exp(-d*T))./(1-g.*exp(-d*T)));
value = value1.*value2.*value3;