function [price] = HestonPricer(S0,T, q, r, v0, K, kappa, nu_H, theta, rho_H, P_C_flag)
% HESTONPRICER  The price under the Heston model using the FFT and 
%   Carr-Madan formula with Simpson's rule.

nu = 0.25;
N = 4096; 
alpha = 1.5;

lambda = 2*pi/(N*nu);
b = N*lambda/2;

j = 1:N;
k = -b + lambda*(j-1);
v = (j-1)*nu;

% If there are duplicate greek letters, e.g., both C-M and Heston use "rho",
% then the variables used in the characteristic function have an additional
% H (for Heston)
rho = (exp(-r.*T).*PhiHeston(v-(alpha+1)*1i, T, S0, r,q, v0, kappa, nu_H, theta, rho_H))./(alpha^2+alpha-v.^2+1i*(2*alpha+1)*v);

delta = [1; zeros(N-1,1)]';
A = exp(1i.*v.*b).*rho.*nu.*( (3+(-1).^j-delta)./3);
FFT = fft(A, N);
a = real(FFT);
calls = 1/pi * exp(-alpha*k).*a; %Price of a call
prices = calls -P_C_flag.*(exp(-q*T)*S0 - exp(-r.*T).*K); %Use Put-Call Parity to obtain the right price depending on de P_C_flag
KK = exp(k);
price = spline(KK, prices, K);