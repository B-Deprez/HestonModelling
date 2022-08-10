function sum = sumofsquares(S0, T, q, r, v0, K, kappa, nu_H, theta, rho_H, P_C_flag, price)
sum = 0;
for t = 1:length(T)
    sum = sum + ( HestonPricer(S0,T(t), q, r(t), v0, K(t), kappa, nu_H, theta, rho_H, P_C_flag(t)) - price(t) )^2;
end
sum = sum/length(T);
