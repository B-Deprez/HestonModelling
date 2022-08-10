%% Illustration of the pay-off
S0 = 92.40;
H = 70;
L = 80; 
ST = 60:0.1:95;

minL = min(L, ST);
Payoff = ST + (S0-minL).*max((minL-H)./abs(minL-H),0);
plot(ST, ST, 'LineWidth', 2)
hold on 
plot(ST, Payoff, 'LineWidth', 2)
xlabel('S_T')
ylabel('Pay-off')
legend('Direct Investment', 'Tempo Certificate')