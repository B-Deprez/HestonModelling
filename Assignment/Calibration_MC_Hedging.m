%% Calibration Heston
% Initialising all parameters
S0 = 92.40; % Initial stock price
q = log((S0+4*1.29)/S0);

% The interest rates
r_t = [1/12, 2/12, 3/12, 6/12,    1,    2,    3,    5]; 
r_y = [0.10, 0.09, 0.10, 0.12, 0.18, 0.20, 0.24, 0.36];
r = @(t) spline(r_t, r_y, t)/100; % function for interest rates at any moment expressed as r*100%

% Plotting the interest
plot(0:0.01:5.5, r(0:0.01:5.5),'LineWidth', 2)
xlabel('Years')
ylabel('Interest')
title('Yield Curve')

data = xlsread('Chevron.xlsx');
K = data(:,1);
T = data(:,4)/365; %Maturities are expressed as days in dataset
P_C_flag= data(:,5);
price = data(:,6);
int = r(T);

%% Calibration
% The resulting point is of the form: kappa, nu, theta, v0, rho
cal = fmincon( @(p) sumofsquares(S0, T, q, int, p(4), K, p(1), p(2), p(3), p(5), P_C_flag, price), [0.5, 0.05, 0.2, 0.05, -0.75], [],[], [],[], [0, 0, 0, 0, -1], [Inf, Inf, Inf, Inf, 1], @(p) mycon(p(1),p(2),p(3), p(4), p(5)));
%% Plot Goodness of Fit
Heston_model = @(K, T, P_C) HestonPricer(S0,T, q, r(T), cal(4), K, cal(1), cal(2), cal(3), cal(5), P_C);
hest = 0; %Here the values according to the model will be stored, 0 is just a dummy 
for i = 1:length(T)
    hest(i) = Heston_model(K(i), T(i), P_C_flag(i));
end

scatter(K, price, 'r*')
hold on 
scatter(K, hest, 'b*') 
for i = 1:length(T)
    plot([K(i), K(i)], [price(i), hest(i)], 'k-')
end
xlabel('K')
ylabel('Price')
title('Calibration: Heston')
legend('Market', 'Model')

%% Monte Carlo
% We first do a (rough) simulation for different durations and barriers, to
% find desirable values for these two. 
% Initialise the variables and matrices
T = 5/365:5/365:3;
H_start = S0*0.70;
H_steps = S0*0.001;
H = H_start:H_steps:S0;
m = 1000;

kappa = cal(1);
nu = cal(2); 
theta = cal(3);
v0 = cal(4);
rho = cal(5);

price_prod = zeros(length(T), length(H));
Y = zeros(length(T), length(H)); %We also incorporate what we can borrow
for t = 1:length(T)
    n = round(365*T(t)); %Round is needed since n is not always an integer due to rounding errors, e.g. n = 6 was once 6.0000, hence floating point
    dt = T(t)/n;
    S = zeros(m, n+1);
    S(:,1) = S0;
    V = zeros(m, n+1);
    V(:,1) = v0;
    
    %Set the interest rate for that period
    int = r(T(t));
    
    % Get the random variables 
    eps1 = randn(m,n); 
    eps2 = randn(m,n);
    eps2 = rho*eps1+ sqrt(1-rho^2)*eps2;

    % M-C simulation using Milstein with partial truncation
    for j = 1:n
        V(:,j+1) = V(:,j) +(kappa*(nu-V(:,j)) -(theta^2)/4)*dt +theta *sqrt(max(0,V(:,j))*dt).*eps2(:,j)+(theta^2)/4 * dt * (eps2(:,j)).^2;
    end

    for i = 1:n
        S(:,i+1) = S(:,i).*(1+(int-q)*dt+sqrt(max(0,V(:,i))*dt).*eps1(:,i));
    end
    
    for h = 1:length(H)
        Barrier = H(h); 
        values = zeros(m,1);
        for i = 1:m
            L = min(S(i,:));
            values(i) =  ( (S0-L)/S0 )*max((L-Barrier)/abs(L-Barrier),0);
        end
        price_prod(t,h) = exp(-int*T(t))*mean(values);
        Y(t,h) = (exp(q*T(t))-1)/exp(int*T(t));
    end
end

%% We plot what we can borrow minus what the product costs
contourf(H, T, Y-price_prod, 'ShowText','on')
xlabel('Barrier')
ylabel('Duration')
title('Profit when Borrowing Maximal Amount')

%% M-C simulation with smaller dt and more steps
t = 2;
n = 365*t;
m = 1000000;
dt = t/n;
S = zeros(m, n+1);
S(:,1) = S0;
V = zeros(m, n+1);
V(:,1) = v0;
q = log((S0+4*1.29)/S0);
%q = log((S0+2*1.29)/S0); % To be conservative with the corona virus
   
%Set the interest rate for that period
int = r(t);
    
% Get the random variables 
eps1 = randn(m,n); 
eps2 = randn(m,n);
eps2 = rho*eps1+ sqrt(1-rho^2)*eps2;

% M-C simulation using Milstein with partial truncation
for j = 1:n
    V(:,j+1) = V(:,j) +(kappa*(nu-V(:,j)) -(theta^2)/4)*dt +theta *sqrt(max(0,V(:,j))*dt).*eps2(:,j)+(theta^2)/4 * dt * (eps2(:,j)).^2;
end

for i = 1:n
S(:,i+1) = S(:,i).*(1+(int-q)*dt+sqrt(max(0,V(:,i))*dt).*eps1(:,i));
end


Barrier = 70; 
values = zeros(m,1);
for i = 1:m
    L = min(S(i,:));
    values(i) =  ( (S0-L)/S0 )*max((L-Barrier)/abs(L-Barrier),0);
end
price = exp(-int*t)*mean(values);

% This next value is also the commisson in prof*100% since it is 
% calculated for an investment of 1.
prof = (exp(q*t)-1)/exp(int*t) - price
price

%% Hedging
price1 = price*S0 %The price from above simulation,
% Needs to be multiplied with S0 since this is taken out of pay-off
% structure (see text)

h = 0.5;

t = 2;
n = 365*t;
m = 1000000;
dt = t/n;
S = zeros(m, n+1);
S(:,1) = S0+h;
V = zeros(m, n+1);
V(:,1) = v0;
   
%Set the interest rate for that period
int = r(t);
    
% Get the random variables 
eps1 = randn(m,n); 
eps2 = randn(m,n);
eps2 = rho*eps1+ sqrt(1-rho^2)*eps2;

% M-C simulation using Milstein with partial truncation
for j = 1:n
    V(:,j+1) = V(:,j) +(kappa*(nu-V(:,j)) -(theta^2)/4)*dt +theta *sqrt(max(0,V(:,j))*dt).*eps2(:,j)+(theta^2)/4 * dt * (eps2(:,j)).^2;
end

for i = 1:n
S(:,i+1) = S(:,i).*(1+(int-q)*dt+sqrt(max(0,V(:,i))*dt).*eps1(:,i));
end


Barrier = 70; 
values = zeros(m,1);
for i = 1:m
    L = min(S(i,:));
    values(i) =  ( (S0+h-L))*max((L-Barrier)/abs(L-Barrier),0);
end
price2 = exp(-int*t)*mean(values)

Delta = (price2-price1)/h
Delta2 = Delta*1000000/S0
Delta_total = 1000000/S0 + Delta2