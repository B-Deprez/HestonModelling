library(readxl)

dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)
source("HestonPricer.R")

## Calibration Heston
# Initialising all parameters
S0 = 92.40 # Initial stock price
q = log((S0+4*1.29)/S0)

# The interest rates
r_t = c(1/12, 2/12, 3/12, 6/12,    1,    2,    3,    5) 
r_y = c(0.10, 0.09, 0.10, 0.12, 0.18, 0.20, 0.24, 0.36)
r = splinefun(r_t, r_y/100) # function for interest rates at any moment expressed as r*100#
  
# Plotting the interest
x_plot = seq(0,5.5, by = 0.01)
plot(x_plot, r(x_plot)*100,
     main = "Yield Curve", 
     xlab = 'Years',
     ylab = 'Interest (%)',
     type = "l",
     lwd = 2, 
     col = "blue")

# Setting up the data
data = read_xlsx('Chevron.xlsx')

K = data[,1][[1]]
Tt = data[,4][[1]]/365 #Maturities are expressed as days in dataset. We need years
P_C_flag= data[,5][[1]]
price = data[,6][[1]]
int = r(Tt)

