%% Parameter
h       = 6.62607015e-34;
h_quer  = h/2/pi;
m       = 1.443e-25; %kg
m       = m /h_quer *(1e-6)^2/1e-3; %m~ = m/(h_quer) *(1mum)�/1ms
omega = 10; %kHz
as = 0.0042; %mum
N = 5000;

par.gamma   = 1e-5;
par.omega   = omega;
par.as      = as;
par.m       = m;
par.N       = N;
par.g         = (2/2.34)^2* 2*N*as*omega;



lambda_0 = 0;
lambda_T = 10;


grid.L = 160;
grid.T = 50;
