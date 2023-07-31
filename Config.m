close all
Physical_parameters;

%% grid
grid.N=4*512;
grid.x = linspace(-grid.L/2,grid.L/2,grid.N)';
grid.dx = grid.x(2)-grid.x(1);
grid.lap4 = gallery( 'toeppen', grid.N, - 1,  16, - 30, 16, - 1 ) / ( 12 * grid.dx ^ 2 );
grid.wav = 2 * pi * ( 0 : grid.N - 1 )' / grid.N;
grid.ilap = - 2 * ( 1 - cos( grid.wav ) ) / grid.dx ^ 2;


%% define potential
Vmax        = 2e3;
sigma       = 2;
% V      = @(lambda) Vmax - Vmax * 1/2 * (1 + erf((grid.x + 50 - lambda)/sigma)) + Vmax * 1/2 * ( 1+ erf((grid.x - 50 + lambda)/sigma));
% dV_dlambda_fun  = @(lambda) 1/(1e-3)*(V(lambda+5e-4)-V(lambda-5e-4));
set_V;


mu_fun = @(lambda_t0,Psi_xt0) 0*4.3242 + 0*3.3822 + real(...
   sum(( -0.5*grid.lap4/par.m + spdiag(V(lambda_t0)+par.g*(abs(Psi_xt0).^2)) )*Psi_xt0.*conj(Psi_xt0) )...
   /sum(abs(Psi_xt0.^2)) );


lambda_T = [15;15];
lambda_0 = [6.2;6.2];
n_l = numel(lambda_0);




%% calculate groundstates

u0 = itp(V(lambda_0),par,grid); 
mu_fun(lambda_0,u0)

uT = itp(V(lambda_T),par,grid);
mu_fun(lambda_T,uT)


mu_fun(lambda_T,uT) - mu_fun(lambda_0,u0)
%% estimate speed of sound
den_0 = abs(u0).^2;
den_T = abs(uT).^2;
figure
plot(grid.x,den_T)
hold on
plot(grid.x,0.1*max(den_T)*ones(size(den_T)),'--')
drawnow

c_s1 = sqrt(par.g*max(abs(uT).^2)/par.m); %estimated speed of sound for final state
c_s0 = sqrt(par.g*max(abs(u0).^2)/par.m); %estimated speed of sound for initial state
dt = 0.1*grid.dx/max(c_s1,c_s0);

%% time grid
grid.Nt = ceil(grid.T/dt);
grid.t = linspace(0,grid.T,grid.Nt);
grid.dt = grid.t(2) - grid.t(1);
[grid.t_mesh,grid.x_mesh] = meshgrid(grid.t,grid.x);


%% pgrid
pgrid.L = grid.L;
pgrid.N = grid.N;
pgrid.x = grid.x;
pgrid.dx = grid.dx;
pgrid.lap4 = grid.lap4;
pgrid.nc_factor = 5;
pgrid.T = grid.T;
pgrid.Nt = ceil(grid.Nt/pgrid.nc_factor);
pgrid.t = linspace(0,pgrid.T,pgrid.Nt);
pgrid.dt = pgrid.t(2) - pgrid.t(1);         
pgrid.lap4_t   =  gallery( 'toeppen', pgrid.Nt, - 1,  16, - 30, 16, - 1 ) / ( 12 * pgrid.dt ^ 2 );
[pgrid.t_mesh,pgrid.x_mesh] = meshgrid(pgrid.t,pgrid.x);

%% cost_fun_fun
ham       = @(lambda_t0)-0.5*grid.lap4/par.m+spdiag(V(lambda_t0));
energy_fun = @(Psi_xt0,lambda_t0) real(trapz(grid.x,conj(Psi_xt0).*(ham(lambda_t0)+par.g/2*spdiag(abs(Psi_xt0).^2))*Psi_xt0));
ET = energy_fun(uT,lambda_T);

state_fun = @(Psi_xT) 1/2*(1-abs(trapz(grid.x,conj(uT).*Psi_xT)).^2);

cost = 'energy';
switch cost
    case 'energy'
        phi = @(Psi_xT) energy_fun(Psi_xT,lambda_T)-ET;
        pT_fun = @(Psi_xT) -2*1i*(ham(lambda_T) - spdiag(mu_fun(lambda_T,Psi_xT)*ones(size(Psi_xT))) + par.g*spdiag(abs(Psi_xT).^2))*Psi_xT;
    case 'state'
        phi = @(Psi_xT) state_fun(Psi_xT);
        pT_fun = @(Psi_xT) 1i * trapz(grid.x,conj(uT).*Psi_xT)*uT;
end

stab = @(lambda_t) par.gamma/2*sum(abs(grid.dt)*trapz((diff(lambda_t,2)/abs(grid.dt)).^2,2));
cost_fun_fun = @(Psi_xT,lambda_t) [phi(Psi_xT) + stab(lambda_t); ...
                                                    phi(Psi_xT); ...
                                                stab(lambda_t)];


int_vec_fun = @(pPsi_xt,pdVdl_xt,p_xt) trapz(pgrid.x,real(conj(pPsi_xt).*pdVdl_xt.*p_xt),1);

%% initializing
% lambda0 = linspace(lambda_0,lambda_T,grid.Nt);


lambda0 = zeros(n_l,grid.Nt);
for i = 1:n_l
    lambda0(i,:) = linspace(lambda_0(i),lambda_T(i),grid.Nt);
end

%% parameters

workspace.grid = grid;
workspace.pgrid = pgrid;
workspace.par = par;
workspace.cost_fun_fun = cost_fun_fun;
workspace.pT_fun = pT_fun;
workspace.int_vec_fun = int_vec_fun;
workspace.V = V;
workspace.dV_dlambda_fun = dV_dlambda_fun;
workspace.u0 = u0;


cf_par = workspace;
gradJ_par = workspace;
