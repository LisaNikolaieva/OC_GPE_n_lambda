function Psi_store = Psi_xt(u0,grid,par,V,lambda)

Psi = u0;

mu_fun = @(lambda_t0,Psi_xt0) 0*4.3242 + 0*3.3822 +real(...
   sum(( -0.5*grid.lap4/par.m + spdiag(V(lambda_t0)+par.g*(abs(Psi_xt0).^2)) )*Psi_xt0.*conj(Psi_xt0) )...
   /sum(abs(Psi_xt0.^2)) );


ham       = @(lambda_t0)-0.5*grid.lap4/par.m+spdiag(V(lambda_t0));
energy_fun = @(Psi_xt0,lambda_t0) real(trapz(grid.x,conj(Psi_xt0).*(ham(lambda_t0)+par.g/2*spdiag(abs(Psi_xt0).^2))*Psi_xt0));



%%
Psi_store = zeros(grid.N,grid.Nt);
Psi_store(:,1) = u0;

%%
for i = 2:grid.Nt
    
Psi = propagate_Psi_SS(grid.dt,grid,par,Psi,V(lambda(:,i)),mu_fun(lambda(:,i),Psi));
Psi = Psi./sqrt(trapz(grid.x,abs(Psi).^2));
Psi_store(:,i) = Psi;
            

end


end