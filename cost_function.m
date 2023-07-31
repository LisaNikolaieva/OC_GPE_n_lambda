function J = cost_function(lambda,cf_par)

% lambda = 1 1 1 1 1    
%          2 2 2 2 2 

cost_fun_fun = cf_par.cost_fun_fun;
grid = cf_par.grid;
par = cf_par.par;
u0 = cf_par.u0;
V = cf_par.V;



[Psi_store] = Psi_xt(u0,grid,par,V,lambda);


J_out = cost_fun_fun(Psi_store(:,end),lambda);
J = J_out(1);
J_E = J_out(2);
J_lambda = J_out(3);
fprintf('J = %f; J_E/S = %f; J_lambda = %f\n',J,J_E,J_lambda)
end