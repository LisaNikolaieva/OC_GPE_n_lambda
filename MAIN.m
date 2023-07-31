close all;
Config;

%%

f_fun = @(X) cost_function(X,cf_par); % f:      X |-> R
gradf_fun = @(X) gradJ(X,gradJ_par);  % gradf:  X |-> x
X0 = lambda0;
X2x = @(X) X2x_matrix(X,grid,pgrid); % X2x:  X |-> x
x2X = @(x) x2X_matrix(x,grid,pgrid); % x2X:  x |-> X
n_steps = 15;

[lambda_store, cost_function_store] = find_min_BFGS2(f_fun,gradf_fun,X0,X2x,x2X,n_steps);



%%
function x = X2x_matrix(X,grid,pgrid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input:            output:
% X = 1 1 1 1 1     x = 1 1 1 
%     2 2 2 2 2         2 2 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = pgrid.Nt;
n_l = size(X,1);

x = zeros(n_l,n);

for i = 1:n_l
   x(i,:) = interp1(grid.t,X(i,:),pgrid.t,'linear','extrap');
end


end




function X = x2X_matrix(x,grid,pgrid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input:         output:
% x = 1 1 1      X = 1 1 1 1 1
%     2 2 2          2 2 2 2 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = grid.Nt;
n_l = size(x,1);

X = zeros(n_l,N);

for i = 1:n_l
   X(i,:) = interp1(pgrid.t,x(i,:),grid.t,'linear','extrap');
end


end