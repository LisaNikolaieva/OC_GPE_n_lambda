function gradJ = gradJ(lambda,gradJ_par)

% lambda = 1 1 1 1 1    
%          2 2 2 2 2 
n_l = size(lambda,1);

grid = gradJ_par.grid;
pgrid = gradJ_par.pgrid;
par = gradJ_par.par;
V = gradJ_par.V;
dV_dlambda_fun = gradJ_par.dV_dlambda_fun;
pT_fun = gradJ_par.pT_fun;
int_vec_fun = gradJ_par.int_vec_fun;
u0 = gradJ_par.u0;


%%
[Psi_store] = Psi_xt(u0,grid,par,V,lambda);

        figure(21)
        plot(grid.x,abs(Psi_store(:,end)).^2)
        set(gca,'XMinorGrid','on');
        set(gca,'YMinorGrid','on');


% plambda = interp1(grid.t,lambda,pgrid.t,'linear','extrap');
plambda = zeros(n_l,pgrid.Nt);
for i = 1:n_l
   plambda(i,:) = interp1(grid.t,lambda(i,:),pgrid.t,'linear','extrap');
end

pPsi_store = interp2(grid.t_mesh,grid.x_mesh,Psi_store,pgrid.t_mesh,pgrid.x_mesh,"nearest",0);                       % psi_p = interp2(grid.t_mesh,grid.x_mesh,Psi_store(:,2:end),pgrid.t_mesh,pgrid.x_mesh,"nearest",0);

pT = pT_fun(Psi_store(:,end));
[p_store] = p_xt(pT,pgrid,par,V,plambda,pPsi_store);

%%
int_vec = zeros(n_l,pgrid.Nt);
pdV_dlambda = zeros(pgrid.N,pgrid.Nt);


for i = 1:n_l
    for j = 1:pgrid.Nt
        pdV_dlambda(:,j) = dV_dlambda_fun(plambda(:,j),i);
    end
    int_vec(i,:) = int_vec_fun(pPsi_store,pdV_dlambda,p_store);
end

        figure(1)
        subplot(4,1,2)
        hold on
        plot(pgrid.t,trapz(grid.x,real(pdV_dlambda),1))
        subplot(4,1,4)
        hold on
        plot(pgrid.t,trapz(grid.x,real(conj(pPsi_store).*pdV_dlambda.*p_store),1))

        subplot(4,1,1)
        plot(pgrid.t,trapz(grid.x,real(conj(pPsi_store)),1))
        hold on
        subplot(4,1,3)
        plot(pgrid.t,trapz(pgrid.x,real(p_store),1))
        hold on

gradJ = par.gamma*plambda' + (pgrid.lap4_t\(int_vec'));
gradJ = gradJ./max(abs(gradJ),[],1);



end