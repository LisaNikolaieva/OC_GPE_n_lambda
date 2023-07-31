function p_store = p_xt(pT,pgrid,par,V,plambda,pPsi_store)

p = pT;

mu_fun = @(lambda_t0,Psi_xt0) 0*4.3242 + 0*3.3822 + real(...
   sum(( -0.5*pgrid.lap4/par.m + spdiag(V(lambda_t0)+par.g*(abs(Psi_xt0).^2)) )*Psi_xt0.*conj(Psi_xt0) )...
   /sum(abs(Psi_xt0.^2)) );




%%
p_store = zeros(pgrid.N,pgrid.Nt);
p_store(:,pgrid.Nt) = pT;

%%
fprintf('start p\n')
for i = pgrid.Nt-1:-1:1
    
            p_kp1 = p;
            pPsi_k = pPsi_store(:,i);
            pPsi_kp1 = pPsi_store(:,i+1);
            V_k = V(plambda(:,i));
            V_kp1 = V(plambda(:,i+1));
            mu_k = mu_fun(plambda(:,i),pPsi_store(:,i));
            mu_kp1 = mu_fun(plambda(:,i+1),pPsi_store(:,i+1));
            
            p = propagate_p_GPE_cn(p_kp1,pPsi_k,pPsi_kp1,abs(pgrid.dt),V_k,V_kp1,mu_k,mu_kp1,pgrid,par);         
            p_store(:,i) = p;
            
end
fprintf('end p\n')


end