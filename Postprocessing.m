Config;
%%
load('lambda_store_cost_function_store2.mat','lambda_store','cost_function_store');


lambda_lin = lambda_store(:,:,1);
lambda_opt = lambda_store(:,:,end);

Psi_store_lin = Psi_xt(u0,grid,par,V,lambda_lin);
Psi_store_opt = Psi_xt(u0,grid,par,V,lambda_opt);

den_opt = abs(Psi_store_opt).^2;
den_lin = abs(Psi_store_lin).^2;
pd = max([den_opt,den_lin],[],'all');



figure
subplot(4,1,1)
plot(grid.t,lambda_lin(1,:))
hold on
grid on
plot(grid.t,lambda_opt(1,:))
title('lambda')

subplot(4,1,2)
plot(grid.t,lambda_lin(2,:))
hold on
grid on
plot(grid.t,lambda_opt(2,:))
title('lambda')

subplot(4,1,3)
imagesc(grid.t,grid.x,den_lin)
caxis([0,pd])
title('linear')

subplot(4,1,4)
imagesc(grid.t,grid.x,den_opt)
caxis([0,pd])
title('optimized')



if 1
%% plot density
    for i = 1:50:grid.Nt
        figure(26)
        subplot(2,1,1)
        plot(grid.x,den_lin(:,i),grid.x,den_T,'--',grid.x,den_0,'--')
        grid on
        ylim([0, pd])
        title(sprintf('linear: t = %.0f ms',grid.t(i)))
        drawnow
        subplot(2,1,2)
        plot(grid.x,den_opt(:,i),grid.x,den_T,'--',grid.x,den_0,'--' )
        grid on
        ylim([0, pd])
        title(sprintf('optimized: t = %.0f ms',grid.t(i)))
        

    end

end




state_cost = @(Psi) 1/2*(1-abs(trapz(grid.x,conj(uT).*Psi)).^2);

state_cost_opt = zeros(1,grid.Nt);
state_cost_lin = zeros(1,grid.Nt);


for i = 1:grid.Nt
    state_cost_lin(i) = state_cost(Psi_store_lin(:,i));
    state_cost_opt(i) = state_cost(Psi_store_opt(:,i));
   
    
end

state_max = max([state_cost_lin,state_cost_opt]);
state_min = min([state_cost_lin,state_cost_opt]);


figure
subplot(1,2,1)
semilogy(grid.t,state_cost_lin);
hold on
grid on
semilogy(grid.t,state_cost_opt)
semilogy(grid.t,state_max*ones(size(grid.t)),'--')
legend('linear','optimized')
title('state cost (t)')
ylim([state_min*0.9,state_max*1.1])

subplot(1,2,2)
semilogy(cost_function_store)
grid on
title('energy cost (i)')

fprintf('cost function reduction factor: %f\n',cost_function_store(1)/cost_function_store(end))






