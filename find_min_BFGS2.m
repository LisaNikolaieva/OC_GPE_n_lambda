function [X_store, f_store] = find_min_BFGS2(f_fun,gradf_fun,X0,X2x,x2X,n_steps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input:
% X0 = 1 1 1 1 1    ->    x0 = 1 1 1  
%      2 2 2 2 2               2 2 2
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = X2x(X0);

N = size(X0,2);    % 5 
n = size(x0,2);    % 3
n_l = size(X0,1);  % 2                   2 2 2 
B0 = repmat(eye(n),[1,1,n_l]); % 3x3x2 1 1 1 2
                               %       1 1 1 2
                               %       1 1 1


X_store = zeros(n_l,N,n_steps+1);
X_store(:,:,1) = X0;

f_store = zeros(1,n_steps+1);
f_store(1) = f_fun(X0);


Xkm1 = X0;
Bkm1 = B0;
skm1 = zeros(n_l,n); %  1 1 1
                     %  2 2 2
is_initial = 1;

for i = 1:n_steps

  tic
 [Bk,sk,Xk,fk] = make_BFGS_step(Bkm1,skm1,Xkm1,f_fun,gradf_fun,x2X,is_initial);
  
 f_store(:,i+1) = fk;
 X_store(:,:,i+1) = Xk;
  
 Bkm1 = Bk;
 skm1 = sk; 
 Xkm1 = Xk;

  is_initial = 0;
  tim = toc();
  fprintf('step %i done in %f s\n',i,tim)
end


end
   
function [Bk,sk,Xk,fk] = make_BFGS_step(Bkm1,skm1,Xkm1,f_fun,gradf_fun,x2X,is_initial)   

n = size(skm1,2);
n_l = size(skm1,1);

grad_J = gradf_fun(Xkm1); %     1 2 
                          %     1 2
                          %     1 2
plot_xgradf(Xkm1,grad_J,n_l);

ykm1 = grad_J;            %     1 2
                          %     1 2
                          %     1 2



if is_initial
    Bk = repmat(eye(n),[1,1,n_l]);
else
    for i = 1:n_l
    Sk = skm1(i,:)'*skm1(i,:);
    SYk = skm1(i,:)'*ykm1(:,i)';
    yksk = ykm1(:,i)'*skm1(i,:)';
    Bk(:,:,i) = (speye(n)-sparse(SYk)/yksk)*Bkm1(:,:,i)*(speye(n)-sparse(SYk')/yksk)+Sk/yksk;
    end
end


pk = zeros(n_l,n);
for i = 1:n_l
pk(i,:) = -Bk(:,:,i)\grad_J(:,i);
pk(i,:) = pk(i,:)/max(abs(pk(i,:)));
end

Pk = x2X(pk);
for i = 1:n_l 
Pk(i,:) = Pk(i,:) - linspace(Pk(i,1),Pk(i,end),size(Pk,2));
pk(i,:) = pk(i,:) - linspace(pk(i,1),pk(i,end),size(pk,2));
end

cost_fun = @(a) f_fun(Xkm1 + a.*Pk);


linsearch_options = optimoptions("fminunc",'Display','iter-detailed');
[a_opt,fk] = fminunc(cost_fun,eps*ones(n_l,1),linsearch_options);
fprintf('optimal a: %g\n',a_opt')

sk = a_opt.*pk;

Xk = Xkm1 + a_opt.*Pk;
end


function plot_xgradf(X,gradJ,n_l)
figure(20);
for i = 1:n_l
subplot(2,n_l,n_l+i)
plot(gradJ(:,i))
xlim([-inf inf])
hold on
subplot(2,n_l,i)
plot(X(i,:))
xlim([-inf inf])
hold on
end
drawnow
end




