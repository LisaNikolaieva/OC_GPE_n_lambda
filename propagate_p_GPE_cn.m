function [p_k] = propagate_p_GPE_cn(p_kp1, psi_k, psi_kp1, dt, V_k, V_kp1, mu_k,mu_kp1, pgrid, par)

m = par.m;
g = par.g;

psi_k_q_real = real(psi_k.^2);
psi_k_q_imag = imag(psi_k.^2);
psi_kp1_q_real = real(psi_kp1.^2);
psi_kp1_q_imag = imag(psi_kp1.^2);

Ak_11 = speye(pgrid.N) + dt/2*g*spdiag(psi_k_q_imag);
Ak_12 = dt/2*(-1/2/m * pgrid.lap4 + spdiag(V_k - mu_k + 2*g*abs(psi_k).^2 - g*psi_k_q_real));
Ak_21 = dt/2*(1/2/m * pgrid.lap4 - spdiag(V_k - mu_k + 2*g*abs(psi_k).^2 + g*psi_k_q_real));
Ak_22 = speye(pgrid.N) - dt/2*g*spdiag(psi_k_q_imag);

Akp1_11 = speye(pgrid.N) - dt/2*g*spdiag(psi_kp1_q_imag);
Akp1_12 = -dt/2*(-1/2/m * pgrid.lap4 + spdiag(V_kp1 - mu_kp1 + 2*g*abs(psi_kp1).^2 - g*psi_kp1_q_real));
Akp1_21 = -dt/2*(1/2/m * pgrid.lap4 - spdiag(V_kp1 - mu_kp1 + 2*g*abs(psi_kp1).^2 + g*psi_kp1_q_real));
Akp1_22 = speye(pgrid.N) + dt/2*g*spdiag(psi_kp1_q_imag);

M_11 = Ak_11-Ak_12*spdiag((1./diag(Ak_22)))*Ak_21;
M_22 = Ak_22-Ak_21*spdiag((1./diag(Ak_11)))*Ak_12;
V_1 = Akp1_11 * real(p_kp1) + Akp1_12 * imag(p_kp1);
V_2 = Akp1_21 * real(p_kp1) + Akp1_22 * imag(p_kp1);

sol =   [M_11\V_1-M_11\(Ak_12*(spdiag(1./diag(Ak_22))*V_2));-M_22\(Ak_21*(spdiag(1./diag(Ak_11))*V_1))+M_22\V_2];


p_k = sol(1:pgrid.N)+1i*sol(pgrid.N+1:2*pgrid.N);

end
