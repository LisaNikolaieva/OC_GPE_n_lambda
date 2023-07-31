lambda_is = 25;
initial_shift = 10;

V_comp_eta      = @(lambda) Vmax - Vmax * 1/2 * (1 + erf((grid.x + 50 - lambda_is + lambda(1))/sigma)) ...
                        + Vmax * 1/2 * ( 1+ erf((grid.x - 50 + lambda_is + lambda(2))/sigma));

V = @(lambda) min(interp1(grid.x-(lambda_is-initial_shift),V_comp_eta(lambda),grid.x,'linear','extrap'),interp1(grid.x+(lambda_is-initial_shift),V_comp_eta(-flip(lambda)),grid.x,'linear','extrap'));
dV_dlambda_fun = @(eta,eta_indx) 1/(1e-3)*(V(eta+[zeros(eta_indx-1,1);5e-4;zeros(numel(eta)-eta_indx,1)])-V(eta-[zeros(eta_indx-1,1);5e-4;zeros(numel(eta)-eta_indx,1)]));

