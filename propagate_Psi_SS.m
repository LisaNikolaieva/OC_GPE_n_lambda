function [Psi] = propagate_Psi_SS(dt,grid,par,Psi,V,mu)

g = par.g;
m = par.m;

kin_en  = ( - 0.5 * grid.ilap / m ); 
eD_hp      = exp( - 1i * kin_en * dt /2 );

[Psi] = sstep(dt,eD_hp,Psi,mu,g,V);

end

function [Psi] = sstep(dt,eD_hp,Psi,mu,g,V)

    Psi = ifftn(eD_hp.*fftn(Psi));
    Psi = exp((V-mu+g*(abs(Psi)).^2)*(dt/1i)).*Psi;
    Psi = ifftn(eD_hp.*fftn(Psi));  

end