function Psi = itp(V,par,grid)

mue_0 = 0;
n = 3000;
Psi = ones(size(grid.x));
normalize = @(Psi)Psi/sqrt(trapz(grid.x,conj(Psi).*Psi));

for i = 1:n
        Psi_next = propagate_Psi_SS(-1i*0.01,grid,par,Psi,V,mue_0);
    Psi_next = normalize(Psi_next);
    
    if trapz(grid.x,abs(Psi_next - Psi))<1e-9
        break;
    end
    Psi = Psi_next;
end
end

