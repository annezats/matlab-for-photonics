% tmatrix for propagation through a film of thickness d
% propagation direction is x
function pmat = tmat_phase(d,Ni,k0,k0z)
    kix = ((Ni*k0)^2 - k0z^2)^0.5;
    
    pmat(1) = exp(-kix*d*1i);
    pmat(2) = exp(kix*d*1i);
end
