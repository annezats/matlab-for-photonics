% fresnel coefficients calculator
function fmat = tmat_fresnel(N0, N1, k0, k0z)
% calculate the fresnel coefficients for a plane wave

% calculate components normal to surface
%k0x = ((N0*k0)^2 - k0z^2)^0.5;
k0x = ((N0*k0)^2 - k0z^2)^0.5;
k1x = ((N1*k0)^2 - k0z^2)^0.5;

% TE polarization (S polarization)
% r12S
fmat(1) = (k0x-k1x)/(k0x+k1x);
% t12S
fmat(2) = (2*k0x)/(k0x+k1x);

% TM polarization (P polarization)
% r12P
fmat(3) = ((N0^2)*k1x-(N1^2)*k0x)/((N0^2)*k1x+(N1^2)*k0x);
% t12P
fmat(4) = (2*N0*N1*k0x)/((N0^2)*k1x+(N1^2)*k0x);

end


