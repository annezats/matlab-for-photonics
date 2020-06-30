% calculate thin film transmission and reflection on an infinite substrate
% x is propagation direction
% z is parallel to surface

clear
close all

% incident angle
theta0degrees = 10;
theta0 = theta0degrees * (pi/180);

% import material coefficients here
%rfi1 = readmatrix('glass-ri.csv');
rfi1 = readmatrix('silicon.csv');
%rfi1 = readmatrix('silver.csv');

% interpolate NaNs
rfi1_corr = fillmissing(rfi1, 'linear');

% wavelengths
lambda1_a = rfi1_corr(:,1)*10^-6;
% real part of refractive index
real1_a = rfi1_corr(:,2);
% real part of refractive index
imag1_a = rfi1_corr(:,3);

% calculation absorption coeff to compare to R,T
% alpha_a = (4*pi.*imag1_a)./lambda1_a;
% calculate dielectric function
% eps_tot_a = (real1_a + (imag1_a.*1i)).^2;
% eps_real_a = real(eps_tot_a);
% eps_imag_a = imag(eps_tot_a);

% air
N0 = 1.00027717 + 0*1i;

% loop through length of array
lindex = numel(lambda1_a);
l = 1; % counting index

while (l <= lindex)  % wavelength
    lambda = lambda1_a(l);
    N1 = real1_a(l) + imag1_a(l)*1i;
    
    % initial wavevector
    % vacuum wavevector
    k0 = 2*pi/lambda;

    % component along surface
    k0z = k0*N0*sin(theta0);

    % get array of fresnel coefficients
    % rijS, tijS, rijP, tijS
    f01 = fresnel(N0,N1,k0,k0z);

    % transfer matrix from 0 to 1
    % S polarization
    T01S = (1/f01(2))*[1 f01(1); f01(1) 1]; 

    % P polarization
    T01P = (1/f01(4))*[1 f01(3); f01(3) 1]; 

    % E field amplitude ratios for forward and backward waves
    rS = T01S(2,1)/T01S(1,1);
    tS = 1/T01S(1,1);
    
    rP = T01P(2,1)/T01P(1,1);
    tP = 1/T01P(1,1);

    % kx, final
    kNx = ((N1*k0)^2 - k0z^2)^0.5;
    % kz, final = k0z = k0*sin(theta0);

    % theta final
    thetaN = real(atan(k0z/kNx));
    thetaNangle = thetaN * (180/pi);
    
    % Power Reflection Coefficient
    Re_S(l) = rS*conj(rS);
    Re_P(l) = rP*conj(rP);

    % Power Transmission Coefficient
    Tr_S(l) = real((real(N1)/real(N0))*(cos(thetaN)/cos(theta0))*tS*conj(tS));
    Tr_P(l) = real((real(N1)/real(N0))*(cos(thetaN)/cos(theta0))*tP*conj(tP));
    
    % Absorption Transmission Coefficient
    %Abs_S(l) = 1 - Re_S(l) - Tr_S(l);
    %Abs_P(l) = 1 - Re_P(l) - Tr_P(l);
    
    % total check
    Tot_S(l) = Re_S(l) + Tr_S(l);
    Tot_P(l) = Re_P(l) + Tr_P(l);
    % increase index
    l= l+1; 
end

subplot(2,1,1), plot(lambda1_a,Re_S,lambda1_a,Re_P), title('Reflection Coeff'), legend('S','P')
subplot(2,1,2), plot(lambda1_a,Tr_S,lambda1_a,Tr_P), title('Transmission Coeff'), legend('S','P')
%subplot(3,1,3), plot(lambda1_a,Abs_S,lambda1_a,Abs_P), title('Absorption Coeff'), legend('S','P')
xlabel('Wavelength') 

%figure()
%plot(lambda1_a,Re_S,lambda1_a,Tr_S,lambda1_a,alpha_a/3e8,'k'),, legend('R','T','alpha')

% figure()
% %plot(lambda1_a,zeros(size(lambda1_a)), 'k'); hold on
% plot(lambda1_a,Re_S,lambda1_a,Tr_S,lambda1_a,eps_real_a), legend('R','T','eps_real'); hold off
% xlim([0.2e-6 0.5e-6])
% ylim([-0.2 1.1])

