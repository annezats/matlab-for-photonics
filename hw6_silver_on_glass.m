
% calculate thin film transmission and reflection on an infinite substrate
% x is propagation direction
% z is parallel to surface

clear
close all

% incident angle
theta = 10;
ti = theta * (pi/180);

% import material coefficients here

% air
N0 = 1.00027717 + 0*1i;

silver = readmatrix('silver.csv');
glass = readmatrix('BK7glass.csv');
%material = readmatrix('silicon.csv');%material=importdata("silicon.csv")
%rfi1 = readmatrix('silver.csv');

% interpolate NaNs
silver = fillmissing(silver, 'linear');
glass = fillmissing(glass, 'linear');
%cut off last few rows that we dont need
glass = glass(1:85,:);

% wavelengths
wll = glass(:,1)*10^-6; %wl=material.data(:,1)
% real part of refractive index
nt1 = silver(:,2);
nt2 = glass(:,2);
% imaginary part of refractive index
kt1 =silver(:,3);
kt2 =glass(:,3);


% calculation absorption coeff to compare to R,T
% alpha_a = (4*pi.*imag1_a)./lambda1_a;
% calculate dielectric function
% eps_tot_a = (real1_a + (imag1_a.*1i)).^2;
% eps_real_a = real(eps_tot_a);
% eps_imag_a = imag(eps_tot_a);

% loop through length of array
lindex = numel(wll);
l = 1; % counting index

while (l <= lindex)  % wavelength
    wl = wll(l);
    N1 = nt1(round(l*1.66)) + kt1(round(l*1.66))*1i; %silver
    N2 = nt2(l) + kt2(l)*1i; %glass

    % initial wavevector 
    % vacuum wavevector 
    k0 = 2*pi/wl;

    % component along surface (k perpendicular)
    k0z = k0*N0*sin(ti);

    % get array of fresnel coefficients
    % rijS, tijS, rijP, tijS
    f01 = fresnel_from_prof(N0,N1,k0,k0z);

    % transfer matrix from 0 to 1
    % S polarization
    T01S = (1/f01(2))*[1 f01(1); f01(1) 1]; 

    % P polarization
    T01P = (1/f01(4))*[1 f01(3); f01(3) 1]; 
    
    %MY ADDITIONS-------------------------------------
    % kx1 (MOVED THIS UP)
    kx1 = ((N1*k0)^2 - k0z^2)^0.5;
    % kz, final = k0z = k0*sin(theta0);
    
    %d1 the thickness of material 1
    d1= 30*10^-9; 
    e=2.71828 ;
    
    %propogation matrix through material 1
    T1=[ e^(j*kx1*d1) 0; 0 e^(-j*kx1*d1)];
    
    % get array of fresnel coefficients
    % rijS, tijS, rijP, tijS
    f12 = fresnel_from_prof(N1,N2,k0,k0z); %pretty sure k0 stays the same?

    % transfer matrix from 1 to 2
    % S polarization
    T12S = (1/f12(2))*[1 f12(1); f12(1) 1]; 

    % P polarization
    T12P = (1/f12(4))*[1 f12(3); f12(3) 1]; 
    
    %Total transfer matrix
    % to add more layers just keep multiplying by *T1*T12S/p
    T02S=T01S*T1*T12S;
    T02P=T01P*T1*T12P;
    
     % kx2
    kx2 = ((N2*k0)^2 - k0z^2)^0.5; %need this for final theta
    % kz, final = k0z = k0*sin(theta0);
    
    %------------------------------------------
    % E field amplitude ratios for forward and backward waves
    rS = T02S(2,1)/T02S(1,1);
    tS = 1/T02S(1,1);
    
    rP = T02P(2,1)/T02P(1,1);
    tP = 1/T02P(1,1);

    % theta final
    %thetaN = real(atan(k0z/kx1));
    thetaN = real(atan(k0z/kx2));%doesnt depend on intermediates
    thetaNangle = thetaN * (180/pi);
    
    % Power Reflection Coefficient
    Re_S(l) = rS*conj(rS);
    Re_P(l) = rP*conj(rP);
    R(l)=(Re_S(l)+Re_P(l))/2;

    % Power Transmission Coefficient
    Tr_S(l) = real((real(N2)/real(N0))*(cos(thetaN)/cos(ti))*tS*conj(tS));
    Tr_P(l) = real((real(N2)/real(N0))*(cos(thetaN)/cos(ti))*tP*conj(tP));
    T(l)=(Tr_S(l)+Tr_P(l))/2;
    
    % Absorption Transmission Coefficient
    Abs_S(l) = 1 - Re_S(l) - Tr_S(l);
    Abs_P(l) = 1 - Re_P(l) - Tr_P(l);
    A(l)=1- R(l)- T(l)%(Abs_S(l)+Abs_P(l))/2;
    
    % total check
    Tot_S(l) = Re_S(l) + Tr_S(l);
    Tot_P(l) = Re_P(l) + Tr_P(l)
    % increase index
    l= l+1; 
end

%subplot(2,1,1), plot(wll,Re_S,wll,Re_P), title('Reflection Coeff'), legend('S','P')
%subplot(2,1,2), plot(wll,Tr_S,wll,Tr_P), title('Transmission Coeff'), legend('S','P')
plot(wll,R,wll,T,wll,A), title('A = 1-R-T'), legend('R','T','A');
xlabel('Wavelength') 

%figure()
%plot(lambda1_a,Re_S,lambda1_a,Tr_S,lambda1_a,alpha_a/3e8,'k'),, legend('R','T','alpha')

% figure()
% %plot(lambda1_a,zeros(size(lambda1_a)), 'k'); hold on
% plot(lambda1_a,Re_S,lambda1_a,Tr_S,lambda1_a,eps_real_a), legend('R','T','eps_real'); hold off
% xlim([0.2e-6 0.5e-6])
% ylim([-0.2 1.1])

