%% defining sample
clear
%close all

% import materials library
load('lotsofdata.mat')

% initial layer
% material name has to be in curly brackets to define it as a "cell"
% then can use function: size(cellfun('length', sample(n).material),2) 
% to determine number of elements in the cell
sample(1).material = {'air'};
sample(1).repeat = 1;
sample(1).lambda = 300*10^-9:100*10^-9:20000*10^-9;
sample(1).N = (1.00027717 + 0*1i).*ones(1,size(sample(1).lambda,2));

n=12;
m=n;
sample(2).material = {'n1', 'n2'};
sample(2).thickness = {110E-9, 80E-9};
sample(2).repeat = m;
sample(2).lambda = {sample(1).lambda, sample(1).lambda};
sample(2).N = {(1.5).*ones(1,size(sample(1).lambda,2)),(2).*ones(1,size(sample(1).lambda,2))};

sample(3).material = {'PMMA'};
sample(3).thickness = 160E-9;
sample(3).repeat = 1;
sample(3).lambda = (M.PMMA.wavelength)*10^-9;
sample(3).N = M.PMMA.n + (M.PMMA.k).*1i;

sample(4).material = sample(2).material;
sample(4).thickness = sample(2).thickness;
sample(4).repeat = n;
sample(4).lambda = sample(2).lambda;
sample(4).N = sample(2).N;

% substrate layer
sample(5).material = {'glass'};
sample(5).repeat = 1;
sample(5).lambda = (M.glass.wavelength)*10^-9;
sample(5).N = M.glass.n + (M.glass.k).*1i;

%% main
clear 'Tr_P'; clear 'Tr_S'; clear 'Re_P'; clear 'Re_S'
% set up lambda and angle variables
% incident angle
%theta0degrees_a = 0:5:50;
theta0degrees_a = 0;
%theta0degrees_a = linspace(10, 20, 5);
theta0_a = theta0degrees_a * (pi/180);

% wavelengths to calculate
lambda_a = linspace(400E-9, 800E-9, 2000);
%lambda_a = 500E-9;

% make sure sample file has already been constructed
% using thin_film_build.m

n_entries = size(sample,2);
i=1;
n_layers=0;

% determine the total number of layers 
while(i <= n_entries)
    n_layers = n_layers + 1*size(cellfun('length', sample(i).material),2)*sample(i).repeat;
    i=i+1;
end

% construct new layer file that exapands out sample... make equivalent to
% case code 

% write first and last layers
%layer(1).N = sample(1).N; 
layer(1).N =interp1(sample(1).lambda, sample(1).N, lambda_a);

%layer(n_layers).N = sample(n_entries).N;
layer(n_layers).N =interp1(sample(n_entries).lambda, sample(n_entries).N,lambda_a);

% write repeating layers
i=2;
% dummy layer index
iln = i;
while(i < n_entries)
    if size(cellfun('length', sample(i).material),2) == 1
        %layer(i).N = sample(i).N;
        layer(iln).N =interp1(sample(i).lambda, sample(i).N,lambda_a);
        layer(iln).thickness = sample(i).thickness;
        iln = iln + 1;
    else
        count = size(cellfun('length', sample(i).material),2);
        repeat = sample(i).repeat;
        k=0;
        while(k<count*repeat)
            j=0;
            while(j < count)
                %layer(i+j+k).N = sample(i).N{1,j+1};
                layer(iln+j+k).N =interp1(sample(i).lambda{1,j+1}, sample(i).N{1,j+1},lambda_a);
                %disp (iln+j+k)
                layer(iln+j+k).thickness = sample(i).thickness{1,j+1};
                j=j+1;
            end
            k=k+j;
        end
        iln = iln+k;
    end
    i=i+1;
end
%return
% loop through angles and wavelengths
i = 1;
j = 1;

imax = size(lambda_a,2);
jmax = size(theta0_a,2);

while(j<=jmax)
    i = 1;
    while(i<=imax)
        lambda = lambda_a(i);
        theta0 = theta0_a(j);
        
        % initial wavevector
        k0 = 2*pi/lambda;
        
        % component along surface
        k0z = k0*layer(1).N(i)*sin(theta0);
        
        % get array of fresnel coefficients
        % work from back to front: save individual and multiply out
        %fmvarname = genvarname((['FM', num2str(n_entries), num2str(n_entries-1)]));
        %pmvarname = genvarname((['PM', num2str(n_entries), num2str(n_entries-1)]));
        
        lnum = n_layers;
%         TtotS = 0*ones(2);
%         TtotP = 0*ones(2);
        
        % calculate last layer (no propagation)
        N0 = layer(lnum-1).N(i);
        N1 = layer(lnum).N(i);
        fm = tmat_fresnel(N0,N1,k0,k0z);
        % S polarization
        T01S = (1/fm(2))*[1 fm(1); fm(1) 1]; 
        % P polarization
        T01P = (1/fm(4))*[1 fm(3); fm(3) 1];
        
        lnum = lnum - 1;
        while(lnum > 1)
            N1 = layer(lnum).N(i);
            N0 = layer(lnum-1).N(i);
            d1 = layer(lnum).thickness;
            % propagation in this layer
            pm = tmat_phase(d1,N1,k0,k0z);
            T01S = ([pm(1) 0; 0 pm(2)]) * T01S;
            T01P = ([pm(1) 0; 0 pm(2)]) * T01P;
            % fresnel coeffs for transfer into this layer
            fm = tmat_fresnel(N0,N1,k0,k0z);
            T01S = ((1/fm(2))*[1 fm(1); fm(1) 1])*T01S;
            T01P = ((1/fm(4))*[1 fm(3); fm(3) 1])*T01P;
            lnum = lnum - 1;
        end
        
        % starting medium - one more sett of fresnel coeffs
%         N0 = layer(1).N(i);
%         N1 = layer(2).N(i);
%         fm = tmat_fresnel(N0,N1,k0,k0z);
%         % S polarization
%         T01S = ((1/fm(2))*[1 fm(1); fm(1) 1])*T01S; 
%         % P polarization
%         T01P = ((1/fm(4))*[1 fm(3); fm(3) 1])*T01P;
           % E field amplitude ratios for forward and backward waves
        rS(i) = T01S(2,1)/T01S(1,1);
        tS(i) = 1/T01S(1,1);
    
        rP(i) = T01P(2,1)/T01P(1,1);
        tP(i) = 1/T01P(1,1);

        % kx, final
        Nf = layer(n_layers).N(i);
        kNx = ((Nf*k0)^2 - k0z^2)^0.5;
        % kz, final = k0z = k0*sin(theta0);

        % theta final
        thetaN = real(atan(k0z/kNx));
        thetaNangle = thetaN * (180/pi);
    
        % Power Reflection Coefficient
        Re_S(i,j) = real(rS(i)*conj(rS(i)));
        Re_P(i,j) = real(rP(i)*conj(rP(i)));

        % Power Transmission Coefficient
        nreal_f = real(layer(n_layers).N(i));
        nreal_i = real(layer(1).N(i));
        Tr_S(i,j) = real((nreal_f/nreal_i)*(cos(thetaN)/cos(theta0))*tS(i)*conj(tS(i)));
        Tr_P(i,j) = real((nreal_f/nreal_i)*(cos(thetaN)/cos(theta0))*tP(i)*conj(tP(i)));
    
        % Power Absorption Coefficient
        Abs_S(i,j) = 1 - Re_S(i,j) - Tr_S(i,j);
        Abs_P(i,j) = 1 - Re_P(i,j) - Tr_P(i,j);
    
    % total check
%     Tot_S(i,j) = Re_S(i,j) + Tr_S(i,j) + Abs_S(i,j);
%     Tot_P(i,j) = Re_P(i,j) + Tr_P(i,j) + Abs_P(i,j);
        i=i+1;
    end
    j=j+1;
end

 clear i; clear j; clear k;
 %clear n_entries; clear n_layers;
 clear imax; clear jmax;
 clear d1; clear fm; 
 clear k0; clear k0z; clear kNx;
 clear lambda; clear lnum; clear N0; clear N1; clear Nf; 
 clear nreal_f; clear nreal_i;
 clear pm; clear rP; clear rS; clear tP; clear tS;
 clear T01P; clear T01S;
 
 %% plot
 % display
% 
% subplot(2,1,1), plot(lambda_a,Re_S,lambda_a,Re_P), title('Reflection Coeff'), legend('S','P')
% ylim([0,1])
% set(gca,'fontsize', 16);
% subplot(2,1,2), plot(lambda_a,Tr_S,lambda_a,Tr_P), title('Transmission Coeff'), legend('S','P')
% ylim([0,1])
% set(gca,'fontsize', 16);

figure
subplot(3,1,1), plot(lambda_a,Re_S,lambda_a,Re_P), title('Reflection Coeff'), legend('S','P')
ylim([0,1])
set(gca,'fontsize', 14);
subplot(3,1,2), plot(lambda_a,Tr_S,lambda_a,Tr_P), title('Transmission Coeff'), legend('S','P')
ylim([0,1])
set(gca,'fontsize', 14);
subplot(3,1,3), plot(lambda_a,Abs_S,lambda_a,Abs_P), title('Absorption Coeff'), legend('S','P')
ylim([0,1])
set(gca,'fontsize', 14);
sgtitle('4')

% lambda_eV = 1240./(lambda_a/10^-9);
% imagesc(theta0degrees_a,lambda_eV, Re_P)
% axis xy
% xlabel('angle')
% ylabel('energy')
% ylim([1.5,3.5])
% set(gca,'fontsize', 16);

%% q factor 

%resonant freq= x when y=0, or ymin
six=1000;%index of lambda=6
seven=1700;%index of lambda=7
[FWHM,x,xin]=FWHMfunc(lambda_a,Tr_P,six,seven);
Q=x/FWHM

%% photon lifetime:
c0=3*10^8;
v=c0./lambda_a;
[vFWHM,vx,vxin]=FWHMfunc(v,Tr_P,six,seven);%dont make a func w the same name as a var
vFWHM=-1*vFWHM
tp=1/(2*pi*vFWHM)%tp=1/(fsr*(1-R^2))%1/(2*pi*fsr)%=Q/(2*pi*c)%photon lifetime


%% round trips
%q=2*(sample(3).thickness)/x; %integer...is not an integer...the relation between res. freq and thickness
% F=Q/q; %Q=qF Finesse
fsr=vFWHM;%fsr=F/FWHM; % F=fsr/FWHM
% n=1;
% c=c0/n %speed of light in cavity= speed of light in vacuum/n
c=v(vxin);
n=c0/c;
theta=0;
fsr=c/(2*n*(sample(3).thickness)*cos(theta)); %full spectral range
%THINGS I DO NOT KNOW:
% c,speed of light in cavity is this not frequency at the cavity mode?
% n of cavity??
%R1 and R2 of the cavity=the reflection of each DBR=R^2
%R=0.02216 %from hw 5 for 4=0.03097 6=0.02216
rt=1/fsr; %1 round trip
trips=tp/rt%=1/2*pi?? %number of round trips =1/2pi*FWHM1^2