%% defining sample
%clear
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

% METAL
sample(2).thickness = 17E-9; %10 is the highest, 17 is most centered on 590
sample(2).repeat = 1;
sample(2).material = {'Ag'};
sample(2).lambda = M.Ag.wavelength*10^-9;
sample(2).N = M.Ag.n + (M.Ag.k).*1i;

sample(3).material = {'2D DJ PEROVSKITE'};
%thickness is in the loop
sample(3).repeat = 1;
sample(3).lambda = wl*10^-9;
sample(3).N = nt+i*kt; %nt+i*kt;

% single intermediate layer example
sample(4).material = sample(2).material;
sample(4).thickness = 100E-9; %less and T starts, more no change
sample(4).repeat = 1;
sample(4).lambda = sample(2).lambda;
sample(4).N = sample(2).N;

% % substrate layer
sample(5).material = {'glass'};
sample(5).repeat = 1;
sample(5).lambda = (M.glass.wavelength)*10^-9;
sample(5).N = M.glass.n + (M.glass.k).*1i;

%% main


clear 'Tr_P'; clear 'Tr_S'; clear 'Re_P'; clear 'Re_S'; clear 'R';    % set up lambda and angle variables
% incident angle
%theta0degrees_a = 0:5:50;
theta0degrees_a = 20;
%theta0degrees_a = linspace(10, 20, 5);
theta0_a = theta0degrees_a * (pi/180);

% wavelengths to calculate
lambda_a = linspace(400E-9, 800E-9, 2000);
%lambda_a = 500E-9;

change=[0, 0.005,0.01,0.02]

for c=1:1:4
    
    cw=73*(1-change(c));  %for m=2 78,205,332,459;for m=3 75,199,323,447; for m=4 93, 234.5,376, 517.5
    sample(3).thickness = cw*10^-9;

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

    % write repeating layers
    i=1;
    while(i < n_entries+1)
        %layer(i).N = sample(i).N;
        layer(i).N =interp1(sample(i).lambda, sample(i).N,lambda_a);
        layer(i).thickness = sample(i).thickness;


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
    
%     finalinfo(c).thickness= cw*10^-9;
%     finalinfo(c).percent=change(c)*100;
%     finalinfo(c).R=Re_S;
%     
%     f=1;
%     if c>1
%         while f<=imax
%             R0=finalinfo(c-1).R(f);
%             Rprime=finalinfo(c).R(f);
%             dR(f)=(Rprime-R0);
%             f=f+1;
%         end
%         finalinfo(c).dR=dR;
%     end

     clear i; clear j; clear k;
     %clear n_entries; clear n_layers;
     clear imax; clear jmax;
     clear d1; clear fm; 
     clear k0; clear k0z; clear kNx;
     clear lambda; clear lnum; clear N0; clear N1; clear Nf; 
     clear nreal_f; clear nreal_i;
     clear pm; clear rP; clear rS; clear tP; clear tS;
     clear T01P; clear T01S;

 
end


