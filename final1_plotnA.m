clear
% angle of incidence (arbitrary constant)
ti=0; %depending on the angle i'm getting different results
Ni= 1.0003 %air is a constant
%-------------------------
material= importdata("n and k of RP and DJ.xlsx", ",",1);
% refractive index of material
%wavelength
ev=material.data.DJ(:,1);
wl= 1240./ev; % := all rows, 1= which column
%wl=1240/eV;
for i=2 %i+1 is which m value to use
    nt=material.data.DJ(:,1+i);
    kt=material.data.DJ(:,4+i);
    1j
    Nt=nt+1j*kt;
    A_coef=4*pi*kt./wl
   
    titl=append('m=',sprintf('%0.5f', i+1))
    figure
    subplot(2,1,1), plot(wl,nt,wl,kt), legend('n','k');
    subplot(2,1,2), plot(wl,A_coef), legend('A'); %blue is the first one, red 2nd
   %all is reflected or absorbed, none reflected
    %fwhm=FWHMfunc(wl,kt,132,171)
end


clear i; clear ans; clear titl; 
clear Rp; clear Rs; 
clear Tp; clear Ts;

