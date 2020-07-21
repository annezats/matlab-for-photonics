clear
% angle of incidence (arbitrary constant)
ti=pi/5 %depending on the angle i'm getting different results
Ni= 1.0003 %air is a constant
%-------------------------
material= importdata("silicon.csv", ",",1);
materialgood= myinterpolate2(material);
%switch name to silicon.csv
%or to silver.csv
% refractive index of material
nt=materialgood(:,2);
kt=materialgood(:,3);
1j
Nt=nt+1j*kt;
%-------------------------
%calculation maagic (see the other file)
[R,T,Rp,Rs,Tp,Ts]=myRT(ti,Ni,Nt);
%R=1-T
%wavelength
wl= material.data(:,1); % := all rows, 1= which column
%the plot thickens..
figure
plot(wl,R,wl, T); %blue is the first one, red 2nd



