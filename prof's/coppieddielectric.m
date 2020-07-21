%------------------------------------------------------------------------%
%----------------- complex dielectric function of a ---------------------%
%------------------------ Lorentz medium --------------------------------%
%------------------------------------------------------------------------%
clc
close all
clear all
%------------------------------------------------------------------------%
hw0   =  1;      % resonant frequency
gam   =  10e-3;    % damping frequency
f     =  8.8e-3;   % oscillator strength
eb    =  11.2;     % background epsilon
%------------------------------------------------------------------------%
om    =  1:0.001:1.6;  % frequency array
ep    =  eb + (  ( f * hw0*hw0 )./ (hw0*hw0-om.*om-1i.*om*gam)); 
%------------------------------------------------------------------------%
figure(1)
subplot(211)
plot ( om, real(ep),'linewidth',2)
xlabel('Energy (eV)','fontsize',14)
ylabel('real(\epsilon(\omega))','fontsize',14)
title('Dielectric function of a Lorentz medium','fontsize',14)
line([hw0 hw0], [10 12],'color','k','Linestyle','--','linewidth',2)
axis([1.0 1.6 10 12])

subplot(212)
plot ( om, imag(ep),'r','linewidth',2)
xlabel('Energy (eV)','fontsize',14)
ylabel('im(\epsilon(\omega))','fontsize',14)
line([hw0 hw0], [0 2],'color','k','Linestyle','--','linewidth',2)
axis([1.0 1.6 0 2])
%------------------------------------------------------------------------%