% input- lorentzian oscilators centered at: 1eV, 1.5eV, 2eV
wo=1
%linewidth=100meV =damping factor
gam=0.1
%plasma frequency(wp) = 10eV (germanium-intrinsic semiconductor)
wp=10
%eoo=1 usually
%x is omega (units of eV)
x=0:0.01:2

%step 1: find dielectric from lorenzian
%ereal=(wp^2).*(wo^2.-x^2)./((wo^2-x^2).^2+x^2.*gam^2)
%eimag=(wp^2*x*gam)./((wo^2-x^2).^2+x^2.*gam^2)
%e=ereal+i*eimag
e = (wp^2)./((wo^2 - x.^2) - 1i.*x.*gam);
etot=1+e
ereal=real(etot)
eimag=imag(etot)

figure

%real and imaginary dielectric function
subplot(1,4,1)
plot(x,ereal,x,eimag)

%real and imaginary refractive index
subplot(1,4,2)
n=etot.^0.5
nreal=real(n)
nimag=imag(n)
plot(x,nreal,x,nimag)

%absorption coefficient
freq=4*10^-14%not sure what to put here
subplot(1,4,3)
alpha=4*pi*nimag*freq
plot(x,alpha)

%imaginary dielectric+imaginary refractive+absorption coefficient
subplot(1,4,4)
plot(x,eimag,x,nimag,x,alpha)