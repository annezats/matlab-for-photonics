clear
clf

% use units of eV
eV = 0:.01:3;

% plasma frequency
wp = 0.2;

% define center and width for first lorentzian function
w1 = 0.5;
w2 = 1.0;
w3= 1.5;
g1 = 0.1;%from 1

% first function
l1 = (wp^2)/((w1^2 - eV.^2) - 1i.*eV.*g1);
l2 = (wp^2)/((w2^2 - eV.^2) - 1i.*eV.*g1);
l3 = (wp^2)/((w3^2 - eV.^2) - 1i.*eV.*g1);
% total dielectric
df = 1 + l1;

% real part
e1 = real(df);

% imaginary part
e2 = imag(df);

%plot
plot(eV, e1, eV, e2)

