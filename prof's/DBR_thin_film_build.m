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

% single intermediate layer example
% sample(2).material = {'Ag'};
% sample(2).thickness = 100E-9;
% sample(2).repeat = 1;
% sample(2).lambda = (M.Ag.wavelength)*10^-9;
% sample(2).N = M.Ag.n + (M.Ag.k).*1i;


% repeating intermediate layer example
% sample(2).material = {'Ag', 'PMMA'};
% sample(2).thickness = {100E-9, 100E-9};
% sample(2).repeat = 3;
% sample(2).lambda = {(M.Ag.wavelength)*10^-9, (M.PMMA.wavelength)*10^-9};
% sample(2).N = {M.Ag.n + (M.Ag.k).*1i, M.PMMA.n + (M.PMMA.k).*1i};

sample(2).material = {'SiO2', 'Si3N4'};
sample(2).thickness = {110E-9, 80E-9};
sample(2).repeat = 6;
sample(2).lambda = {(M.SiO2.wavelength)*10^-9, (M.Si3N4.wavelength)*10^-9};
%sample(2).N = {M.SiO2.n + (M.SiO2.k).*1i, M.Si3N4.n + (M.Si3N4.k).*1i};
sample(2).N = {M.SiO2.n, M.Si3N4.n};

sample(3).material = {'PMMA'};
sample(3).thickness = 160E-9;
sample(3).repeat = 1;
sample(3).lambda = (M.PMMA.wavelength)*10^-9;
sample(3).N = M.PMMA.n + (M.PMMA.k).*1i;

% sample(4).material = {'Si3N4','SiO2'};
% sample(4).thickness = {80E-9,110E-9};
% sample(4).repeat = 6;
% sample(4).lambda = {(M.Si3N4.wavelength)*10^-9,(M.SiO2.wavelength)*10^-9};
% %sample(4).N = {M.SiO2.n + (M.SiO2.k).*1i, M.Si3N4.n + (M.Si3N4.k).*1i};
% sample(4).N = {M.Si3N4.n, M.SiO2.n};

sample(4).material = {'SiO2', 'Si3N4'};
sample(4).thickness = {110E-9, 80E-9};
sample(4).repeat = 6;
sample(4).lambda = {(M.SiO2.wavelength)*10^-9, (M.Si3N4.wavelength)*10^-9};
%sample(4).N = {M.SiO2.n + (M.SiO2.k).*1i, M.Si3N4.n + (M.Si3N4.k).*1i};
sample(4).N = {M.SiO2.n, M.Si3N4.n};

% substrate layer
sample(5).material = {'glass'};
sample(5).repeat = 1;
sample(5).lambda = (M.glass.wavelength)*10^-9;
sample(5).N = M.glass.n + (M.glass.k).*1i;


