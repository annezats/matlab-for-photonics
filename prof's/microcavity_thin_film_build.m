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
sample(2).thickness = 30E-9;
sample(2).repeat = 1;
sample(2).material = {'Ag'};
sample(2).lambda = M.Ag.wavelength*10^-9;
sample(2).N = M.Ag.n + (M.Ag.k).*1i;

% repeating intermediate layer example
% sample(2).material = {'Ag', 'PMMA'};
% sample(2).thickness = {100E-9, 100E-9};
% sample(2).repeat = 3;
% sample(2).lambda = {(M.Ag.wavelength)*10^-9, (M.PMMA.wavelength)*10^-9};
% sample(2).N = {M.Ag.n + (M.Ag.k).*1i, M.PMMA.n + (M.PMMA.k).*1i};

sample(3).material = {'PMMA'};
sample(3).thickness = 160E-9;
sample(3).repeat = 1;
sample(3).lambda = (M.PMMA.wavelength)*10^-9;
sample(3).N = M.PMMA.n + (M.PMMA.k).*1i;

% single intermediate layer example
sample(4).material = sample(2).material;
sample(4).thickness = 100E-9;
sample(4).repeat = 1;
sample(4).lambda = sample(2).lambda;
sample(4).N = sample(2).N;

% % substrate layer
sample(5).material = {'glass'};
sample(5).repeat = 1;
sample(5).lambda = (M.glass.wavelength)*10^-9;
sample(5).N = M.glass.n + (M.glass.k).*1i;


