%DtN-TDG solver for Helmholtz equation on periodic grating

close all; addpath src

%-----------------------------------
%Parameters definition
%-----------------------------------
%Problem parameters
param.K=5; %wavenumber
param.theta=-pi/3; %incident angle
param.alp=param.K*cos(param.theta); %quasi-periodicity parameter

%Discretization parameters
param.h=1.5; %mesh width
param.alpha=1/2; param.beta=1/2; param.delta=1/2; %TDG flux coefficients 
param.M=20; %number of Fourier modes

%-----------------------------------
%Mesh definition
%-----------------------------------
%Select domain
domain = 'double_rectangle_abs';
[mesh,param] = GenerateMesh(param,domain);

%-----------------------------------
%Numerical solution with DtN-TDG 
%-----------------------------------
param.nd=24; %number of plane wave directions
%define the plane wave direction vectors
param.d=zeros(2,param.nd); 
for l=1:param.nd
    param.d(:,l)=[cos((2*pi*l)/param.nd); sin((2*pi*l)/param.nd)];
end

disp(['DtN-TDG approximation on the domain ', domain, ' with k=', num2str(param.K),...
    ', p=', num2str(param.nd), ', h=', num2str(param.h),...
    ', with ', num2str(param.nd.*size(mesh.t,1)), ' basis functions' ])

if param.K.*param.L/2 > param.M
    warning('The value of M is too small, consider increasing for better stability')
end

%basis functions
phi = @(x1,x2,d,k) exp(1i*k.*(x1.*d(1)+x2.*d(2)));

A = MatrixDtNTDG(mesh,param); %system matrix
b = rhsDtNTDG(mesh,param); %system rhs
u=A\b; %solve the system
disp('Solved linear system')


%-----------------------------------
%Solution plot
%-----------------------------------
PlotSolution(mesh,param,u,phi); %plots real part and absolute value
