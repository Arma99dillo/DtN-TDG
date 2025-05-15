%DtN-TDG solver for Helmholtz equation on periodic grating
clearvars;

%-----------------------------------
%Parameters definition
%-----------------------------------
%Problem parameters
param.K=20; %wavenumber
param.theta=-pi/3; %incident angle
param.alp=param.K*cos(param.theta); %quasi-periodicity parameter

%Discretization parameters
param.h=0.5; %mesh width
param.alpha=1/2; param.beta=1/2; param.delta=1/2; %TDG flux coefficients 
param.M=50; %number of Fourier modes


%-----------------------------------
%Mesh definition
%-----------------------------------
%Select domain
domain = 'u_shape';
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

%basis functions
phi = @(x1,x2,d,k) exp(1i*k.*(x1.*d(1)+x2.*d(2)));

A = MatrixDtNTDG(mesh,param); %system matrix
b = rhsDtNTDG(mesh,param); %system rhs
u=A\b; %solve the system


%-----------------------------------
%Solution plot
%-----------------------------------
PlotSolution(mesh,param,u,phi); %plots real part and absolute value
