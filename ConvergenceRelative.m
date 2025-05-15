% DtN-TDG solver for Helmholtz equation on periodic grating 
% Relative error convergence test

clearvars; addpath quadtriangle\

%-----------------------------------
%Parameters definition
%-----------------------------------
%Problem parameters
param.K=5; %wavenumber
param.theta=-pi/4; %incident angle
param.alp=param.K*cos(param.theta); %quasi-periodicity parameter

%Discretization parameters
param.h=0.75; %mesh width
param.alpha=1/2; param.beta=1/2; param.delta=1/2; %TDG flux coefficients 
param.M=100; %number of Fourier modes
nd_raff=20; %number of directions for the refined solution


%-----------------------------------
%Mesh definition
%-----------------------------------
domain = 'dir_square'; %select domain
[mesh,param] = GenerateMesh(param,domain); %generate mesh 

%%
%-----------------------------------
%Refined solution for relative error
%-----------------------------------
%basis functions and derivatives
phi = @(x1,x2,d,k) exp(1i*k.*(x1.*d(1)+x2.*d(2)));
grad_phi = @(x1,x2,d,k) 1i*k.*d.*exp(1i*k.*(x1.*d(1)+x2.*d(2)));

param_raff=param;
param_raff.nd=nd_raff;

%define the plane wave direction vectors
param_raff.d=zeros(2,param_raff.nd); 
for l=1:param_raff.nd
    param_raff.d(:,l)=[cos((2*pi*l)/param_raff.nd); sin((2*pi*l)/param_raff.nd)];
end

%matrix, rhs and solution for the refined parameters
A_raff = MatrixDtNTDG(mesh,param_raff);
b_raff = rhsDtNTDG(mesh,param_raff);
u_raff=A_raff\b_raff;

%%
%-----------------------------------
%Cycle on number of directions
%-----------------------------------
ni=3; nf=nd_raff-1; %min and max value
L2Error=zeros(nf-ni+1,1); H1Error=zeros(nf-ni+1,1); v=1; %error vectors
for nd=ni:nf %cycle on number of PW directions

    param.nd=nd; %number of plane wave directions
    %define the plane wave direction vectors
    param.d=zeros(2,param.nd); 
    for l=1:param.nd
        param.d(:,l)=[cos((2*pi*l)/param.nd); sin((2*pi*l)/param.nd)];
    end
    
    A = MatrixDtNTDG(mesh,param); %system matrix
    b = rhsDtNTDG(mesh,param); %system rhs
    u=A\b; %solve the system

    %Error computation
    [err2,err1] = SolErrRel(mesh,param,u,param_raff,u_raff,phi,grad_phi);
    L2Error(v) = err2;
    H1Error(v) = err1;

    v=v+1;
end

%%
%-----------------------------------
%Convergence plot
%-----------------------------------
PlotConvergence(ni,nf,L2Error,H1Error);

%%
%-----------------------------------
%Solution and scatterer plot
%-----------------------------------
PlotErrRelLog(mesh,param,u,param_raff,u_raff,phi);
PlotScatterer(param,domain)
