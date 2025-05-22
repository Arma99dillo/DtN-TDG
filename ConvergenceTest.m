% DtN-TDG solver for Helmholtz equation on periodic grating 
% Convergence test

close all; addpath quadtriangle\; addpath src\

%-----------------------------------
%Parameters definition
%-----------------------------------
%Problem parameters
param.K=5; %wavenumber
param.theta=-pi/4; %incident angle
param.alp=param.K*cos(param.theta); %quasi-periodicity parameter

%Discretization parameters
param.h=1.5; %mesh width
param.alpha=1/2; param.beta=1/2; param.delta=1/2; %TDG flux coefficients 
param.M=20; %number of Fourier modes


%-----------------------------------
%Mesh definition
%-----------------------------------
%Select domain - only double and triple rectangle available
domain = 'double_rectangle';
%generate mesh and compute exact solution
[mesh,param,uex,uexdx,uexdy] = GenerateMeshSol(param,domain);

disp(['DtN-TDG convergence test on the domain ', domain, ' with k=', num2str(param.K)])

if param.K.*param.L/2 > param.M
    warning('The value of M is too small, consider increasing for better stability')
end

%-----------------------------------
%Cycle on number of directions
%-----------------------------------
%basis functions and derivatives
phi = @(x1,x2,d,k) exp(1i*k.*(x1.*d(1)+x2.*d(2)));
grad_phi = @(x1,x2,d,k) 1i*k.*d.*exp(1i*k.*(x1.*d(1)+x2.*d(2)));

ni=3; nf=30; %min and max value
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
    [err2,err1] = SolErr(mesh,param,u,phi,grad_phi,uex,uexdx,uexdy);
    L2Error(v) = err2;
    H1Error(v) = err1;
    disp([ 'Computed error for p=', num2str(param.nd) ])


    v=v+1;
end


%-----------------------------------
%Convergence plot
%-----------------------------------
PlotConvergence(ni,nf,L2Error,H1Error);


%-----------------------------------
%Solution and error plot
%-----------------------------------
PlotErrorLog(mesh,param,u,phi,uex);
