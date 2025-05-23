% DtN-TDG solver for Helmholtz equation on periodic grating 
% Convergence test

close all; addpath quadtriangle; addpath src

%-----------------------------------
%Parameters definition
%-----------------------------------
%Problem parameters
param.K=5; %wavenumber
theta_0=-1.563806234657490;

%Discretization parameters
param.h=1.5; %mesh width
param.alpha=1/2; param.beta=1/2; param.delta=1/2; %TDG flux coefficients 
param.M=100; %number of Fourier modes

param.nd=30; %number of plane wave directions
%define the plane wave direction vectors
param.d=zeros(2,param.nd);
for l=1:param.nd
    param.d(:,l)=[cos((2*pi*l)/param.nd); sin((2*pi*l)/param.nd)];
end

%-----------------------------------
%Mesh definition
%-----------------------------------
%Select domain - only double and triple rectangle available
domain = 'triple_rectangle';

disp(['Theta-convergence test on the domain ', domain, ' with k=', num2str(param.K)...
    ', p=', num2str(param.nd), ', h=', num2str(param.h),...
    ', with ', num2str(param.nd.*size(mesh.t,1)), ' basis functions' ])

%-----------------------------------
%Cycle on theta
%-----------------------------------
%basis functions and derivatives
phi = @(x1,x2,d,k) exp(1i*k.*(x1.*d(1)+x2.*d(2)));
grad_phi = @(x1,x2,d,k) 1i*k.*d.*exp(1i*k.*(x1.*d(1)+x2.*d(2)));

v=1; %error vectors
for theta=theta_0-1e-4:1e-5:theta_0+1e-4 %cycle on theta
    
    param.theta=theta; %incident angle
    param.alp=param.K*cos(param.theta); %quasi-periodicity parameter

    %generate mesh and compute exact solution
    [mesh,param,uex,uexdx,uexdy] = GenerateMeshSol(param,domain);
    
    A = MatrixDtNTDG(mesh,param); %system matrix
    b = rhsDtNTDG(mesh,param); %system rhs
    u=A\b; %solve the system

    %Error computation
    [err2] = SolErr(mesh,param,u,phi,grad_phi,uex,uexdx,uexdy);
    L2Error(v) = err2;
    disp([ 'Computed error for incident angle theta=', num2str(param.theta) ])

    v=v+1;
end


%-----------------------------------
%Convergence plot
%-----------------------------------
figure()
semilogy(theta_0-1e-4:1e-5:theta_0+1e-4,L2Error,'*-','LineWidth',1.2); grid
xlim([theta_0-1.1e-4 theta_0+1.1e-4])
LL = legend('$L^2$ error','FontSize', 14);
set(LL, 'Interpreter', 'latex');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14,'TickLabelInterpreter', 'latex')
xticks([theta_0-1e-4, theta_0, theta_0+1e-4])
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',14,'TickLabelInterpreter', 'latex')
xlabel('Incident angle','FontSize',18, 'Interpreter','latex')
ylabel('Error','FontSize',18, 'Interpreter','latex')
