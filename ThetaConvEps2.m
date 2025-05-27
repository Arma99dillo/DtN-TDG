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
param.M=20; %number of Fourier modes

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
domain = 'triple_rectangle_eps2';

disp(['Theta-convergence test for epsilon=2, with k=', num2str(param.K)...
    ', p=', num2str(param.nd), ', h=', num2str(param.h)])

%-----------------------------------
%Cycle on theta
%-----------------------------------
%basis functions and derivatives
phi = @(x1,x2,d,k) exp(1i*k.*(x1.*d(1)+x2.*d(2)));
grad_phi = @(x1,x2,d,k) 1i*k.*d.*exp(1i*k.*(x1.*d(1)+x2.*d(2)));

thetavect=theta_0-1e-4:1e-5:theta_0+1e-4;
L2Error=zeros(size(thetavect));   %error vectors
for v=1:size(thetavect,2) %cycle on theta
    
    param.theta=thetavect(v); %incident angle
    param.alp=param.K*cos(param.theta); %quasi-periodicity parameter

    %generate mesh and compute exact solution
    [mesh,param,uex,uexdx,uexdy] = GenerateMeshSol(param,domain);
    
    disp([ 'Linear system assembly for incident angle theta=', num2str(param.theta) ])
    A = MatrixDtNTDG(mesh,param); %system matrix
    b = rhsDtNTDG(mesh,param); %system rhs
    u=A\b; %solve the system

    %Error computation
    [err2] = SolErr(mesh,param,u,phi,grad_phi,uex,uexdx,uexdy);
    L2Error(v) = err2;

end


%-----------------------------------
%Convergence plot
%-----------------------------------
figure()
semilogy(thetavect,L2Error,'*-','LineWidth',1.2); grid
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
