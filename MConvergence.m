% DtN-TDG solver for Helmholtz equation on periodic grating
% M-convergence error test

close all; addpath quadtriangle; addpath src

%-----------------------------------
%Parameters definition
%-----------------------------------
%Problem parameters
param.theta=-pi/3; %incident angle

%Discretization parameters
param.h=0.75; %mesh width
param.alpha=1/2; param.beta=1/2; param.delta=1/2; %TDG flux coefficients
param.nd=20; %number of directions
M_raff=100; %number of Fourier modes for the refined solution


%-----------------------------------
%Mesh definition
%-----------------------------------
domain = 'u_shape'; %select domain
[mesh,param] = GenerateMesh(param,domain); %generate mesh

%basis functions and derivatives
phi = @(x1,x2,d,k) exp(1i*k.*(x1.*d(1)+x2.*d(2)));
grad_phi = @(x1,x2,d,k) 1i*k.*d.*exp(1i*k.*(x1.*d(1)+x2.*d(2)));

%define the plane wave direction vectors
param.d=zeros(2,param.nd);
for l=1:param.nd
    param.d(:,l)=[cos((2*pi*l)/param.nd); sin((2*pi*l)/param.nd)];
end


%-----------------------------------
%Convergence test
%-----------------------------------

KK=[4,6,8,10]; %choose different wavenumbers
Mi=2; Mf=50; %min and max value of M
L2Error=zeros((Mf-Mi)/2+1,size(KK,2)); %error vectors

for j=1:size(KK,2)
    
    param.K=KK(j);
    param.alp=param.K*cos(param.theta); %quasi-periodicity parameter

    disp(['M-convergence test on the domain ', domain, ' with k=', num2str(param.K), ...
        ', p=', num2str(param.nd), ', h=', num2str(param.h),...
        ', with ', num2str(param.nd.*size(mesh.t,1)), ' basis functions' ])

    %Refined solution for relative error
    param_raff=param;
    param_raff.M=M_raff;

    disp(['Started computing refined solution (M=', num2str(M_raff), ')' ])

    %matrix, rhs and solution for the refined parameters
    A_raff = MatrixDtNTDG(mesh,param_raff);
    b_raff = rhsDtNTDG(mesh,param_raff);
    u_raff=A_raff\b_raff;


    %Cycle on number of Fourier modes
    v=1;
    for M=Mi:2:Mf

        param.M=M;

        A = MatrixDtNTDG(mesh,param); %system matrix
        b = rhsDtNTDG(mesh,param); %system rhs
        u=A\b; %solve the system

        %Error computation
        [err2,~] = SolErrRel(mesh,param,u,param_raff,u_raff,phi,grad_phi);
        L2Error(v,j) = err2;
        disp([ 'Computed error for k=', num2str(param.K) , ', M=', num2str(param.M), ' Fourier modes truncation' ])

        v=v+1;
    end

end


%-----------------------------------
%Convergence plot
%-----------------------------------
figure()
semilogy(Mi:2:Mf,L2Error(:,1),'o-',Mi:2:Mf,L2Error(:,2),'*-', ...
    Mi:2:Mf,L2Error(:,3),'^-',Mi:2:Mf,L2Error(:,4),'s-','LineWidth',1.2); grid
xlim([Mi Mf])
LL = legend('$k=4$','$k=6$','$k=8$','$k=10$','FontSize', 14);
set(LL, 'Interpreter', 'latex');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14,'TickLabelInterpreter', 'latex')
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',14,'TickLabelInterpreter', 'latex')
xlabel('Number of Fourier modes','FontSize',18, 'Interpreter','latex')
ylabel('$L^2$ Error','FontSize',18, 'Interpreter','latex')

