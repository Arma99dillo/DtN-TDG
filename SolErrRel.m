function [L2Err,H1Err,NormL2] = SolErrRel(mesh,param,u,param_raff,u_raff,phi,grad_phi)
%Given the refined solution, computes the L2 and H1 error

%get parameters
t=mesh.t; p=mesh.p; nd=param.nd; d=param.d;
K=param.K; epsilon=param.epsilon; E=mesh.E;
nd_raff=param_raff.nd; d_raff=param_raff.d; 

L2Err=0; H1Err=0; NormL2=0; NormH1=0;
%cycle on mesh elements
for T=1:size(t,1)

    k=K*sqrt(epsilon(E(T))); %wavenumber in the current element

    %solution coefficients corresponding to the element
    coeff=u((T-1)*nd+1:T*nd);
    coeff_raff=u_raff((T-1)*nd_raff+1:T*nd_raff);
    
    %get Duffy quadrature points
    Tr = [p(t(T,1),:); p(t(T,2),:); p(t(T,3),:)];
    Qp = quadtriangle(15,'domain',Tr);
    
    %compute refined and approximate solution in the quadrature points
    u_app=zeros(size(Qp.Points,1),1); u_app_r=zeros(size(Qp.Points,1),1); z=zeros(size(Qp.Points,1),1);
    g_app=zeros(2,size(Qp.Points,1)); g_app_r=zeros(2,size(Qp.Points,1)); w=zeros(size(Qp.Points,1),1);
    ml2=zeros(size(Qp.Points,1),1); mh1=zeros(size(Qp.Points,1),1);
    for l=1:size(Qp.Points,1)
        for j=1:nd
            dj=d(:,j);
            u_app(l)=u_app(l)+coeff(j).*phi(Qp.Points(l,1),Qp.Points(l,2),dj,k);
            g_app(:,l)=g_app(:,l)+coeff(j).*grad_phi(Qp.Points(l,1),Qp.Points(l,2),dj,k);
        end
        for j=1:nd_raff
            dj=d_raff(:,j);
            u_app_r(l)=u_app_r(l)+coeff_raff(j).*phi(Qp.Points(l,1),Qp.Points(l,2),dj,k);
            g_app_r(:,l)=g_app_r(:,l)+coeff_raff(j).*grad_phi(Qp.Points(l,1),Qp.Points(l,2),dj,k);
        end
        z(l) = abs(u_app_r(l)-u_app(l))^2; 
        ml2(l)=abs(u_app_r(l))^2;
        w(l) = abs(g_app_r(1,l)-g_app(1,l))^2 + abs(g_app_r(2,l)-g_app(2,l))^2;
        mh1(l) = abs(g_app_r(1,l))^2 + abs(g_app_r(2,l))^2;
    end
    
    %compute errors
    L2Err=L2Err+dot(Qp.Weights,z);
    NormL2=NormL2+dot(Qp.Weights,ml2);
    H1Err=H1Err+dot(Qp.Weights,z)+dot(Qp.Weights,w);
    NormH1=NormH1+dot(Qp.Weights,ml2)+dot(Qp.Weights,mh1);
end

L2Err=sqrt(L2Err)/sqrt(NormL2);
H1Err=sqrt(H1Err)/sqrt(NormH1);
NormL2=sqrt(NormL2);

end