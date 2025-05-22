function A = MatrixDtNTDG(mesh,param)
%define the system matrix for the DTN-TDG method

tic()
disp('Started linear system assembly')

%get the parameters
t=mesh.t; p=mesh.p; B=mesh.B; K=param.K; LI=mesh.LI; alpha=param.alpha; 
beta=param.beta; delta=param.delta; epsilon=param.epsilon; eps_p=param.epsilon(1);
eps_m=param.epsilon(end); alp=param.alp; E=mesh.E; LU=mesh.LU; LD=mesh.LD; LR=mesh.LR;
M=param.M; H=param.H; S=param.L; nd=param.nd; d=param.d;

%sparse matrix allocation
m = size(t,1);
A=spalloc(m*nd,m*nd,10*m*nd);

R =[0,1;-1,0]; %rotation matrix (used for normal vectors)
t = [t t(:,1)]; %modify t for easily deal with the edges

%-----------------------------------
%Internal elements
%-----------------------------------

for T = 1:m 
    %cicle on all the elements to deal with matrix terms associated to
    %basis functions defined on the same element

    k=K*sqrt(epsilon(E(T))); %wavenumber inside the element

    %identify the internal edges
    %include the left and right boundary between the internal edges
    A_aux =zeros(nd,nd); %auxiliary matrix
    for u=1:3
        if IsInt(T,u,B,t,p,S) %check if it's an internal edge
            %cycle on PW directions: this way I include all basis functions
            for l = 1:nd
                for j= 1:nd
                    dl = d(:,l); dj = d(:,j); diff=dl-dj; %PW directions and difference
                    p1=p(t(T,u),:)'; p2=p(t(T,u+1),:)'; %endpoints
                    n = (R*(p2-p1))/norm(p2-p1); %outward normal
                    A_aux(j,l) = A_aux(j,l) + 1i*k.*(-(1/2).*(dot(dj,n)+dot(dl,n))-beta.*dot(dj,n).*dot(dl,n)-alpha).*phi_int(diff,p1,p2,k);
                end
            end
        end
    end    
    A((T-1)*nd+1:T*nd,(T-1)*nd+1:T*nd)=A((T-1)*nd+1:T*nd,(T-1)*nd+1:T*nd)+A_aux; %update matrix entries
end

%cycle on internal edges to deal with basis functions associated to two
%adjiacent elements
for L=1:size(LI,1)

    t1=LI(L,3); %first element T1
    t2=LI(L,4); %second element T2
    p1=p(LI(L,1),:)'; p2=p(LI(L,2),:)'; %endpoints
    n = (R*(p2-p1))/norm(p2-p1); %outwad normal to T1
    k1=K*sqrt(epsilon(E(t1))); k2=K*sqrt(epsilon(E(t2))); %wavenumbers inside T! and T2
    
    %first cycle: consider T1 for test and T2 for trial   
    A_aux =zeros(nd,nd); %auxiliary matrix
    for l=1:nd %cycle on PW directions
        for j=1:nd
        dl = d(:,l); dj = d(:,j); diff=k1.*dl-k2.*dj; %PW directions and difference
        A_aux(j,l) = A_aux(j,l) + (0.5*(1i*k2*dot(dj,n)+1i*k1*dot(dl,n))+alpha*1i*k2+beta*1i*k1*dot(dl,n)*dot(dj,n)).*phi_int(diff,p1,p2,1);           
        end
    end
    A((t2-1)*nd+1:t2*nd,(t1-1)*nd+1:t1*nd)=A((t2-1)*nd+1:t2*nd,(t1-1)*nd+1:t1*nd)+A_aux; %update matrix entries
    
    %second cycle: consider T2 for test and T1 for trial 
    %change the normal and the wavenumbers
    n=-n; k1=K*sqrt(epsilon(E(t2))); k2=K*sqrt(epsilon(E(t1)));  
    A_aux =zeros(nd,nd); %auxiliary matrix
    for l=1:nd %cycle on PW directions
        for j=1:nd
        dl = d(:,l); dj = d(:,j); diff=k1.*dl-k2.*dj; %PW directions and difference
        A_aux(j,l) = A_aux(j,l) + (0.5*(1i*k2*dot(dj,n)+1i*k1*dot(dl,n))+alpha*1i*k2+beta*1i*k1*dot(dl,n)*dot(dj,n)).*phi_int(diff,p1,p2,1);           
        end
    end
    A((t1-1)*nd+1:t1*nd,(t2-1)*nd+1:t2*nd)=A((t1-1)*nd+1:t1*nd,(t2-1)*nd+1:t2*nd)+A_aux; %update matrix entries 

end

%-----------------------------------
%Boundary elements
%-----------------------------------

%first I deal with computing terms which are useful for both DtN boundaries

%first compute all the beta_n^+ and beta_n^- terms for the Fourier 
%expansion, since they are used for all the basis functions
v=1; alpha_n=zeros(2*M+1,1); beta_n_p=zeros(2*M+1,1); beta_n_m=zeros(2*M+1,1);
for n=-M:M
    alpha_n(v)=alp+2*pi*n/S; 
    beta_n_p(v)=sqrt(K^2*eps_p-alpha_n(v)^2);
    beta_n_m(v)=sqrt(K^2*eps_m-alpha_n(v)^2);
    v=v+1;
end
%get all the Fourier coefficients for upper and lower boundary
kp=K*sqrt(eps_p); km=K*sqrt(eps_m);
[FCp,FCm]=FourierCoeff(p,LU,LD,S,H,kp,km,d,nd,M,alpha_n); 

%Upper boundary 
k=K*sqrt(eps_p); %wavenumber on the upper boundary
%cicle on all the upper elements to deal with global terms containg the DtN operator 
for L1=1:size(LU,1)
    for L2=1:size(LU,1)
        t1=LU(L1,3); %first element T1
        t2=LU(L2,3); %second element T2
        l1=[LU(L1,1:2)]; %endpoints of T1
        l2=[LU(L2,1:2)]; %endpoints of T2

        A_aux =zeros(nd,nd); %auxiliary matrix
        for l=1:nd %PW dircetions cycle
            for j=1:nd
                dl = d(:,l); dj = d(:,j); %PW directions
                [val1,val2,val3]=DtNInt(p,dl,dj,nd,k,l1,l2,L1,L2,j,l,FCp,M,alpha_n,beta_n_p,H,S); %get integral values
                A_aux(j,l) = A_aux(j,l) - delta*dl(2)*val2 + (-1+delta*dj(2))*val1 + (delta/(1i*k))*val3;
            end
        end
        A((t2-1)*nd+1:t2*nd,(t1-1)*nd+1:t1*nd)=A((t2-1)*nd+1:t2*nd,(t1-1)*nd+1:t1*nd)+A_aux; %update matrix
    end
end

%cycle on the upper elements for local terms
for L=1:size(LU,1)
    t1=LU(L,3); %adjacent element
    p1=p(LU(L,1),:)'; p2=p(LU(L,2),:)'; %endpoints
    n = [0;1]; %outward normal
    A_aux =zeros(nd,nd); %auxiliary matrix
    for l = 1:nd %PW directions cycle
        for j= 1:nd
            dl = d(:,l); dj = d(:,j); diff=dl-dj; %PW directions and difference
            A_aux(j,l) = A_aux(j,l) - 1i*k.*dot(dj,n).*(delta*dot(dl,n)+1).*phi_int(diff,p1,p2,k);
        end
    end
    A((t1-1)*nd+1:t1*nd,(t1-1)*nd+1:t1*nd)=A((t1-1)*nd+1:t1*nd,(t1-1)*nd+1:t1*nd)+A_aux; %update matrix
end


%Lower boundary
k=K*sqrt(eps_m); %wavenumber on the lower boundary
%cicle on all the lower elements to deal with global terms containg the DtN operator 
for L1=1:size(LD,1)
    for L2=1:size(LD,1)
        t1=LD(L1,3); %first element T1
        t2=LD(L2,3); %second element T2
        l1=[LD(L1,1:2)]; %endpoints of T1
        l2=[LD(L2,1:2)]; %endpoints of T2
        A_aux =zeros(nd,nd); %auxiliary matrix
        for l=1:nd %PW directions cycle
            for j=1:nd
                dl = d(:,l); dj = d(:,j); %PW directions
                [val1,val2,val3]=DtNInt_m(p,dl,dj,nd,k,l1,l2,L1,L2,j,l,FCm,M,alpha_n,beta_n_m,H,S); %get integral values
                A_aux(j,l) = A_aux(j,l) + delta*dl(2)*val2 -(1+delta*dj(2))*val1 + (delta/(1i*k))*val3;
            end
        end
        A((t2-1)*nd+1:t2*nd,(t1-1)*nd+1:t1*nd)=A((t2-1)*nd+1:t2*nd,(t1-1)*nd+1:t1*nd)+A_aux; %update matrix
    end
end

%cycle on the lower elements for local terms
for L=1:size(LD,1)
    t1=LD(L,3); %adjiacent element
    p1=p(LD(L,1),:)'; p2=p(LD(L,2),:)'; %endpoints
    n = [0;-1]; %outward normal
    A_aux =zeros(nd,nd); %auxiliary matrix
    for l = 1:nd %PW directions cycle
        for j= 1:nd
            dl = d(:,l); dj = d(:,j); diff=dl-dj; %PW directions and difference
            A_aux(j,l) = A_aux(j,l) - 1i*k.*dot(dj,n).*(delta*dot(dl,n)+1).*phi_int(diff,p1,p2,k);
        end
    end
    A((t1-1)*nd+1:t1*nd,(t1-1)*nd+1:t1*nd)=A((t1-1)*nd+1:t1*nd,(t1-1)*nd+1:t1*nd)+A_aux; %update matrix
end


%Left and right boundary
%deal with all the terms associated to quasi-periodicity: her we consider
%the left and right boundary as the same boundary
c=exp(1i*alp*S); %quasi-periodicity constant
for L=1:size(LR,1)
    t1=LR(L,3); %first element T1 (right)
    t2=LR(L,4); %second element T2 (left)
    p1=p(LR(L,1),:)'; p2=p(LR(L,2),:)'; %endpoints
    n = (R*(p2-p1))/norm(p2-p1); %outward normal to T1
    k=K*sqrt(epsilon(E(t1))); %wavenumber in the elements (it's the same)

    %first cycle: consider T1 for test and T2 for trial   
    A_aux =zeros(nd,nd); %auxiliary matrix
    for l=1:nd %PW directions cycle
        for j=1:nd
        dl = d(:,l); dj = d(:,j); %PW directions
        A_aux(j,l) = A_aux(j,l) - c'.*1i*k.*(-(1/2)*(dot(dj,n)+dot(dl,n))-beta*dot(dj,n)*dot(dl,n)-alpha).*phi_int_bound(dl,dj,p1,p2,k,1,S);            
        end
    end
    A((t2-1)*nd+1:t2*nd,(t1-1)*nd+1:t1*nd)=A((t2-1)*nd+1:t2*nd,(t1-1)*nd+1:t1*nd)+A_aux; %update matrix

    n = -n; %change the normal to be outward to T2
    %second cycle: consider T2 for test and T1 for trial  
    A_aux =zeros(nd,nd); %auxiliary matrix
    for l=1:nd %PW directions cycle
        for j=1:nd
        dl = d(:,l); dj = d(:,j); %PW directions
        A_aux(j,l) = A_aux(j,l) - c.*1i*k.*(-(1/2)*(dot(dj,n)+dot(dl,n))-beta*dot(dj,n)*dot(dl,n)-alpha).*phi_int_bound(dl,dj,p1,p2,k,0,S);          
        end
    end
    A((t1-1)*nd+1:t1*nd,(t2-1)*nd+1:t2*nd)=A((t1-1)*nd+1:t1*nd,(t2-1)*nd+1:t2*nd)+A_aux; %update matrix
end

%cycle on Dirichlet boundary if necessary
if mesh.DirObst
    LDir=mesh.LDir;
    for index=1:size(LDir,1)
        t1=LDir(index,3); %adjiacent element
        k=K*sqrt(epsilon(E(t1))); %wavenumber in the element
        p1=p(LDir(index,1),:)'; p2=p(LDir(index,2),:)'; %endpoints
        n = (R*(p2-p1))/norm(p2-p1); %outward normal
        A_aux =zeros(nd,nd); %auxiliary matrix
        for l=1:nd %PW directions cycle
            for j=1:nd
                dl = d(:,l); dj = d(:,j); diff=dl-dj; %PW directions and difference
                A_aux(j,l) = A_aux(j,l) + 1i*k.*(-alpha-dot(dl,n)).*phi_int(diff,p1,p2,k);
            end
        end
        A((t1-1)*nd+1:t1*nd,(t1-1)*nd+1:t1*nd)=A((t1-1)*nd+1:t1*nd,(t1-1)*nd+1:t1*nd)+A_aux; %update matrix
    end
end

return