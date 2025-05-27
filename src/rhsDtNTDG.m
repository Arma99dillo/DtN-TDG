function b = rhsDtNTDG(mesh,param)
%computes the rhs for the DtN-TDG method

%get parameters
t=mesh.t; p=mesh.p; K=param.K; delta=param.delta; 
epsilon=param.epsilon; alp=param.alp; LU=mesh.LU; H=param.H; L=param.L; 
nd=param.nd; d=param.d; bet=param.K*sin(param.theta);

%rhs definition
m = size(t,1);
b = zeros(m*nd,1);
di=[cos(param.theta);sin(param.theta)]; %incident direction

toll=1e-6;  
k=K*sqrt(epsilon(1)); %wavenumber on the upper boundary

for F=1:size(LU,1) %cycle on the upper boundary
    t1=LU(F,3); %adjiacent element
    p1=p(LU(F,1),:)'; p2=p(LU(F,2),:)'; %endpoints
    x1=p2(1); x2=p1(1); %segment endpoints

    for j=1:nd %PW directions cycle
        dj = d(:,j); diff=di-dj; %PW directions

        %0-th Fourier coefficients of phi_j
        if abs(k*dj(1)-alp) <= toll
            phi_j_0=(1/L).*exp(1i*k*dj(2)*H).*(x2-x1);
        else
            phi_j_0=(1/L).*(exp(1i*k*dj(2)*H)./(1i*(k*dj(1)-alp))).*(exp(1i*(k*dj(1)-alp)*x2)-exp(1i*(k*dj(1)-alp)*x1));
        end
        %update rhs
        b((t1-1)*nd+j)=b((t1-1)*nd+j) + 2*1i*bet*(1-delta*dj(2))*phi_int(diff,p1,p2,k) - (1i*delta*2*L*abs(bet)^2*phi_j_0'/k).*exp(1i*bet*H);
    end
end

toc()

return