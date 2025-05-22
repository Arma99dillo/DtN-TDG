function [val1,val2,val3]=DtNInt(p,dl,dj,nd,k,l1,l2,L1,L2,j,l,coeff,M,alpha_n,beta_n,H,S)

toll=1e-10;

x1=p(l1(2),1); x2=p(l1(1),1); %endpoints phi_l
x1t=p(l2(2),1); x2t=p(l2(1),1); %endpoints phi_j

%get Fourier coefficients from the vector
coeffl=coeff((2*M+1)*nd*(L1-1)+(l-1)*(2*M+1)+1 : (2*M+1)*nd*(L1-1)+l*(2*M+1));
coeffj=coeff((2*M+1)*nd*(L2-1)+(j-1)*(2*M+1)+1 : (2*M+1)*nd*(L2-1)+j*(2*M+1));
 
%compute integrals of T(phi_l)*phi_j', phi_l*T(phi_j)' and T(phi_l)*T(phi_j)'
v=1; val1=0; val2=0; val3=0;
for n=-M:M
    if abs(alpha_n(v)-k*dj(1)) <= toll
        val1=val1+(1i*coeffl(v)*beta_n(v)).*exp(-1i*k*dj(2)*H).*(x2t-x1t);
    else
        val1=val1+(1i*coeffl(v)*beta_n(v)).*(exp(-1i*k*dj(2)*H)./(1i*(alpha_n(v)-k*dj(1)))).*(exp(1i*(alpha_n(v)-k*dj(1))*x2t)-exp(1i*(alpha_n(v)-k*dj(1))*x1t));
    end
    if abs(alpha_n(v)-k*dl(1)) <= toll
        val2=val2+(1i*coeffj(v)*beta_n(v))'.*exp(1i*k*dl(2)*H).*(x2-x1);
    else
        val2=val2+(1i*coeffj(v)*beta_n(v))'.*(exp(1i*k*dl(2)*H)./(-1i*(alpha_n(v)-k*dl(1)))).*(exp(-1i*(alpha_n(v)-k*dl(1))*x2)-exp(-1i*(alpha_n(v)-k*dl(1))*x1));
    end
    val3=val3+S.*abs(beta_n(v)).^2.*coeffl(v).*(coeffj(v))';
    v=v+1;
end

return