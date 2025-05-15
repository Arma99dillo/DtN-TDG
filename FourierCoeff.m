function [coeffp,coeffm]=FourierCoeff(p,LU,LD,S,H,kp,km,d,nd,M,alpha_n)
%calcolo coefficienti di Fourier delle funzioni sui lati DtN

toll=1e-6;

coeffp=zeros((2*M+1)*nd*size(LU,1),1); coeffm=zeros((2*M+1)*nd*size(LU,1),1); v=1;
for L=1:size(LU,1) %cicle on the elements
    x1=p(LU(L,2),1)'; x2=p(LU(L,1),1)'; %get the endpoints of the upper edge
    x1t=p(LD(L,1),1)'; x2t=p(LD(L,2),1)'; %and lower edge
    for l=1:nd %cicle on directions
        dl = d(:,l); u=1;
        for n=-M:M %cicle on Fourier modes
            if abs(kp*dl(1)-alpha_n(u)) <= toll
                coeffp(v)=(1/(S)).*exp(1i*kp*dl(2)*H).*(x2-x1);
            else
                coeffp(v)=(1/(S)).*(exp(1i*kp*dl(2)*H)./(1i*(kp*dl(1)-alpha_n(u)))).*(exp(1i*(kp*dl(1)-alpha_n(u))*x2)-exp(1i*(kp*dl(1)-alpha_n(u))*x1));
            end
            if abs(km*dl(1)-alpha_n(u)) <= toll
                coeffm(v)=(1/(S)).*exp(-1i*km*dl(2)*H).*abs(x2t-x1t);
            else
                coeffm(v)=(1/(S)).*(exp(-1i*km*dl(2)*H)./(1i*(km*dl(1)-alpha_n(u)))).*(exp(1i*(km*dl(1)-alpha_n(u))*x2t)-exp(1i*(km*dl(1)-alpha_n(u))*x1t));
            end
            v=v+1; u=u+1;
        end
    end
end

return