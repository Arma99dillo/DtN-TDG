function val = phi_int_bound(dl,dj,p1,p2,k,l,L)
%computes the integral of exp(i*k*dot(x,d)) on the right vertical segment [p1,p2]

y1=p1(2); y2=p2(2);

if l==1 %right test, left trial
    cost=exp(1i*k*L*dl(1));
else %left test, right trial
    cost=exp(-1i*k*L*dj(1));
end

if abs(dl(2)-dj(2))<=1e-6
    val = cost*(y2-y1);
else
    val = (cost/(1i*k*(dl(2)-dj(2)))).*(exp(1i*k*(dl(2)-dj(2))*y2)-exp(1i*k*(dl(2)-dj(2))*y1));
end

return