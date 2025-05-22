function val = phi_int(d,p1,p2,k)
%computes exactly the integral of exp(i*k*dot(x,d)) on the segment [p1,p2]

if norm(d)<=1e-10 %integrand is 1
    val = norm(p2-p1);
elseif abs(dot(d,p2-p1)) <= 1e-10 %dividend is zero
    val= norm(p2-p1).*exp(1i*k*dot(p1,d));
else
    val = (norm(p2-p1)/(1i*k*dot((p2-p1),d))).*exp(1i*k*dot(p1,d)).*expm1(1i*k*dot(p2-p1,d));
end

return