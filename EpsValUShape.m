function E = EpsValUShape(mesh,param)
%returns a vector indicating the region where each triangle is located

%get parameters
p=mesh.p; t=mesh.t; H=param.H; L=param.L;

E=zeros(size(t,1),1);
for j=1:size(t,1)
    bar = (p(t(j,1),:) + p(t(j,2),:) + p(t(j,3),:))/3; %triangle baricenter
    if bar(2) < -0.5*H
        E(j)=2;
    elseif bar(2) > 0.5*H
        E(j)=1;
    else
        if bar(1)> L/4 && bar(1) < 3*L/4
            E(j)=2;
        else
            E(j)=1;
        end
    end
end

return