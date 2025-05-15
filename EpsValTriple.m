function E = EpsValTriple(mesh,param)
%returns a vector indicating the region where each triangle is located

%get parameters
p=mesh.p; t=mesh.t; L1=param.L1; L2=param.L2;

E=zeros(size(t,1),1);
for k=1:size(t,1)
    bar = (p(t(k,1),:) + p(t(k,2),:) + p(t(k,3),:))/3; %triangle baricenter
    if bar(2) < L2 %lower region
        E(k)=3;
    elseif bar(2) > L1 %upper region
        E(k)=1;
    else %middle region
        E(k)=2;
    end
end

return