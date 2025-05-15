function E = EpsValDouble(mesh)
%returns a vector indicating the region where each triangle is located

%get parameters
p=mesh.p; t=mesh.t;

E=zeros(size(t,1),1);
for j=1:size(t,1)
    bar = (p(t(j,1),:) + p(t(j,2),:) + p(t(j,3),:))/3; %triangle baricenter
    if bar(2) < 0 %lower region
        E(j)=2;
    else %upper region
        E(j)=1;
    end
end

return