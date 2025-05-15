function [mesh,param] = GenerateMesh(param, domain)
%Define the triangulation on the domain specified 

if strcmp(domain,'double_rectangle')
    param.H=3; param.L=2*pi; %truncation parameter H and period L
    param.epsilon=[1, 1.5]; %vector of relative permittivity values inside the domain
    %Make sure that eps_+ is the first and eps_- is the last
    mesh.DirObst = false; %no Dirichlet obstacle

    %define triangulation - MUST BE PERIODIC
    [ mesh.p,mesh.t,mesh.I,mesh.B ] = MeshDouble(param);

    %define matrices for the element edges and divide them in
    %internal, upper, lower and left/right edges
    [mesh.LI,mesh.LU,mesh.LD,mesh.LR] = MeshEdges(mesh,param);

    %assign to every mesh element the corresponding epsilon value
    mesh.E = EpsValDouble(mesh);

elseif strcmp(domain,'triple_rectangle')
    param.H=5; param.L=2*pi; %truncation parameter H and period L
    param.L1=2; param.L2=-2; %adjust the heigth of the rectangles
    param.epsilon=[1, 2, 1]; %vector of relative permittivity values inside the domain
    %Make sure that eps_+ is the first and eps_- is the last
    mesh.DirObst = false; %no Dirichlet obstacle

    %define triangulation - MUST BE PERIODIC
    [ mesh.p,mesh.t,mesh.I,mesh.B ] = MeshTriple(param);

    %define matrices for the element edges and divide them in
    %internal, upper, lower and left/right edges
    [mesh.LI,mesh.LU,mesh.LD,mesh.LR] = MeshEdges(mesh,param);

    %assign to every mesh element the corresponding epsilon value
    mesh.E = EpsValTriple(mesh,param);

elseif strcmp(domain,'dir_square')
    param.H=5; param.L=2*pi; %truncation parameter H and period L
    param.height=1; %height of the square
    param.L1=3; param.L2=-3; %decide where epsilon varies
    param.epsilon=[1, 1, 1]; %vector of relative permittivity values inside the domain
    %Make sure that eps_+ is the first and eps_- is the last
    mesh.DirObst = true; %there is a Dirichlet obstacle

    %define triangulation - MUST BE PERIODIC
    [ mesh.p,mesh.t,mesh.I,mesh.B ] = MeshDirSquare(param);

    %define matrices for the element edges and divide them in
    %internal, upper, lower and left/right edges
    [mesh.LI,mesh.LU,mesh.LD,mesh.LR,mesh.LDir] = MeshEdgesDir(mesh,param);

    %assign to every mesh element the corresponding epsilon value
    mesh.E = EpsValTriple(mesh,param);

elseif strcmp(domain,'u_shape')
    param.H=2; param.L=2*pi; %truncation parameter H and period L
    param.epsilon=[1, (1.27+0.1*1i)^2]; %vector of relative permittivity values inside the domain
    %Make sure that eps_+ is the first and eps_- is the last
    mesh.DirObst = false; %no Dirichlet obstacle

    %define triangulation - MUST BE PERIODIC
    [ mesh.p,mesh.t,mesh.I,mesh.B ] = MeshUShape(param);

    %define matrices for the element edges and divide them in
    %internal, upper, lower and left/right edges
    [mesh.LI,mesh.LU,mesh.LD,mesh.LR] = MeshEdges(mesh,param);

    %assign to every mesh element the corresponding epsilon value
    mesh.E = EpsValUShape(mesh,param);

elseif strcmp(domain,'layers_rectangle')
    param.H=3; param.L=2*pi; %truncation parameter H and period L
    param.epsilon=[1, 1.49^2, 2.13^2, 2.02^2, 1.453^2]; %vector of relative permittivity values inside the domain
    %Make sure that eps_+ is the first and eps_- is the last
    mesh.DirObst = false; %no Dirichlet obstacle

    %define triangulation - MUST BE PERIODIC
    [ mesh.p,mesh.t,mesh.I,mesh.B ] = MeshLayersRect(param);

    %define matrices for the element edges and divide them in
    %internal, upper, lower and left/right edges
    [mesh.LI,mesh.LU,mesh.LD,mesh.LR] = MeshEdges(mesh,param);

    %assign to every mesh element the corresponding epsilon value
    mesh.E = EpsValLayersRect(mesh,param);

else
    error('Wrong domain!')
end

return