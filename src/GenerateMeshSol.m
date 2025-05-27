function [mesh,param,uex, uexdx, uexdy] = GenerateMeshSol(param, domain)
%Define the triangulation on the domain specified 

if strcmp(domain,'double_rectangle_abs')
    param.H=3; param.L=2*pi; %truncation parameter H and period L
    param.epsilon=[1, (1.25+0.1*1i)^2]; %vector of relative permittivity values inside the domain
    %Make sure that eps_+ is the first and eps_- is the last
    mesh.DirObst = false; %no Dirichlet obstacle

    %define triangulation - MUST BE PERIODIC
    [ mesh.p,mesh.t,mesh.I,mesh.B ] = MeshDouble(param);

    %define matrices for the element edges and divide them in
    %internal, upper, lower and left/right edges
    [mesh.LI,mesh.LU,mesh.LD,mesh.LR] = MeshEdges(mesh,param);

    %assign to every mesh element the corresponding epsilon value
    mesh.E = EpsValDouble(mesh);

    %compute exact solution and gradient
    param.d1=cos(param.theta); param.d2=sin(param.theta);
    param.gamma=sqrt(param.epsilon(2)-param.d1^2);
    param.R=(param.d2+param.gamma)/(param.d2-param.gamma);
    param.T=(2*param.d2)/(param.d2-param.gamma);
    
    uex = {@(x) exp(1i*param.K.*(x(1).*param.d1+x(2).*param.d2)) + ...
        param.R.*exp(1i*param.K.*(x(1).*param.d1-x(2).*param.d2)), ...
        @(x) param.T.*exp(1i*param.K.*(x(1).*param.d1-x(2).*param.gamma))};

    uexdx = {@(x) 1i*param.K*param.d1.*exp(1i*param.K.*(x(1).*param.d1+x(2).*param.d2)) + ...
        1i*param.K*param.d1*param.R.*exp(1i*param.K.*(x(1).*param.d1-x(2).*param.d2)), ...
        @(x) 1i*param.K*param.d1*param.T.*exp(1i*param.K.*(x(1).*param.d1-x(2).*param.gamma))};
    
    uexdy = {@(x) 1i*param.K*param.d2.*exp(1i*param.K.*(x(1).*param.d1+x(2).*param.d2)) - ...
        1i*param.K*param.d2*param.R.*exp(1i*param.K.*(x(1).*param.d1-x(2).*param.d2)), ...
        @(x) -1i*param.K*param.gamma*param.T.*exp(1i*param.K.*(x(1).*param.d1-x(2).*param.gamma))};

elseif strcmp(domain,'double_rectangle_no_abs')
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

    %compute exact solution and gradient
    param.d1=cos(param.theta); param.d2=sin(param.theta);
    param.gamma=sqrt(param.epsilon(2)-param.d1^2);
    param.R=(param.d2+param.gamma)/(param.d2-param.gamma);
    param.T=(2*param.d2)/(param.d2-param.gamma);
    
    uex = {@(x) exp(1i*param.K.*(x(1).*param.d1+x(2).*param.d2)) + ...
        param.R.*exp(1i*param.K.*(x(1).*param.d1-x(2).*param.d2)), ...
        @(x) param.T.*exp(1i*param.K.*(x(1).*param.d1-x(2).*param.gamma))};

    uexdx = {@(x) 1i*param.K*param.d1.*exp(1i*param.K.*(x(1).*param.d1+x(2).*param.d2)) + ...
        1i*param.K*param.d1*param.R.*exp(1i*param.K.*(x(1).*param.d1-x(2).*param.d2)), ...
        @(x) 1i*param.K*param.d1*param.T.*exp(1i*param.K.*(x(1).*param.d1-x(2).*param.gamma))};
    
    uexdy = {@(x) 1i*param.K*param.d2.*exp(1i*param.K.*(x(1).*param.d1+x(2).*param.d2)) - ...
        1i*param.K*param.d2*param.R.*exp(1i*param.K.*(x(1).*param.d1-x(2).*param.d2)), ...
        @(x) -1i*param.K*param.gamma*param.T.*exp(1i*param.K.*(x(1).*param.d1-x(2).*param.gamma))};

elseif strcmp(domain,'triple_rectangle_eps2')
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

    %compute exact solution and gradient
    param.d1=cos(param.theta); param.d2=sin(param.theta);
    param.gamma=sqrt(param.epsilon(2)-param.d1^2);
    P = [exp(-1i*param.K*param.L1*param.d2), -exp(-1i*param.K*param.L1*param.gamma), -exp(1i*param.K*param.L1*param.gamma), 0;
        -param.d2*exp(-1i*param.K*param.L1*param.d2), param.gamma*exp(-1i*param.K*param.L1*param.gamma), -param.gamma*exp(1i*param.K*param.L1*param.gamma), 0;
        0, exp(1i*param.K*param.L1*param.gamma), exp(-1i*param.K*param.L1*param.gamma), -exp(-1i*param.K*param.L1*param.d2);
        0, -param.gamma*exp(1i*param.K*param.L1*param.gamma), param.gamma*exp(-1i*param.K*param.L1*param.gamma), -param.d2*exp(-1i*param.K*param.L1*param.d2)];
    l = [-exp(1i*param.K*param.L1*param.d2); -param.d2*exp(1i*param.K*param.L1*param.d2); 0; 0];
    
    S = P\l; param.R=S(1); param.T1=S(2); param.T2=S(3); param.T3=S(4);

    uex = {@(x) exp(1i*param.K.*(x(1).*param.d1+x(2).*param.d2)) + ...
        param.R.*exp(1i*param.K.*(x(1).*param.d1-x(2).*param.d2)), ...
        @(x) param.T1.*exp(1i*param.K.*(x(1).*param.d1-x(2).*param.gamma)) + ...
        param.T2.*exp(1i*param.K.*(x(1).*param.d1+x(2).*param.gamma)), ...
        @(x) param.T3.*exp(1i*param.K.*(x(1).*param.d1+x(2).*param.d2))};
    
    uexdx = {@(x) 1i*param.K*param.d1.*exp(1i*param.K.*(x(1).*param.d1+x(2).*param.d2)) + ...
        1i*param.K*param.d1.*param.R.*exp(1i*param.K.*(x(1).*param.d1-x(2).*param.d2)), ...
        @(x) 1i*param.K*param.d1.*param.T1.*exp(1i*param.K.*(x(1).*param.d1-x(2).*param.gamma)) + ...
        1i*param.K*param.d1.*param.T2.*exp(1i*param.K.*(x(1).*param.d1+x(2).*param.gamma)), ...
        @(x) 1i*param.K*param.d1.*param.T3.*exp(1i*param.K.*(x(1).*param.d1+x(2).*param.d2))};

    uexdy = {@(x) 1i*param.K*param.d2.*exp(1i*param.K.*(x(1).*param.d1+x(2).*param.d2)) - ...
        1i*param.K*param.d2.*param.R.*exp(1i*param.K.*(x(1).*param.d1-x(2).*param.d2)), ...
        @(x) -1i*param.K*param.gamma.*param.T1.*exp(1i*param.K.*(x(1).*param.d1-x(2).*param.gamma)) + ...
        1i*param.K*param.gamma.*param.T2.*exp(1i*param.K.*(x(1).*param.d1+x(2).*param.gamma)), ...
        @(x) 1i*param.K*param.d2.*param.T3.*exp(1i*param.K.*(x(1).*param.d1+x(2).*param.d2))};

elseif strcmp(domain,'triple_rectangle_eps10')
    param.H=5; param.L=2*pi; %truncation parameter H and period L
    param.L1=2; param.L2=-2; %adjust the heigth of the rectangles
    param.epsilon=[1, 10, 1]; %vector of relative permittivity values inside the domain
    %Make sure that eps_+ is the first and eps_- is the last
    mesh.DirObst = false; %no Dirichlet obstacle

    %define triangulation - MUST BE PERIODIC
    [ mesh.p,mesh.t,mesh.I,mesh.B ] = MeshTriple(param);

    %define matrices for the element edges and divide them in
    %internal, upper, lower and left/right edges
    [mesh.LI,mesh.LU,mesh.LD,mesh.LR] = MeshEdges(mesh,param);

    %assign to every mesh element the corresponding epsilon value
    mesh.E = EpsValTriple(mesh,param);

    %compute exact solution and gradient
    param.d1=cos(param.theta); param.d2=sin(param.theta);
    param.gamma=sqrt(param.epsilon(2)-param.d1^2);
    P = [exp(-1i*param.K*param.L1*param.d2), -exp(-1i*param.K*param.L1*param.gamma), -exp(1i*param.K*param.L1*param.gamma), 0;
        -param.d2*exp(-1i*param.K*param.L1*param.d2), param.gamma*exp(-1i*param.K*param.L1*param.gamma), -param.gamma*exp(1i*param.K*param.L1*param.gamma), 0;
        0, exp(1i*param.K*param.L1*param.gamma), exp(-1i*param.K*param.L1*param.gamma), -exp(-1i*param.K*param.L1*param.d2);
        0, -param.gamma*exp(1i*param.K*param.L1*param.gamma), param.gamma*exp(-1i*param.K*param.L1*param.gamma), -param.d2*exp(-1i*param.K*param.L1*param.d2)];
    l = [-exp(1i*param.K*param.L1*param.d2); -param.d2*exp(1i*param.K*param.L1*param.d2); 0; 0];
    
    S = P\l; param.R=S(1); param.T1=S(2); param.T2=S(3); param.T3=S(4);

    uex = {@(x) exp(1i*param.K.*(x(1).*param.d1+x(2).*param.d2)) + ...
        param.R.*exp(1i*param.K.*(x(1).*param.d1-x(2).*param.d2)), ...
        @(x) param.T1.*exp(1i*param.K.*(x(1).*param.d1-x(2).*param.gamma)) + ...
        param.T2.*exp(1i*param.K.*(x(1).*param.d1+x(2).*param.gamma)), ...
        @(x) param.T3.*exp(1i*param.K.*(x(1).*param.d1+x(2).*param.d2))};
    
    uexdx = {@(x) 1i*param.K*param.d1.*exp(1i*param.K.*(x(1).*param.d1+x(2).*param.d2)) + ...
        1i*param.K*param.d1.*param.R.*exp(1i*param.K.*(x(1).*param.d1-x(2).*param.d2)), ...
        @(x) 1i*param.K*param.d1.*param.T1.*exp(1i*param.K.*(x(1).*param.d1-x(2).*param.gamma)) + ...
        1i*param.K*param.d1.*param.T2.*exp(1i*param.K.*(x(1).*param.d1+x(2).*param.gamma)), ...
        @(x) 1i*param.K*param.d1.*param.T3.*exp(1i*param.K.*(x(1).*param.d1+x(2).*param.d2))};

    uexdy = {@(x) 1i*param.K*param.d2.*exp(1i*param.K.*(x(1).*param.d1+x(2).*param.d2)) - ...
        1i*param.K*param.d2.*param.R.*exp(1i*param.K.*(x(1).*param.d1-x(2).*param.d2)), ...
        @(x) -1i*param.K*param.gamma.*param.T1.*exp(1i*param.K.*(x(1).*param.d1-x(2).*param.gamma)) + ...
        1i*param.K*param.gamma.*param.T2.*exp(1i*param.K.*(x(1).*param.d1+x(2).*param.gamma)), ...
        @(x) 1i*param.K*param.d2.*param.T3.*exp(1i*param.K.*(x(1).*param.d1+x(2).*param.d2))};

else
    error('Wrong domain!')
end

return