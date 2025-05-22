function [p_new,t_new,I,B] = MeshLayersRect(param)

%Get the parameters
H=param.H; h=param.h; L=param.L; 

% Create a PDE model
model = createpde();

% Define the geometry 
% rectangle with vertices v = [v1; v2; v3; v4] 
v1 = [0,-H; 0,-2*H/3; L,-2*H/3; L,-H];  
g1 = [3;4;v1(:,1);v1(:,2)]; % [3,4,x1,x2,x3,x4,y1,y2,y3,y4]'
v2 = [0,-2*H/3; 0,-H/3; L,-H/3; L,-2*H/3]; 
g2 = [3;4;v2(:,1);v2(:,2)]; % [3,4,x1,x2,x3,x4,y1,y2,y3,y4]'
v3 = [0,-H/3; 0,0; L,0; L,-H/3]; 
g3 = [3;4;v3(:,1);v3(:,2)]; % [3,4,x1,x2,x3,x4,y1,y2,y3,y4]'
v4 = [0,0; 0,H; L,H; L,0]; 
g4 = [3;4;v4(:,1);v4(:,2)]; % [3,4,x1,x2,x3,x4,y1,y2,y3,y4]'
v5 = [L/4,0; L/4,2*H/3; 3*L/4,2*H/3; 3*L/4,0]; 
g5 = [3;4;v5(:,1);v5(:,2)]; % [3,4,x1,x2,x3,x4,y1,y2,y3,y4]'

geom_matrix = [g1 g2 g3 g4 g5];
set_formula = 'g1+g2+g3+g4+g5';     

% Add the geometry to the PDE model
geometryFromEdges(model, decsg(geom_matrix, set_formula, char('g1', 'g2', 'g3','g4','g5')'));

% Generate the mesh with Hmax set to h
mesh = generateMesh(model, 'Hmax', h,'GeometricOrder','linear');

% % Plot the mesh
pdeplot(model);
axis equal;

p = mesh.Nodes'; % (n_nodes x 2) coordinates of the vertices
t = mesh.Elements'; % (n_t x 2) indices of the vertices

% order triangles and nodes putting the internal nodes 
% before the boundary ones
fd = @(p) min(min(min(H+p(:,2),H-p(:,2)),p(:,1)),L-p(:,1)); %function used to get boundary elements
%fd has to be zero on the external boundary of Omega

n_nodes = size(p,1);
n_triangles = size (t,1);


% anti-clockwise orienting of vertices
 for i_t=1:n_triangles
            v=p(t(i_t,:),:);
            if det([v(:,1) v(:,2) ones(3,1)]) < 0
                beep
                t(i_t,[1 2 3])= t(i_t,[2 1 3]);
            end
 end
 
% nodes reordering
old2new = zeros(n_nodes,1);
new2old = zeros(n_nodes,1);
tol = 10^-3;
i_int = 1;
i_bound = n_nodes;
for i = 1:n_nodes
    if abs(fd(p(i,:))) < tol
        old2new(i)=i_bound;
        new2old(i_bound)=i;
        i_bound = i_bound -1;
    else
        old2new(i)=i_int;
        new2old(i_int)=i;
        i_int = i_int +1;
    end
end

% ordered mesh
p_new = p(new2old,:); 
t_new = old2new(t);
B=(i_bound+1):n_nodes; %boundary indices
I=1:(i_int-1); %internal indices

end