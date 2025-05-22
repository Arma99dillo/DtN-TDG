function PlotSolution(mesh,param,u,phi)
%plots real part and absolute value of the approximate solution

tic()
disp('Started solution plot')

%get parameters
t1=mesh.t; p1=mesh.p; nd=param.nd; d=param.d; K=param.K; 
E=mesh.E; epsilon=param.epsilon;

figure()
%cycle on all the mesh elements
for T=1:size(t1,1)

    k=K*sqrt(epsilon(E(T))); %wavenumber in the element
    
    %solution coefficients corresponding to the element
    coeff=u((T-1)*nd+1:T*nd);
    
    %mesh on the element to plot the solution
    pv = [p1(t1(T,1),:); p1(t1(T,2),:); p1(t1(T,3),:)]; %vertices
    g = [2;3;pv(:,1);pv(:,2)]; % [3,4,x1,x2,x3,x4,y1,y2,y3,y4]'
    model = createpde();
    geometryFromEdges(model, decsg(g));
    m = generateMesh(model, 'Hmax', 0.01,'GeometricOrder','linear');
    p = m.Nodes'; t = m.Elements'; 

    %compute approximate solution on the mesh vertices
    u_app=zeros(size(p,1),1);
    for v=1:size(p,1)
        for j=1:nd
            dj=d(:,j); %plane wave direction
            u_app(v)=u_app(v)+coeff(j).*phi(p(v,1),p(v,2),dj,k);
        end
    end
    
    %determine colorbar limits from data
    if T==1
        minColorLimit1 = min(min(real(u_app))); 
        maxColorLimit1 = max(max(abs(u_app)));
    else
        minColorLimit1 = min(min(real(u_app)),minColorLimit1);   
        maxColorLimit1 = max(max(abs(u_app)),maxColorLimit1);
    end
    
    %plot real part and absolute value
    subplot(1,2,1)
    trisurf(t,p(:,1),p(:,2),real(u_app)); 
    hold on

    subplot(1,2,2)
    trisurf(t,p(:,1),p(:,2),abs(u_app));
    hold on
    
end

subplot(1,2,1)
grid off; colorbar; axis equal; shading flat; view(2); set(gca,'fontsize',12)
title('Real part of numerical solution')

subplot(1,2,2)
grid off; colorbar; axis equal; shading flat; view(2); set(gca,'fontsize',12)
title('Absolute value of numerical solution')

disp('Ended solution plot')
toc()
end