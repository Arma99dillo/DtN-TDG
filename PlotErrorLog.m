function PlotErrorLog(mesh,param,u,phi,uex)
%plots real part and absolute value of the approximate solution
%and real part and absolute value of numerical error

t1=mesh.t; p1=mesh.p; nd=param.nd; d=param.d; 
K=param.K; epsilon= param.epsilon; E=mesh.E;

figure()
%cycle on mesh elements
for T=1:size(t1,1)

    k=K*sqrt(epsilon(E(T))); %wavenumber in the current element
    
    %solution coefficients corresponding to the element
    coeff=u((T-1)*nd+1:T*nd);
    
    %mesh on the element to plot the solution
    pv = [p1(t1(T,1),:); p1(t1(T,2),:); p1(t1(T,3),:)]; %vertices
    g = [2;3;pv(:,1);pv(:,2)]; % [3,4,x1,x2,x3,x4,y1,y2,y3,y4]'
    model = createpde();
    geometryFromEdges(model, decsg(g));
    m = generateMesh(model, 'Hmax', 0.01,'GeometricOrder','linear');
    p = m.Nodes'; t = m.Elements'; 
    
    %get exact solution expression in the current element
    u_ex=uex{E(T)};

    %compute approximate and exact solution on the mesh vertices
    u_app=zeros(size(p,1),1); u_exact=zeros(size(p,1),1); diff=zeros(size(p,1),1);
    for v=1:size(p,1)
        for j=1:nd
            dj=d(:,j);
            u_app(v)=u_app(v)+coeff(j).*phi(p(v,1),p(v,2),dj,k);
        end
        u_exact(v)=u_ex(p(v,:)');
        diff(v)=u_exact(v)-u_app(v); %difference
    end

    if T==1
        minColorLimit1 = min(min(real(u_app)));   % determine colorbar limits from data
        maxColorLimit1 = max(max(abs(u_app)));
        minColorLimit2 = min(min((log10(abs(diff)))));
        maxColorLimit2 = max(max((log10(abs(diff)))));
    else
        minColorLimit1 = min(min(real(u_app)),minColorLimit1);   % determine colorbar limits from data
        maxColorLimit1 = max(max(abs(u_app)),maxColorLimit1);
        minColorLimit2 = min(min((log10(abs(diff)))),minColorLimit2);
        maxColorLimit2 = max(max((log10(abs(diff)))),maxColorLimit2);
    end

    subplot(1,3,1)
    trisurf(t,p(:,1),p(:,2),real(u_app)); shading flat; 
    clim([minColorLimit1,maxColorLimit1]); axis equal; hold on
    view(2) 
    set(gca,'fontsize',12)

    subplot(1,3,2)
    trisurf(t,p(:,1),p(:,2),abs(u_app)); shading flat; 
    clim([minColorLimit1,maxColorLimit1]); colorbar; axis equal; hold on
    view(2)
    set(gca,'fontsize',12)

    subplot(1,3,3)
    trisurf(t,p(:,1),p(:,2),log10(abs(diff))); shading flat; 
    cH = colorbar;
    cH.Ticks = linspace(floor(minColorLimit2),ceil(maxColorLimit2),ceil(maxColorLimit2)-floor(minColorLimit2)+2);
    cH.TickLabelInterpreter = 'tex';
    for ii = 1:(ceil(maxColorLimit2)-floor(minColorLimit2)+2)
        cH.TickLabels{ii} = [sprintf('10^{%d}',floor(minColorLimit2)+ii-1)];
    end
    axis equal; hold on
    view(2)
    set(gca,'fontsize',12)
    
end

end