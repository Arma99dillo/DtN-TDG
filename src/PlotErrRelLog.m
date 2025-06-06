function PlotErrRelLog(mesh,param,u,param_raff,u_raff,phi)
%plots absolute value of the approximate solution and numerical error

tic()

t1=mesh.t; p1=mesh.p; nd=param.nd; d=param.d; 
K=param.K; epsilon= param.epsilon; E=mesh.E;
nd_raff=param_raff.nd; d_raff=param_raff.d; 

figure()
%cycle on mesh elements
for T=1:size(t1,1)

    k=K*sqrt(epsilon(E(T))); %wavenumber in the current element
    
    %solution coefficients corresponding to the element
    coeff=u((T-1)*nd+1:T*nd);
    coeff_raff=u_raff((T-1)*nd_raff+1:T*nd_raff);
    
    %mesh on the element to plot the solution
    pv = [p1(t1(T,1),:); p1(t1(T,2),:); p1(t1(T,3),:)]; %vertices
    g = [2;3;pv(:,1);pv(:,2)]; % [3,4,x1,x2,x3,x4,y1,y2,y3,y4]'
    model = createpde();
    geometryFromEdges(model, decsg(g));
    m = generateMesh(model, 'Hmax', 0.01,'GeometricOrder','linear');
    p = m.Nodes'; t = m.Elements'; 

    %compute approximate and refined solution on the mesh vertices
    u_app=zeros(size(p,1),1); u_app_r=zeros(size(p,1),1); diff=zeros(size(p,1),1);
    for v=1:size(p,1)
        for j=1:nd
            dj=d(:,j);
            u_app(v)=u_app(v)+coeff(j).*phi(p(v,1),p(v,2),dj,k);
        end
        for j=1:nd_raff
            dj=d_raff(:,j);
            u_app_r(v)=u_app_r(v)+coeff_raff(j).*phi(p(v,1),p(v,2),dj,k);
        end
        diff(v)=u_app_r(v)-u_app(v); %difference
    end

    if T==1
        minColorLimit2 = min(min((log10(abs(diff)))));
        maxColorLimit2 = max(max((log10(abs(diff)))));
    else
        minColorLimit2 = min(min((log10(abs(diff)))),minColorLimit2);
        maxColorLimit2 = max(max((log10(abs(diff)))),maxColorLimit2);
    end

    subplot(1,2,1)
    trisurf(t,p(:,1),p(:,2),abs(u_app)); hold on

    subplot(1,2,2)
    trisurf(t,p(:,1),p(:,2),(log10(abs(diff)))); hold on
    
end

subplot(1,2,1)
grid off; colorbar; axis equal; shading flat; view(2); set(gca,'fontsize',12)

subplot(1,2,2)
shading flat;
cH = colorbar;
cH.Ticks = linspace(floor(minColorLimit2),ceil(maxColorLimit2),ceil(maxColorLimit2)-floor(minColorLimit2)+2);
cH.TickLabelInterpreter = 'tex';
for ii = 1:(ceil(maxColorLimit2)-floor(minColorLimit2)+2)
    cH.TickLabels{ii} = [sprintf('10^{%d}',floor(minColorLimit2)+ii-1)];
end
axis equal; grid off; hold on
view(2)
set(gca,'fontsize',12)

toc()

end