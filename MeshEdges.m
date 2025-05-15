function [LI,LU,LD,LR] = MeshEdges(mesh,param)
%create internal and boundary edges matrices
%for every row, the matrix contains the vertices of the edge and the 
%triangle(s) adjacent to that edge

%get parameters
p=mesh.p; t=mesh.t; B=mesh.B; H=param.H; L=param.L;

%I have to distinguish between internal and boundary edges
control = 0; ni=0; nb=0; c=0; toll=1e-3;
for j = 1:size(t,1)
    
    Edges = [ t(j,2) t(j,3); t(j,3) t(j,1); t(j,1) t(j,2) ]; %triangle edges
    
    for l = 1:3  %cicle on the edges
        if ni+nb>=3
            li = [Edges(l,2) Edges(l,1)]; %flipped edge
            [~, index]=ismember(LI(:,1:2),li,'rows');
            if sum(index)== 1 %I check if I already included the edge
                %if it already appeared, it's an internal edge
                control = 1;
                %update LI matrix with the new triangle
                v=find(index);
                LI(v,4) = j;
            end
        end

        if control == 0 %if it is a new edge I add it to the matrix
            %check if it's a boundary edge
            if(ismember(Edges(l,1),B) && ismember(Edges(l,2),B))
                c=c+1;
            end
            if (abs(p(Edges(l,1),1)-p(Edges(l,2),1))<=toll || abs(p(Edges(l,1),2)-p(Edges(l,2),2))<=toll)
                c=c+1;
            end
            if c==2
                %boundary edge if both the nodes are on the boundary and
                %the edge is horizonthal or vertical
                nb=nb+1;
                LB(nb,1:2)= Edges(l,:);
                LB(nb,3)=j;
            else               
                ni = ni+1;
                LI(ni,1:2) = Edges(l,:);
                LI(ni,3) = j;
            end
        end
        control = 0;
        c=0;
    end
end

%divide boundary edges in left, right, upper and lower
nl=0; nr=0; nu=0; nd=0;
for j=1:size(LB,1)
    if abs(p(LB(j,1),1))<=toll && abs(p(LB(j,2),1))<=toll %left
        nl=nl+1;
        Left(nl)=j;
    elseif abs(p(LB(j,1),1)-L)<=toll && abs(p(LB(j,2),1)-L)<=toll %right
        nr=nr+1;
        Right(nr)=j;
    elseif abs(p(LB(j,1),2)-H)<=toll && abs(p(LB(j,2),2)-H)<=toll %upper
        nu=nu+1;
        LU(nu,:)=LB(j,:);
    else %lower
        nd=nd+1;
        LD(nd,:)=LB(j,:);
    end
end

%for the left and right edges I couple them with respect to the x2
%coordinate, creating the LR matrix
nlr=0;
for j=1:nr
    l1=[p(LB(Right(j),1),2),p(LB(Right(j),2),2)]; %x2 coordinates
    for k=1:nl
        l2=[p(LB(Left(k),2),2),p(LB(Left(k),1),2)]; %flipped x2 coordinates
        if(abs(l1-l2)<=toll)
            nlr=nlr+1;
            LR(nlr,1:3) = LB(Right(j),:);
            LR(nlr,4) = LB(Left(k),3);
        end
    end
end


return