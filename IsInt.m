function int = IsInt(T,u,B,t,p,L)
%returns true if in internal element, false otherwise 
toll=1e-6;

if(ismember(t(T,u),B) && ismember(t(T,u+1),B)) && (abs(p(t(T,u),2)-p(t(T,u+1),2))<=toll)
    %boundary element
    int=false;
elseif(ismember(t(T,u),B) && ismember(t(T,u+1),B)) &&  (abs(p(t(T,u),1)-p(t(T,u+1),1))<=toll) && (abs(p(t(T,u),1))>toll && abs(p(t(T,u),1)-L)>toll)
    %boundary element
    int=false;
else
    int=true;
end

return