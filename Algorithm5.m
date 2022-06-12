function [ UO,UI,Uinv ] = Algorithm5( UO,UI,Uinv,L )
%Given UO(mxr-matrix),UI,Uinv,L, extract a subspace dimension such that
%UO is mxr-1 matrix
k=0;
tol=10^(-1);
% flag=0;
r=length(UI);
[~,lengthUO]=size(UO);

if abs(L(end,end)) > tol && r<lengthUO
   
   UI=UI*L(1:end-1,1:end-1);
   Linv=L(1:r,1:r)'-L(end,1:end-1)'*L(1:end-1,end)'/L(end,end);
   Uinv=Linv*Uinv;
   UO=UO(:,1:end-1)+UO(:,end)*L(end,1:end-1)*Uinv;
   
        
else
    
    if r<lengthUO
        zerovector=zeros(r,1);
        UI=[UI zerovector;
            zerovector' 1]*L;
        Uinv=L'*[Uinv zerovector;
                zerovector' 1];
    else
        
        UI=UI*L;
        Uinv=L'*Uinv;
    end
    
        [~,k]=max(Uinv(end,:));
        if abs(UI(k,end)) < 10^-4
            
            k=0;
            
        end
    
    if k==0
    for j=1:r
        if abs(UI(j,end)) > 10^(-2) && abs(Uinv(end,j)) > 10^(-2)
           
            k=j;
            break
        end
    end
    end
    if k==0
        
    for j=1:r
        if abs(UI(j,end)) > 10^(-12) && abs(Uinv(end,j)) > 10^(-12)
            k=j;
            break
        end
    end
    end
        
    if k==0        
        flag=1
        return;
    end
        
    UI(:,end)=[];
    u=UI(k,:);
    UI(k,:)=[];
    Uinvkcol=Uinv(:,k);
    Uinvendend=Uinv(end,k);
    Uinvkcol(end)=[];
    Uinvlastrow=Uinv(end,:);
    Uinvlastrow(k)=[];
    Uinv(:,k)=[];
    Uinv(end,:)=[];
    Uinv=Uinv-Uinvkcol*Uinvlastrow/Uinvendend;
    UOkcol=UO(:,k);    
    UO(:,k)=[];    
    UO=UO+UOkcol*(u*Uinv);
end

    
    
end

