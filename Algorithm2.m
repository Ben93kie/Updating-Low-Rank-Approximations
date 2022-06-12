function [ U,S,V ] = Algorithm2( U,S,V,C,D )
%Given r-LRA USV^T, compute r-LRA USV^T+CD^T
[m,~]=size(U);
r=min(size(S));
n=length(V);
[~,c]=size(C);
rplusc=r+c;
if r==m || r==n
    K=S+U'*C*D'*V;
    [Unew,S,Vnew]=svd(K);
    U=U*Unew(:,1:r);
    V=V*Vnew(:,1:r);
    S=S(1:r,1:r);
else 
    [P,RX] = qr((eye(m)-U*U')*C,0);
    [Q,RY] = qr((eye(n)-V*V')*D,0);
    if rplusc>min(m,n)
        K=zeros(min(m,n),min(m,n));
        K(1:r,1:r)=S;
        P=P(:,1:min(m,n)-r);
        RX=RX(1:min(m,n)-r,:);
        Q=Q(:,1:min(m,n)-r);
        RY=RY(1:min(m,n)-r,:);
    else
        K=zeros(rplusc,rplusc);
        K(1:r,1:r)=S;
    end
        K=K+[U'*C;RX]*[V'*D;RY]';

[Unew,S,Vnew]=svd(K);
S=S(1:r,1:r);
U=[U,P]*Unew(:,1:r);
V=[V,Q]*Vnew(:,1:r);
end
end
