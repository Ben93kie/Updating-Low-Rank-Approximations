function [ U,S,V] = Algorithm3( U,S,V,C)
%Given r-LRA USV^T compute [r-LRA USV^T,C]
m=length(U);
n=length(V);
r=min(size(S));
[~,c]=size(C);
[P,RC]=qr((eye(m)-U*U')*C,0);
K=zeros(r+c,r+c);
K(1:r,1:r)=S;
K(1:r,r+1:r+c)=U'*C;
K(r+1:r+c,r+1:r+c)=RC;
[Unew,S,Vnew]=svd(K);
Unew=Unew(:,1:r);
Vnew=Vnew(:,1:r);
S=S(1:r,1:r);
U=[U,P]*Unew;
Vz=zeros(n+c,r+c);
Vz(1:n,1:r)=V;
Vz(n+1:n+c,r+1:r+c)=eye(c,c);
V=Vz*Vnew;
end
