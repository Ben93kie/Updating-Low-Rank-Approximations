function [ U,S,V] = Algorithm4( U,S,V,c)
%Given r-LRA USV^T, delete last c columns of LRA and give r-LRA back
%Possibly has c many zero singular values
m=length(U);
n=length(V);
D=zeros(n,c);
D(n-c+1:end,:)=eye(c,c);
r=min(size(S));
VD=V'*D;
[Q,RY]=qr(D-V*VD,0);
K=zeros(r,r+c);
K(1:r,1:r)=S-S*VD*D'*V;
K(1:r,r+1:r+c)=-S*VD*RY';
[Unew,S,Vnew]=svd(K);
Unew=Unew(:,1:r);
Vnew=Vnew(:,1:r);
S=S(1:r,1:r);
U=U*Unew;
V=[V,Q];
V=V*Vnew;
V=V(1:n-c,:);
end
