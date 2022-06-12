function [ UO,S,VO ] = Algorithm8_2( UO,S,VO,DeltaA,tol,maxrank )
%Given r-LRA USV^T compute r-LRA of USV^T+DeltaA
%using truncation criterion

r=length(S);
[~,n]=size(DeltaA);
UI=eye(r,r);
VI=UI;
Uinv=UI;
Vinv=UI;
for j=1:n
[ UO,VO,UI,VI,Uinv,Vinv,S ] = Algorithm6_2( UO,VO,UI,VI,Uinv,Vinv,S,DeltaA(:,j),j,tol,maxrank );

end
UO=UO*UI;
VO=VO*VI;

