function [ UO,S,VO ] = Algorithm9_1( A,tol,maxrank )
%Given A compute r-LRA of A
%using stopping criterion
[~,n]=size(A);
UI=1;
VI=UI;
Uinv=UI;
Vinv=UI;
VO=UI;
UO=A(:,1);
normUO=norm(UO);
S=normUO;
if normUO<10^-14
    normUO=0;%dont expect zero column as first column of A
else
    UO=UO/normUO;
for j=2:n    
[ UO,VO,UI,VI,Uinv,Vinv,S ] = Algorithm7_1( UO,VO,UI,VI,Uinv,Vinv,S,A(:,j),tol,maxrank);
end
UO=UO*UI;
VO=VO*VI;
end

