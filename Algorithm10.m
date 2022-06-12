function [ U,S,V ] = Algorithm10( U,S,V,DeltaA )
%Given r-DLRA USV^T compute r-DLRA USV^T+DeltaA
K=U*S;
DeltaAV=DeltaA*V;
K=K+DeltaAV;
[U,S]=qr(K,0);
S=S-U'*DeltaAV;
L=V*S'+DeltaA'*U;
[V,S]=qr(L,0);
S=S';
end

