function [ UO,VO,UI,VI,Uinv,Vinv,S ] = Algorithm6_1( UO,VO,UI,VI,Uinv,Vinv,S,x,j,tol,maxrank )
%Given extended r-LRA UOUISVI^TVO^T compute extended r-LRA
%UOUISVI^TVO^T+xe_j with stopping criterion (set p=0 after r>maxrank)
rows=length(UO);
n=length(VO);
ej=zeros(n,1);
ej(j)=1;
m=UI'*(UO'*x);
r=length(m);
p=x-UO*(UI*m);
normp=norm(p);
n=VI'*(VO'*ej);
q=ej-VO*(VI*n);
normq=norm(q);
if ~exist('maxrank','var')
     % maxrank does not exist, so default is sqrt(number of rows)
      maxrank = sqrt(rows);
 end
if r+1>maxrank && normp>tol && normq>tol
    normp=0;
end
if normp > tol && normq > tol
    zerovector=zeros(r,1);
    K=[S+m*n' m*normq;
        n'*normp normp*normq];
    [L,S,R]=svd(K);
    UO=[UO,p/normp];
    VO=[VO,q/normq];
    UI=[UI zerovector;
        zerovector' 1]*L;
    VI=[VI zerovector;
        zerovector' 1]*R;
    Uinv=L'*[Uinv zerovector;
            zerovector' 1];
    Vinv=R'*[Vinv zerovector;
            zerovector' 1];
elseif normp > tol && normq < tol
    K=[S+m*n';n'*normp];
    [L,S,R]=svd(K);
    UO=[UO,p/normp];
    [UO,UI,Uinv]=Algorithm5( UO,UI,Uinv,L );
    VI=VI*R;
    Vinv=R'*Vinv;
    S=S(1:r,:);
elseif normp < tol && normq > tol
    K=[S+m*n' normq*m];
    [L,S,R]=svd(K);
    VO=[VO,q/normq];
    [VO,VI,Vinv]=Algorithm5( VO,VI,Vinv,R );
    UI=UI*L;
    Uinv=L'*Uinv;
    S=S(:,1:r);
else 
    K=S+m*n';
    [L,S,R]=svd(K);
    [UO,UI,Uinv]=Algorithm5( UO,UI,Uinv,L );
    [VO,VI,Vinv]=Algorithm5( VO,VI,Vinv,R );
    S=S(1:end-1,1:end-1);
end

        

end
