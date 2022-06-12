function [ UO,VO,UI,VI,Uinv,Vinv,S] = Algorithm7_2( UO,VO,UI,VI,Uinv,Vinv,S,x,tol,maxrank)
%Given extended r-LRA UOUISVI^TVO^T compute extended r-LRA
%[UOUISVI^TVO^T,x] with truncation criterion (truncate largest triplet)

n=length(VO);
m=UI'*(UO'*x);
r=length(m);
p=x-UO*(UI*m);
 
normp=norm(p);
if ~exist('maxrank','var')
     % maxrank does not exist, so default is sqrt(number of rows)
      maxrank = sqrt(rows);
 end
if normp > tol 
    zerovector=zeros(r,1);
    K=[S m;
        zeros(1,r) normp];
    [L,S,R]=svd(K);
    if r+1<=maxrank
    UO=[UO,p/normp];
    VO=[VO zeros(n,1);
        zeros(1,r) 1];
    UI=[UI zerovector;
        zerovector' 1]*L;
    VI=[VI zerovector;
        zerovector' 1]*R;
    Uinv=L'*[Uinv zerovector;
            zerovector' 1];
    Vinv=R'*[Vinv zerovector;
            zerovector' 1];
    else 

        [ UO,UI,Uinv ] = Algorithm5( [UO,p/normp],UI,Uinv,L );
        [ VO,VI,Vinv ] = Algorithm5( [VO zeros(n,1);
                                        zeros(1,r) 1],VI,Vinv,R );

        S=S(1:r,1:r);
    end
        

else 
    K=[S m];
    [L,S,R]=svd(K);
    VI=VI*R(1:r,1:r);
    Rsubinv=R(1:r,1:r)'-R(end,1:end-1)'*R(1:end-1,end)'/R(end,end);
    Vinv=Rsubinv*Vinv;
    VO=[VO;R(end,1:end-1)*Vinv];
    UI=UI*L;
    Uinv=L'*Uinv;
    S=S(:,1:r);

end

        



end
