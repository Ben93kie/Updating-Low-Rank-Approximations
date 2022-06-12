%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 8a)
clear all
dimension=100;
ranki=10;
tol=10^-12;
d=zeros(dimension,1);
for l=1:dimension
    d(l)=2^(-l/10);
end
D=diag(d);
T1=rand(dimension,dimension);
T2=rand(dimension,dimension);
T1=T1-diag(diag(T1));
T2=T2-diag(diag(T2));
T1=T1-T1';
T1=T1*1/5;
T2=T2-T2';
T2=T2*1/5;
Q1=eye(dimension,dimension);
Q2=Q1;
intervallength=1;
h=1/100;
A=D;
normend1=zeros(96,1);
normend2=zeros(96,1);
normend3=zeros(96,1);
normend4=zeros(96,1);
normend5=zeros(96,1);
normend6=zeros(96,1);


[U,S,V]=svd(A);
U=U(:,1:ranki);
S=S(1:ranki,1:ranki);
V=V(:,1:ranki);

U3=U;U4=U;
S3=S;S4=S;
V3=V;V4=V;
norm1=zeros(21,1);
norm1(1)=norm(U*S*V'-A,'fro');
norm2=zeros(21,1);

norm2(1,1)=norm(U*S*V'-A,'fro');


norm3=norm2;
norm4=norm2;
norm5=norm2;
norm6=norm2;
normg=norm1;
Aold=A;
for j=1:100  
    Q1=((eye(dimension,dimension)-h*T1))\Q1;
    Q2=((eye(dimension,dimension)-h*T2))\Q2;
    
    A=Q1*(exp(j/100)*D)*Q2';
    if mod(j,5)==0
        testdelta=zeros(dimension,dimension);
        testdelta(1,:)=ones(1,dimension);
        DeltaA=A-Aold;       
        [ U,S,V ] = Algorithm10( U,S,V,DeltaA );           
        [ U1,S1,V1 ] = Algorithm9_2( A,tol,ranki );      
        [ U2,S2,V2 ] = Algorithm9_1( A,tol,ranki );       
        [ U3,S3,V3 ] = Algorithm8_2( U3,S3,V3,DeltaA,tol,ranki );       
        [ U4,S4,V4 ] = Algorithm8_1( U4,S4,V4,DeltaA,tol,ranki ); 
        
        norm2(j/5+1,1)=norm(U1*S1*V1'-A,'fro');
        norm3(j/5+1,1)=norm(U2*S2*V2'-A,'fro');
        norm4(j/5+1,1)=norm(U3*S3*V3'-A,'fro');
        norm5(j/5+1,1)=norm(U4*S4*V4'-A,'fro');
        
        [ug,sg,vg]=svd(A);
        normg(j/5+1)=norm(ug(:,1:ranki)*sg(1:ranki,1:ranki)*vg(:,1:ranki)'-A,'fro');
        norm1(j/5+1)=norm(U*S*V'-A,'fro'); 
        
        Aold=A;
    end
end

t=[0:0.05:1];

plot(t,normg,'m-square');
hold on
plot(t,norm5,'Color','black','LineStyle','-');
hold on
plot(t,norm4,'Color','cyan','LineStyle','--','Marker','d');
hold on
plot(t,norm3,'Color','green','LineStyle','--','Marker','+');
hold on
plot(t,norm2,'b--o');
hold on
plot(t,norm1,'Color','red','LineStyle',':','Marker','*');
xlabel('t')
ylabel('error')
legend('$||X-A||$','$||Z^{8.1}-A|||$','$||Z^{8.1}-A|||$','$||Z^{9.1}-A|||$','$||Z^{9.2}-A|||$','$||Y-A||$')


%Figure 8 a ends here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 8 b
%% 

clear all
dimension=100;
tol=10^-12;
d=zeros(dimension,1);
for l=1:dimension
    d(l)=2^(-l/10);
end
D=diag(d);
T1=rand(dimension,dimension);
T2=rand(dimension,dimension);
T1=T1-diag(diag(T1));
T2=T2-diag(diag(T2));
T1=T1-T1';
T1=T1*1/5;
T2=T2-T2';
T2=T2*1/5;
Q1=eye(dimension,dimension);
Q2=Q1;
intervallength=1;
h=1/100;
A=D;
normend1=zeros(96,1);
normend2=zeros(96,1);
normend3=zeros(96,1);
normend4=zeros(96,1);
normend5=zeros(96,1);
normend6=zeros(96,1);

for ranki=3:97
ranki
[U,S,V]=svd(A);
U=U(:,1:ranki);
S=S(1:ranki,1:ranki);
V=V(:,1:ranki);
U3=U;U4=U;
S3=S;S4=S;
V3=V;V4=V;
norm1=zeros(21,1);
norm1(1)=norm(U*S*V'-A,'fro');
norm2=zeros(21,1);
norm2(1,1)=norm(U*S*V'-A,'fro');

norm3=norm2;
norm4=norm2;
norm5=norm2;
norm6=norm2;
normg=norm1;

Aold=A;

for j=1:100  
    Q1=((eye(dimension,dimension)-h*T1))\Q1;
    Q2=((eye(dimension,dimension)-h*T2))\Q2;
    
    A=Q1*(exp(j/100)*D)*Q2';
    if mod(j,5)==0
        testdelta=zeros(dimension,dimension);
        testdelta(1,:)=ones(1,dimension);
        DeltaA=A-Aold;       
        [ U,S,V ] = Algorithm10( U,S,V,DeltaA );           
        [ U1,S1,V1 ] = Algorithm9_2( A,tol,ranki );      
        [ U2,S2,V2 ] = Algorithm9_1( A,tol,ranki );       
        [ U3,S3,V3 ] = Algorithm8_2( U3,S3,V3,DeltaA,tol,ranki );       
        [ U4,S4,V4 ] = Algorithm8_1( U4,S4,V4,DeltaA,tol,ranki ); 
        
        norm2(j/5+1,1)=norm(U1*S1*V1'-A,'fro');
        norm3(j/5+1,1)=norm(U2*S2*V2'-A,'fro');
        norm4(j/5+1,1)=norm(U3*S3*V3'-A,'fro');
        norm5(j/5+1,1)=norm(U4*S4*V4'-A,'fro');
        
        [ug,sg,vg]=svd(A);
        normg(j/5+1)=norm(ug(:,1:ranki)*sg(1:ranki,1:ranki)*vg(:,1:ranki)'-A,'fro');
        norm1(j/5+1)=norm(U*S*V'-A,'fro'); 
        
        Aold=A;
    end
    
end

normend1(ranki)=normg(end);
normend2(ranki)=norm1(end);
normend3(ranki)=norm2(end);
normend4(ranki)=norm3(end);
normend5(ranki)=norm4(end);
normend6(ranki)=norm5(end);

end

semilogy(3:3:97,normend1(1:3:end-1,1),'m-square','MarkerSize',3)
hold on
semilogy(3:3:97,normend6(1:3:end-1,1),'Color','black','LineStyle','-','MarkerSize',3)
hold on
semilogy(3:3:97,normend5(1:3:end-1,1),'Color','cyan','LineStyle','--','Marker','d')
hold on
semilogy(3:3:97,normend4(1:3:end-1,1),'Color','green','LineStyle','--','Marker','+')
hold on
semilogy(3:3:97,normend3(1:3:end-1,1),'b--o','MarkerSize',3)
hold on
semilogy(3:3:97,normend2(1:3:end-1,1),'Color','red','LineStyle',':','Marker','*','MarkerSize',3)
hold on
xlabel('rank')
ylabel('error')
legend('$||X-A||$','$||Z^{8.1}-A||$','$||Z^{8.2}-A||$','$||Z^{9.1}-A||$','$||Z^{9.2}-A||$','$||Y-A||$','Interpreter','latex')
%%Figure 8b ends here