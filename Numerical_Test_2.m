%%Figure 6 and 7. See below
clear all
dimension=100;
ranki=10;
approxrank1=5;
approxrank2=20;
tol=10^-12;
averaging=1;
A1=zeros(dimension,dimension);
A2=A1;
A1(1:ranki,1:ranki)=eye(ranki,ranki)+rand(ranki,ranki)*0.5;
A2(1:ranki,1:ranki)=eye(ranki,ranki)+rand(ranki,ranki)*0.5;
A1=A1+rand(dimension,dimension)*10^-2;
A2=A2+rand(dimension,dimension)*10^-2;
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
intervallength=10;
h=1/100;
interval=[0:h:intervallength];
intervalveclength=length(interval);
repeattimes=20;
everyxthstep=intervalveclength/20;
A=A1+A2;
[U,S,V]=svd(A);
UY1=U(:,1:approxrank1);
SY1=S(1:approxrank1,1:approxrank1);
VY1=V(:,1:approxrank1);
UY2=U(:,1:approxrank1);
SY2=S(1:approxrank1,1:approxrank1);
VY2=V(:,1:approxrank1);
UY3=U(:,1:20);
SY3=S(1:20,1:20);
VY3=V(:,1:20);
UZ815=U(:,1:approxrank1);
UZ81=U(:,1:10);
UZ81re=U(:,1:5);
SZ815=S(1:5,1:5);
SZ81=S(1:10,1:10);
SZ81re=S(1:5,1:5);
VZ815=V(:,1:5);
VZ81=V(:,1:10);
VZ81re=V(:,1:5);
UZ825=U(:,1:approxrank1);
SZ825=S(1:approxrank1,1:approxrank1);
VZ825=V(:,1:approxrank1);
UZ82=U(:,1:10);
SZ82=S(1:10,1:10);
VZ82=V(:,1:10);
UZ82re=UY1;
SZ82re=SY1;
VZ82re=VY1;

singvecsexact=zeros(intervalveclength,approxrank2);
singvecY1=zeros(intervalveclength,approxrank1);
singvecY2=zeros(intervalveclength,approxrank1);
singvecZ91=zeros(intervalveclength,approxrank1);
singvecZ9110=zeros(intervalveclength,10);
singvecZ92=zeros(intervalveclength,approxrank1);
singvecZ9210=zeros(intervalveclength,10);
singvecZ825=zeros(intervalveclength,approxrank1);
singvecZ82=zeros(intervalveclength,10);
singvecZ82re=zeros(intervalveclength,approxrank1);
singvecZ815=zeros(intervalveclength,5);
singvecZ81=zeros(intervalveclength,10);
singvecZ81re=zeros(intervalveclength,5);

singvecsexact(1,:)=diag(S(1:approxrank2,1:approxrank2));
singvecY1(1,:)=singvecsexact(1,1:approxrank1);
singvecY2(1,:)=singvecsexact(1,1:approxrank1);
singvecY3=zeros(intervalveclength,20);
singvecY3(1,:)=diag(S(1:20,1:20));
singvecZ91(1,:)=singvecsexact(1,1:approxrank1);
singvecZ9110(1,:)=singvecsexact(1,1:10);
singvecZ92(1,:)=singvecsexact(1,1:approxrank1);
singvecZ9210(1,:)=singvecsexact(1,1:10);
singvecZ825(1,:)=singvecsexact(1,1:approxrank1);
singvecZ82(1,:)=singvecsexact(1,1:10);
singvecZ82re(1,:)=singvecsexact(1,1:approxrank1);
singvecZ815(1,:)=singvecsexact(1,1:5);
singvecZ81(1,:)=singvecsexact(1,1:10);
singvecZ81re(1,:)=singvecsexact(1,1:5);
Aold=A;
normX=zeros(1001,1);
normX10=normX;
normX(1)=norm(UY1*SY1*VY1'-A,'fro');
normX10(1)=norm(U(:,1:10)*S(1:10,1:10)*V(:,1:10)'-A,'fro');
normX20(1)=norm(U(:,1:20)*S(1:20,1:20)*V(:,1:20)'-A,'fro');
normY1=normX;
normY2=normX;
normY3=normX20;
normZ81=normX10;
normZ815=normX;
normZ81re=normX;
normZ82=normX10;
normZ825=normX;
normZ82re=normX;
normZ91=normX;
normZ9110=normX10;
normZ92=normX;
normZ9210=normX10;



for j=1:intervalveclength-1
   
    Q1=((eye(dimension,dimension)-h*T1))\Q1;
    Q2=((eye(dimension,dimension)-h*T2))\Q2;
    j
    A=Q1*(A1+cos(j*h)*A2)*Q2';
    DeltaA=A-Aold;    
    [UX,SX,VX]=svd(A);
   [ UY1,SY1,VY1 ]      = Algorithm10( UY1,SY1,VY1,DeltaA );   
   if mod(j,30)==0
            UY2=UX(:,1:5);
            SY2=SX(1:5,1:5);
            VY2=VX(:,1:5);
   end
   [ UY2,SY2,VY2 ]      = Algorithm10( UY2,SY2,VY2,DeltaA );
   [ UY3,SY3,VY3 ]      = Algorithm10( UY3,SY3,VY3,DeltaA ); 
   [ UZ91,SZ91,VZ91 ]   = Algorithm9_1( A,tol,approxrank1 );
   [ UZ9110,SZ9110,VZ9110 ]   = Algorithm9_1( A,tol,10 );
   [ UZ92,SZ92,VZ92 ]   = Algorithm9_2( A,tol,approxrank1 );
   [ UZ9210,SZ9210,VZ9210 ]   = Algorithm9_2( A,tol,10 );
   [ UZ815,SZ815,VZ815 ]= Algorithm8_1( UZ815,SZ815,VZ815,DeltaA,tol,approxrank1 );
   [ UZ81,SZ81,VZ81 ]= Algorithm8_1( UZ81,SZ81,VZ81,DeltaA,tol,10 );
   if mod(j,30)==0
            UZ81re=UX(:,1:5);
            SZ81re=SX(1:5,1:5);
            VZ81re=VX(:,1:5);
   end
   [ UZ81re,SZ81re,VZ81re ]= Algorithm8_1( UZ81re,SZ81re,VZ81re,DeltaA,tol,5 );
   [ UZ825,SZ825,VZ825 ]= Algorithm8_2( UZ825,SZ825,VZ825,DeltaA,tol,approxrank1 );
   [ UZ82,SZ82,VZ82 ]   = Algorithm8_2( UZ82,SZ82,VZ82,DeltaA,tol,10 );
   if mod(j,30)==0
            UZ82re=UX(:,1:5);
            SZ82re=SX(1:5,1:5);
            VZ82re=VX(:,1:5);
   end
   [ UZ82re,SZ82re,VZ82re ]   = Algorithm8_2( UZ82re,SZ82re,VZ82re,DeltaA,tol,approxrank1 );
   singY1=svd(SY1);
   singY1=singY1(1:5);
   singY2=svd(SY2);
   singY3=svd(SY3);
   singY3=singY3(1:20);
   singZ91=diag(SZ91(1:5,1:5)); 
   singZ9110=diag(SZ9110(1:10,1:10)); 
   singZ92=diag(SZ92(1:5,1:5)); 
   singZ9210=diag(SZ9210(1:10,1:10)); 
   singZ815=diag(SZ815(1:5,1:5));
   singZ81=diag(SZ81(1:10,1:10));
   singZ81re=diag(SZ81re(1:5,1:5));
   singZ825=diag(SZ825(1:5,1:5));
   singZ82=diag(SZ82(1:10,1:10));
   singZ82re=diag(SZ82re(1:5,1:5));
        
   singvecsexact(j+1,:)=singvecsexact(j+1,:)+diag(SX(1:approxrank2,1:approxrank2))';
   singvecY1(j+1,:)=singvecY1(j+1,:)+singY1';
   singvecY2(j+1,:)=singvecY2(j+1,:)+singY2(1:approxrank1)';
   singvecY3(j+1,:)=singvecY3(j+1,:)+singY3(1:20)';
   singvecZ91(j+1,:)=singvecZ91(j+1,:)+singZ91(1:approxrank1)';
   singvecZ9110(j+1,:)=singvecZ9110(j+1,:)+singZ9110(1:10)';
   singvecZ92(j+1,:)=singvecZ92(j+1,:)+singZ92(1:approxrank1)';
   singvecZ9210(j+1,:)=singvecZ9210(j+1,:)+singZ9210(1:10)';
   singvecZ815(j+1,:)=singvecZ815(j+1,:)+singZ815(1:approxrank1)';
   singvecZ81(j+1,:)=singvecZ81(j+1,:)+singZ81(1:10)';
   singvecZ81re(j+1,:)=singvecZ81re(j+1,:)+singZ81re(1:approxrank1)';
   singvecZ825(j+1,:)=singvecZ825(j+1,:)+singZ825(1:approxrank1)';
   singvecZ82(j+1,:)=singvecZ82(j+1,:)+singZ82(1:10)';
   singvecZ82re(j+1,:)=singvecZ82re(j+1,:)+singZ82re(1:approxrank1)';
   Aold=A;
   
   normX(j+1)       =norm(UX(:,1:5)*SX(1:5,1:5)*VX(:,1:5)'-A,'fro');
   normX10(j+1)     =norm(UX(:,1:10)*SX(1:10,1:10)*VX(:,1:10)'-A,'fro');
   normY1(j+1)      =norm(UY1*SY1*VY1'-A,'fro');
   normY2(j+1)      =norm(UY2*SY2*VY2'-A,'fro');
   normY3(j+1)      =norm(UY3*SY3*VY3'-A,'fro');
   normZ81(j+1)     =norm(UZ81*SZ81*VZ81'-A,'fro');
   normZ815(j+1)    =norm(UZ815*SZ815*VZ815'-A,'fro');
   normZ81re(j+1)   =norm(UZ81re*SZ81re*VZ81re'-A,'fro');
   normZ82(j+1)     =norm(UZ82*SZ82*VZ82'-A,'fro');
   normZ825(j+1)    =norm(UZ825*SZ825*VZ825'-A,'fro');
   normZ82re(j+1)   =norm(UZ82re*SZ82re*VZ82re'-A,'fro');
   normZ91(j+1)     =norm(UZ91*SZ91*VZ91'-A,'fro');
   normZ9110(j+1)   =norm(UZ9110*SZ9110*VZ9110'-A,'fro');
   normZ92(j+1)     =norm(UZ92*SZ92*VZ92'-A,'fro');
   normZ9210(j+1)   =norm(UZ9210*SZ9210*VZ9210'-A,'fro');
end
%%%%%%%
%Figure 6
t=[0:0.1:10];
figure(14)
plot(t,normX(1:10:end),'m-square')
hold on
plot(t,normZ815(1:10:end),'Color','black','LineStyle','-')
hold on
plot(t,normZ825(1:10:end),'Color','cyan','LineStyle','--','Marker','d')
hold on
plot(t, normZ91(1:10:end),'Color','green','LineStyle','--','Marker','+')
hold on
plot(t,normZ92(1:10:end),'b--o')
hold on
plot(t,normY1(1:10:end),'Color','red','LineStyle',':','Marker','*')
hold on
xlabel('t')
legend('$||X-A||$','$||Z^{8.1}-A||$','$||Z^{8.2}-A||$','$||Z^{9.1}-A||$','$||Z^{9.2}-A||$','$||Y-A||$','Interpreter','latex')

figure(15)
plot(t,normX(1:10:end),'m-square')
hold on
plot(t,normZ81re(1:10:end),'Color','black','LineStyle','-')
hold on
plot(t,normZ82re(1:10:end),'Color','cyan','LineStyle','--','Marker','d')
hold on
plot(t, normZ91(1:10:end),'Color','green','LineStyle','--','Marker','+')
hold on
plot(t,normZ92(1:10:end),'b--o')
hold on
plot(t,normY2(1:10:end),'Color','red','LineStyle',':','Marker','*')
hold on
xlabel('t')
legend('$||X-A||$','$||Z^{8.1}-A||$','$||Z^{8.2}-A||$','$||Z^{9.1}-A||$','$||Z^{9.2}-A||$','$||Y-A||$','Interpreter','latex')

figure(16)
plot(t,normX10(1:10:end),'m-square')
hold on
plot(t,normZ81(1:10:end),'Color','black','LineStyle','-')
hold on
plot(t,normZ82(1:10:end),'Color','cyan','LineStyle','--','Marker','d')
hold on
plot(t, normZ9110(1:10:end),'Color','green','LineStyle','--','Marker','+')
hold on
plot(t,normZ9210(1:10:end),'b--o')
hold on
plot(t,normY3(1:10:end),'Color','red','LineStyle',':','Marker','*')
hold on
xlabel('t')
legend('$||X-A||$','$||Z^{8.1}-A||$','$||Z^{8.2}-A||$','$||Z^{9.1}-A||$','$||Z^{9.2}-A||$','$||Y-A||$','Interpreter','latex')

%Figure 6 ends here
%%%%%%%%%%%%%%%%%%%%
%Figure 7
%% 

figure(1)
for f=1:5
plot(interval,singvecsexact(1:end,f),'m')
hold on
plot(interval(1:10:end),singvecY1(1:10:end,f),'-ro','MarkerSize',1.5)
hold on
end
%xlabel('t')
%ylabel('Magnitude of Singular Values')
legend('$\sigma_j$, $j=1,\dots,5$ of $X$','$\sigma_j$, $j=1,\dots,5$ of rank-5 $Y$')

figure(2)
for f=1:5
plot(interval,singvecsexact(1:end,f),'m')
hold on
plot(interval(1:10:end),singvecZ92(1:10:end,f),'-bo','MarkerSize',1.5)
hold on
end
%xlabel('t')
%ylabel('Magnitude of Singular Values')
legend('$\sigma_j$, $j=1,\dots,5$ of $X$','$\sigma_j$, $j=1,\dots,5$ of rank-5 $Z^{9.2}$')

figure(3)
for f=1:5
plot(interval,singvecsexact(1:end,f),'m')
hold on
plot(interval(1:10:end),singvecY2(1:10:end,f),'-ro','MarkerSize',1.5)
hold on
end
%xlabel('t')
%ylabel('Magnitude of Singular Values')
legend('$\sigma_j$, $j=1,\dots,5$ of $X$','$\sigma_j$, $j=1,\dots,5$ of rank-5 $Y$ with restart')

figure(4)
for f=1:20
plot(interval,singvecsexact(1:end,f),'m')
hold on
plot(interval(1:10:end),singvecY3(1:10:end,f),'-ro','MarkerSize',2)
hold on
end    
%xlabel('t')
%ylabel('Magnitude of Singular Values')
legend('$\sigma_j$, $j=1,\dots,20$ of $X$','$\sigma_j$, $j=1,\dots,20$ of rank-20 $Y$')

figure(5)
for f=1:10
plot(interval,singvecsexact(1:end,f),'m')
hold on
plot(interval(1:10:end),singvecZ82(1:10:end,f),'-co','MarkerSize',2)
hold on
end    
%xlabel('t')
%ylabel('Magnitude of Singular Values')
legend('$\sigma_j$, $j=1,\dots,10$ of $X$','$\sigma_j$, $j=1,\dots,10$ of rank-10 $Z^{8.2}$')

figure(6)
for f=1:5
plot(interval,singvecsexact(1:end,f),'m')
hold on
plot(interval(1:10:end),singvecZ82re(1:10:end,f),'-co','MarkerSize',2)
hold on
end    
%xlabel('t')
%ylabel('Magnitude of Singular Values')
legend('$\sigma_j$, $j=1,\dots,5$ of $X$','$\sigma_j$, $j=1,\dots,5$ of rank-5 $Z^{8.2}$ with restart')

figure(7)
for f=1:5
plot(interval,singvecsexact(1:end,f),'m')
hold on
plot(interval(1:10:end),singvecZ825(1:10:end,f),'-co','MarkerSize',2)
hold on
end    
%xlabel('t')
%ylabel('Magnitude of Singular Values')
legend('$\sigma_j$, $j=1,\dots,5$ of $X$','$\sigma_j$, $j=1,\dots,5$ of rank-5 $Z^{8.2}$')

figure(8)
for f=1:10
plot(interval,singvecsexact(1:end,f),'m')
hold on
plot(interval(1:10:end),singvecZ81(1:10:end,f),'-ko','MarkerSize',2)
hold on
end    
%xlabel('t')
%ylabel('Magnitude of Singular Values')
legend('$\sigma_j$, $j=1,\dots,10$ of $X$','$\sigma_j$, $j=1,\dots,10$ of rank-10 $Z^{8.1}$')

figure(9)
for f=1:5
plot(interval,singvecsexact(1:end,f),'m')
hold on
plot(interval(1:10:end),singvecZ81re(1:10:end,f),'-ko','MarkerSize',2)
hold on
end    
%xlabel('t')
%ylabel('Magnitude of Singular Values')
legend('$\sigma_j$, $j=1,\dots,5$ of $X$','$\sigma_j$, $j=1,\dots,5$ of rank-5 $Z^{8.1}$ with restart')

figure(10)
for f=1:5
plot(interval,singvecsexact(1:end,f),'m')
hold on
plot(interval(1:10:end),singvecZ815(1:10:end,f),'-ko','MarkerSize',2)
hold on
end    
%xlabel('t')
%ylabel('Magnitude of Singular Values')
legend('$\sigma_j$, $j=1,\dots,5$ of $X$','$\sigma_j$, $j=1,\dots,5$ of rank-5 $Z^{8.1}$')

figure(11)
for f=1:5
plot(interval,singvecsexact(1:end,f),'m')
hold on
plot(interval(1:10:end),singvecZ91(1:10:end,f),'-go','MarkerSize',1.5)
hold on
end
%xlabel('t')
%ylabel('Magnitude of Singular Values')
legend('$\sigma_j$, $j=1,\dots,5$ of $X$','$\sigma_j$, $j=1,\dots,5$ of rank-5 $Z^{9.1}$')

figure(12)
for f=1:10
plot(interval,singvecsexact(1:end,f),'m')
hold on
plot(interval(1:10:end),singvecZ9110(1:10:end,f),'-go','MarkerSize',1.5)
hold on
end
%xlabel('t')
%ylabel('Magnitude of Singular Values')
legend('$\sigma_j$, $j=1,\dots,10$ of $X$','$\sigma_j$, $j=1,\dots,10$ of rank-10 $Z^{9.1}$')

figure(13)
for f=1:10
plot(interval,singvecsexact(1:end,f),'m')
hold on
plot(interval(1:10:end),singvecZ9210(1:10:end,f),'-bo','MarkerSize',1.5)
hold on
end
%xlabel('t')
%ylabel('Magnitude of Singular Values')
legend('$\sigma_j$, $j=1,\dots,10$ of $X$','$\sigma_j$, $j=1,\dots,10$ of rank-10 $Z^{9.2}$')

