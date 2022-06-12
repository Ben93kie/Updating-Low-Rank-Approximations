%First Example that produces all the plots and the table
%It is structured such that every section produces a specific plot/table
%Run and advance runs the current section and jumps to the next

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create figure 4
clear all
dimension=100;
approxrank=10;
A1=zeros(dimension,dimension);
A2=A1;
A1(1:approxrank,1:approxrank)=eye(approxrank,approxrank)+rand(approxrank,approxrank)*0.5;
A2(1:approxrank,1:approxrank)=eye(approxrank,approxrank)+rand(approxrank,approxrank)*0.5;
A1=A1+rand(dimension,dimension)*10^-3;
A2=A2+rand(dimension,dimension)*10^-3;
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
A=A1+A2;
Aplots=zeros(dimension,dimension,3);
Aplots(:,:,1)=A;
for j=1:100
    Q1=((eye(dimension,dimension)-h*T1))\Q1;
    Q2=((eye(dimension,dimension)-h*T2))\Q2;
    A=Q1*(A1+exp(j/100)*A2)*Q2';
     if mod(j,50)==0
         Aplots(:,:,j/50+1)=A;
     end
end

for k=1:3
    figure(k)
    surf(log10(abs(Aplots(:,:,k))))
    axis([0 100 0 100 -6 2])
    hold on
end
%Create figure 4 ends here
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create figure 5

clear all
dimension=100;
approxrank=10;
tol=10^-12;
table=zeros(6,6);
A1=zeros(dimension,dimension);
A2=A1;
A1(1:approxrank,1:approxrank)=eye(approxrank,approxrank)+rand(approxrank,approxrank)*0.5;
A2(1:approxrank,1:approxrank)=eye(approxrank,approxrank)+rand(approxrank,approxrank)*0.5;
A1=A1+rand(dimension,dimension)*10^-3;
A2=A2+rand(dimension,dimension)*10^-3;
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
A=A1+A2;
[UY,SY,VY]=svd(A);
UY=UY(:,1:approxrank);
SY=SY(1:approxrank,1:approxrank);
VY=VY(:,1:approxrank);
UZ82=UY;UZ81=UY;
SZ82=SY;SZ81=SY;
VZ82=VY;VZ81=VY;
normX=zeros(21,1);
normX(1)=norm(UY*SY*VY'-A,'fro');
normZ92=normX;
normZ91=normZ92;
normZ82=normZ92;
normZ81=normZ92;
normY=normX;
Aold=A;
for j=1:100 
    Q1=((eye(dimension,dimension)-h*T1))\Q1;
    Q2=((eye(dimension,dimension)-h*T2))\Q2;     
    A=Q1*(A1+exp(j/100)*A2)*Q2';
    
    if mod(j,5)==0
       DeltaA=A-Aold;   
        [UX,SX,VX]=svd(A);
        [ UZ81,SZ81,VZ81 ] = Algorithm8_1( UZ81,SZ81,VZ81,DeltaA,tol,approxrank ); 
        [ UZ82,SZ82,VZ82 ] = Algorithm8_2( UZ82,SZ82,VZ82,DeltaA,tol,approxrank ); 
        [ UZ91,SZ91,VZ91 ] = Algorithm9_1( A,tol,approxrank );    
        [ UZ92,SZ92,VZ92 ] = Algorithm9_2( A,tol,approxrank );   
        [ UY,SY,VY ]       = Algorithm10 ( UY,SY,VY,DeltaA );        
        
        normX(j/5+1)=norm(UX(:,1:approxrank)*SX(1:approxrank,1:approxrank)*VX(:,1:approxrank)'-A,'fro');
        normZ81(j/5+1)=norm(UZ81*SZ81*VZ81'-A,'fro');
        normZ82(j/5+1)=norm(UZ82*SZ82*VZ82'-A,'fro');
        normZ91(j/5+1)=norm(UZ91*SZ91*VZ91'-A,'fro');
        normZ92(j/5+1)=norm(UZ92*SZ92*VZ92'-A,'fro');
        normY(j/5+1)=norm(UY*SY*VY'-A,'fro');
               
        Aold=A;
    end
    
end
t=[0:0.05:1];
plot(t,normX,'m-square');
hold on
plot(t,normZ81,'Color','black','LineStyle','-');
hold on
plot(t,normZ82,'Color','cyan','LineStyle','--','Marker','d');
hold on
plot(t,normZ91,'Color','green','LineStyle','--','Marker','+');
hold on
plot(t,normZ92,'b--o');
hold on
plot(t,normY,'Color','red','LineStyle',':','Marker','*');
xlabel('t')
ylabel('error')
legend('$||X-A||$','$||Z^{8.1}-A||$','$||Z^{8.2}-A||$','$||Z^{9.1}-A||$','$||Z^{9.2}-A||$','$||Y-A||$','Interpreter','latex')

%Create figure 5 ends here
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create table 2
clear all
epsvec=[10^-1,10^-2,10^-3,10^-4,10^-5,0];
table=zeros(6,7);
for k=1:length(epsvec)
dimension=100;
approxrank=10;
tol=10^-12;
A1=zeros(dimension,dimension);
A2=A1;
A1(1:approxrank,1:approxrank)=eye(approxrank,approxrank)+rand(approxrank,approxrank)*0.5;
A2(1:approxrank,1:approxrank)=eye(approxrank,approxrank)+rand(approxrank,approxrank)*0.5;
A1=A1+rand(dimension,dimension)*epsvec(k);
A2=A2+rand(dimension,dimension)*epsvec(k);
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
A=A1+A2;
[UY,SY,VY]=svd(A);
UY=UY(:,1:approxrank);
SY=SY(1:approxrank,1:approxrank);
VY=VY(:,1:approxrank);
UZ82=UY;UZ81=UY;
SZ82=SY;SZ81=SY;
VZ82=VY;VZ81=VY;
normX=zeros(21,1);
normX(1)=norm(UY*SY*VY'-A,'fro');
normZ92=normX;
normZ91=normZ92;
normZ82=normZ92;
normZ81=normZ92;
normY=normX;
Aold=A;
for j=1:100 
    Q1=((eye(dimension,dimension)-h*T1))\Q1;
    Q2=((eye(dimension,dimension)-h*T2))\Q2;     
    A=Q1*(A1+exp(j/100)*A2)*Q2';
    
    if mod(j,5)==0
       DeltaA=A-Aold;   
        [UX,SX,VX]=svd(A);
        [ UZ81,SZ81,VZ81 ] = Algorithm8_1( UZ81,SZ81,VZ81,DeltaA,tol,approxrank ); 
        [ UZ82,SZ82,VZ82 ] = Algorithm8_2( UZ82,SZ82,VZ82,DeltaA,tol,approxrank ); 
        [ UZ91,SZ91,VZ91 ] = Algorithm9_1( A,tol,approxrank );    
        [ UZ92,SZ92,VZ92 ] = Algorithm9_2( A,tol,approxrank );   
        [ UY,SY,VY ]       = Algorithm10 ( UY,SY,VY,DeltaA );        
        
        normX(j/5+1)=norm(UX(:,1:approxrank)*SX(1:approxrank,1:approxrank)*VX(:,1:approxrank)'-A,'fro');
        normZ81(j/5+1)=norm(UZ81*SZ81*VZ81'-A,'fro');
        normZ82(j/5+1)=norm(UZ82*SZ82*VZ82'-A,'fro');
        normZ91(j/5+1)=norm(UZ91*SZ91*VZ91'-A,'fro');
        normZ92(j/5+1)=norm(UZ92*SZ92*VZ92'-A,'fro');
        normY(j/5+1)=norm(UY*SY*VY'-A,'fro');
               
        Aold=A;
    end    
end
table(k,1)=epsvec(k);
table(k,2)=normX(end);
table(k,3)=normZ81(end);
table(k,4)=normZ82(end);
table(k,5)=normZ91(end);
table(k,6)=normZ92(end);
table(k,7)=normY(end);
end
table

%Create table 2 ends here