%Plot 9 a. 20 minutes computation time
load('docbyterm.tfidf.norm.txt')
tol=10^-12;
A=zeros(7095,5896);
for j=1:247159
    A(docbyterm_tfidf_norm(j,1),docbyterm_tfidf_norm(j,2))=docbyterm_tfidf_norm(j,3);
end
timeg=zeros(191,1);
timed=timeg;
Asparse=sparse(A);
for ranki=30:30
    ranki
[Ud,Sd,Vd]=svds(Asparse,84);
Ua=Ud;
Sa=Sd;
Va=Vd;
Uas=Ua;
Sas=Sd;
Vas=Va;
Ug=Ud;
Sg=Sd;
Vg=Vd;
normg=zeros(11,1);
normg(1)=norm(Ud*Sd*Vd'-Asparse,'fro')/norm(Asparse,'fro');
normd=normg;
normt=normg;
norma=normg;
normas=normg;
normst=norma;
for k=2:11
    k
    DeltaA=sprand(7095,5896,2.4027e-04);
    Asparse=Asparse+DeltaA;
    A=A+DeltaA;
    [ Ug,Sg,Vg ] = svds( Asparse,84);
    [ Uas,Sas,Vas ] = Algorithm8_1( Uas,Sas,Vas,DeltaA,tol,84 );
    [ Ua,Sa,Va ] = Algorithm8_2( Ua,Sa,Va,DeltaA,tol,84 );
    [ Ust,Sst,Vst ] = Algorithm9_1( A,tol,84 );
    [ Ut,St,Vt ] = Algorithm9_2( A,10^-9,84 );
    [ Ud,Sd,Vd ] = Algorithm10( Ud,Sd,Vd,DeltaA);
     
     normd(k)=norm(Ud*Sd*Vd'-Asparse,'fro')/norm(Asparse,'fro');
     normg(k)=norm(Ug*Sg*Vg'-Asparse,'fro')/norm(Asparse,'fro');
     normt(k)=norm(Ut*St*Vt'-Asparse,'fro')/norm(Asparse,'fro');
     norma(k)=norm(Ua*Sa*Va'-Asparse,'fro')/norm(Asparse,'fro');
     normas(k)=norm(Uas*Sas*Vas'-Asparse,'fro')/norm(Asparse,'fro');
     normst(k)=norm(Ust*Sst*Vst'-Asparse,'fro')/norm(Asparse,'fro');
    
end
end


plot([0:10],normg(1:11),'m-square')
hold on
plot([0:10],normas(1:11),'Color','black','LineStyle','-')
hold on
plot([0:10],norma(1:11),'Color','cyan','LineStyle','--','Marker','d')
hold on
plot([0:10],norma(1:11),'Color','green','LineStyle','--','Marker','+')
hold on
plot([0:10],normt(1:11),'b--o')
hold on
plot([0:10],normd(1:11),'Color','red','LineStyle',':','Marker','*')
xlabel('t')
ylabel('error')
legend('$\frac{||X-A||}{||A||}$','$\frac{||Z^{8.1}-A||}{||A||}$','$\frac{||Z^{8.2}-A||}{||A||}$','$\frac{||Z^{9.1}-A||}{||A||}$','$\frac{||Z^{9.2}-A||}{||A||}$','$\frac{||Y-A||}{||A||}$','Interpreter','latex')


%Plot 9 a ends here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot 9 b
%% Plot takes long (1 hour)

load('docbyterm.tfidf.norm.txt')
tol=10^-12;
A=zeros(7095,5896);
for j=1:247159
    A(docbyterm_tfidf_norm(j,1),docbyterm_tfidf_norm(j,2))=docbyterm_tfidf_norm(j,3);
end
timeg=zeros(101,1);
timed=timeg;
Asparse=sparse(A);
[Ud,Sd,Vd]=svds(Asparse,84);
Ug=Ud;
Sg=Sd;
Vg=Vd;
for ranki=30:130
    ranki

for k=2:11
    
    DeltaA=sprand(7095,5896,2.4027e-04);
    Asparse=Asparse+DeltaA;
    A=A+DeltaA;
    tic
    [ Ug,Sg,Vg ] = svds( Asparse,ranki);
    timeg(ranki-29)=timed(ranki-29)+toc;
    tic
    [ Ud,Sd,Vd ] = Algorithm10( Ug,Sg,Vg,DeltaA);
    timed(ranki-29)=timed(ranki-29)+toc;
    
    
end
end

plot([30:130],timed,'Color','red','LineStyle',':','Marker','*')
hold on
plot([30:130],timeg,'m-square')
xlabel('rank')
ylabel('CPU time')
legend('$||X-A||$','$||Y-A||$','Interpreter','latex')


