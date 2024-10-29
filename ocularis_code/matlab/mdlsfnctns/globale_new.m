function [parametri,residui]=globale_new(rif,current)

%registra current su rif, ottiene la matrice che applicata al ref port
%a al curr. Sono dati corrispondenti

ind=~isnan(current(:,1));
rif=rif(ind,:);
current=current(ind,:);

%non-linear least-squares initialization (translation)
par(1,4:6)=mean(current(:,1:3))-mean(rif(:,1:3));

%non-linear least-squares angles determination
options=optimset('lsqnonlin');
options=optimset(options,'Display','off');
rif(:,4)=1;
current(:,4)=1;
fun=inline('rif*transpose(costrrmat(par))-current','par','rif','current');

[parametri,Res1,Res2]=lsqnonlin(fun,par,-Inf,+Inf,options,rif,current);
residui=Res2(:,1:3);
