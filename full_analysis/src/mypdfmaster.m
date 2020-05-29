function [xtotal,xtestout]=mypdfmaster(sigmar,mur,sigmatheta,thetazero,xtest)


epsilon=0.005; %0.05
dL=0.01; %0.01
deltaL=log(epsilon):-dL:-10;
dthet=0.00005; %% change this for better acc. %0.01 
xtotal=zeros(numel(xtest),1);
deltatheta1=acos(-1):-dthet:acos(-epsilon);
deltatheta2=acos(epsilon):-dthet:acos(1);
lower=-0.01; %change these to smaller if you have a tight singularity at 0. Can be zero. 
higher=-lower;

for i=1:numel(xtest)
    
if xtest(i)~=0    

x1=mypdflog(sigmar,mur,sigmatheta,thetazero,xtest(i)./(exp(deltaL)),exp(deltaL));
x2=mypdflog(sigmar,mur,sigmatheta,thetazero,-xtest(i)./(exp(deltaL)),-exp(deltaL));
x3=-mypdfcos(sigmar,mur,sigmatheta,thetazero,xtest(i)./cos(deltatheta1),cos(deltatheta1));
x4=mypdfcos(sigmar,mur,sigmatheta,thetazero,xtest(i)./cos(deltatheta2),cos(deltatheta2));

xtotal(i)=dL*(sum(x1)+sum(x2))-dthet*(sum(x3)+sum(x4));
%sum(x1)+sum(x2)
else

    
x1=mypdflog(sigmar,mur,sigmatheta,thetazero,lower./(exp(deltaL)),exp(deltaL));
x2=mypdflog(sigmar,mur,sigmatheta,thetazero,-lower./(exp(deltaL)),-exp(deltaL));
x3=-mypdfcos(sigmar,mur,sigmatheta,thetazero,lower./cos(deltatheta1),cos(deltatheta1));
x4=mypdfcos(sigmar,mur,sigmatheta,thetazero,lower./cos(deltatheta2),cos(deltatheta2));

xtotallower=dL*(sum(x1)+sum(x2))-dthet*(sum(x3)+sum(x4));

x1=mypdflog(sigmar,mur,sigmatheta,thetazero,higher./(exp(deltaL)),exp(deltaL));
x2=mypdflog(sigmar,mur,sigmatheta,thetazero,-higher./(exp(deltaL)),-exp(deltaL));
x3=-mypdfcos(sigmar,mur,sigmatheta,thetazero,higher./cos(deltatheta1),cos(deltatheta1));
x4=mypdfcos(sigmar,mur,sigmatheta,thetazero,higher./cos(deltatheta2),cos(deltatheta2));

xtotalhigher=dL*(sum(x1)+sum(x2))-dthet*(sum(x3)+sum(x4));
    
xtotal(i)=0.5*(xtotalhigher+xtotallower);
end


xtestout=xtest;

end