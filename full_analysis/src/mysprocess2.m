function [Xf,x,zero]=mysprocess2(sigmar,mur,sigmatheta,thetazero,t,x,varargin)


c=0;
if isempty(varargin)
nJ=120;
else
nJ=varargin{1};    
end


nC=numel(x);
deltaS=1/nJ;


Xmax=(mur+5*sigmar); %5 used for MC runs
Xmin=-Xmax;
deltaC=(Xmax-Xmin)/nC;
dtheta=0.1;  %0.01 used for MC runs irrelevant if t=1
thta=[-3.14:dtheta:3.14];

zero=0;
%x=[Xmin:deltaC:Xmax];
xInt=zeros(1,nC);
Xf=zeros(1,nC);
xTotal=Xf;


S=[deltaS:deltaS:1-deltaS];
F=1./S;
if t==1
      theta=thetazero;
      thF=1;
else
   thF=numel(thta); 
end
      
      
for th=1:thF
    if t~=1
    theta=thta(th);
    
 % numel(thta)-th;
    end
  
  
for i=1:numel(x)
 
 
C=x(i);
S2=C./S;


if abs(C)<0.5

[xtotal,xtestout]=mypdfmaster(sigmar,mur,sigmatheta,thetazero,C);
xInt(i)=xtotal;
else

A=mypdf(sigmar,mur,sigmatheta,theta,S2,S);
B=mypdf(sigmar,mur,sigmatheta,theta,-S2,-S);

xInt(i)=sum(A-B)*deltaS;

end

end


if t~=1
xTotal=xTotal+xInt.*mynormal(theta,sigmatheta*(t-1)^0.5,thetazero)*dtheta;
else
xTotal=xInt;    
end
end

Xf=xTotal;


end

 