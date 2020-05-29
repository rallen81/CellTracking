function x=mypdfcos(sigmar,mur,sigmatheta,thetazero,r,cos)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

a=-1./(cos);
kboundup=abs(ceil(thetazero+sigmatheta*3));
kbound=abs(ceil(thetazero-sigmatheta*3));
Xt=0;

if sigmatheta<3.5
for k=-kbound:kboundup
Xt=Xt+a.*(mynormal(-acos(cos)+2*pi*k,sigmatheta,thetazero)+mynormal(acos(cos)+2*pi*k,sigmatheta,thetazero));
end
end



%work out normal with extra precision. Scale? High resolution standard
%normal?

if sigmatheta>=3.5
Xt=a/pi;   
end

Rt=mynormal(r,sigmar,mur);
x=Xt.*Rt;


%if abs(r)<0.5
%Rt=1;
%else
%Rt=0;    
%end    

%if abs(cos)<0.5
%Xt=1;
%else
%Xt=0;    ZX
%end    

%x=Xt*Rt;
end

