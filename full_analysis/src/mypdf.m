function x=mypdf(sigmar,mur,sigmatheta,thetazero,r,cos)

a=1./(cos.*(1-cos.^2).^0.5);
kboundup=abs(ceil(thetazero+sigmatheta*3));
kbound=abs(ceil(thetazero-sigmatheta*3));
Xt=0;

if sigmatheta<3.5
    for k=-kbound:kboundup
        Xt=Xt+(mynormal(-acos(cos)+2*pi*k,sigmatheta,thetazero)+mynormal(acos(cos)+2*pi*k,sigmatheta,thetazero));
    end
Xt=Xt.*a;

end

if sigmatheta>=3.5
Xt=a/pi;   
end
Rt=mynormal(r,sigmar,mur);
x=Xt.*Rt;


end

