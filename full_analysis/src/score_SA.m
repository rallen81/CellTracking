function [score,xout,Xf]=score_SA(p,xFitted,nFit,yFitted,nFity,int_acc,yflag)
    
alpha=p(1);
mur=p(2);
sigmar=p(3);
sigmatheta=p(4);
mur2=p(5);
sigmar2=p(6);
sigmatheta2=p(7);

%int_acc=1000;
%%change the integration here to match the SA in the function file
xout=xFitted;
nx=nFit;
yout=yFitted;
ny=nFity;

Deltaxout=0.5*(xout(2)-xout(1));
t=1;

%[ny,yout]=scalehist(cellnewx(:,2),lower);

xmin=min(abs(xout));
imintemp=find(abs(xout)==xmin);

imin=imintemp(1);

nx(imin)=nx(imin)*(xout(imin)-xout(imin-1));

[Xf,x]=mysprocess(sigmar,mur,sigmatheta,0,t,xout,int_acc);
[Xf2,x2]=mysprocess(sigmar2,mur2,sigmatheta2,0,t,xout,int_acc);

[Xfplus,x]=mysprocess(sigmar,mur,sigmatheta,0,t,xout(imin)+Deltaxout,int_acc);
[Xf2plus,x2]=mysprocess(sigmar2(1),mur2(1),sigmatheta2(1),0,t,xout(imin)+Deltaxout,int_acc);
[Xfminus,x]=mysprocess(sigmar(1),mur(1),sigmatheta(1),0,t,xout(imin)-Deltaxout,int_acc);
[Xf2minus,x2]=mysprocess(sigmar2(1),mur2(1),sigmatheta2(1),0,t,xout(imin)-Deltaxout,int_acc);


sumX=Xf*alpha(1)+Xf2*(1-alpha(1));
sumXplus=Xfplus*alpha(1)+Xf2plus*(1-alpha(1));
sumXminus=Xfminus*alpha(1)+Xf2minus*(1-alpha(1));

sumX(imin)=0.5*(sumX(imin)+sumXplus)*Deltaxout+0.5*(sumX(imin)+sumXminus)*Deltaxout;
scoreOut=nx-sumX;
scoreOutX=scoreOut.^2;


if yflag==1
%%%% Fit the Y-histogram too. 


Deltayout=0.5*(yout(2)-yout(1));


ymin=min(abs(yout));
imintemp=find(abs(yout)==ymin);

iminy=imintemp(1);

ny(iminy)=ny(iminy)*(yout(iminy)-yout(iminy-1));

[Yf,y]=mysprocess(sigmar(1),mur(1),sigmatheta(1),pi/2,t,yout);
[Yf2,y2]=mysprocess(sigmar2(1),mur2(1),sigmatheta2(1),pi/2,t,yout);

[Yf,y]=mysprocess(sigmar2(1),mur2(1),sigmatheta2(1),pi/2,t,yout);
[Yf2,y2]=mysprocess(sigmar2(1),mur2(1),sigmatheta2(1),pi/2,t,yout);

[Yfplus,y]=mysprocess(sigmar(1),mur(1),sigmatheta(1),pi/2,t,yout(iminy)+Deltayout);
[Yf2plus,y2]=mysprocess(sigmar2(1),mur2(1),sigmatheta2(1),pi/2,t,yout(iminy)+Deltayout);
[Yfminus,y]=mysprocess(sigmar(1),mur(1),sigmatheta(1),pi/2,t,yout(iminy)-Deltayout);
[Yf2minus,y2]=mysprocess(sigmar2(1),mur2(1),sigmatheta2(1),pi/2,t,yout(iminy)-Deltayout);

sumY=Yf*alpha(1)+Yf2*(1-alpha(1));
sumYplus=Yfplus*alpha(1)+Yf2plus*(1-alpha(1));
sumYminus=Yfminus*alpha(1)+Yf2minus*(1-alpha(1));

sumY(iminy)=0.5*(sumY(iminy)+sumYplus)*Deltayout+0.5*(sumY(iminy)+sumYminus)*Deltayout;

scoreOutY=ny-sumY;

scoreOutY=scoreOutY.^2;
%}
score=sum(scoreOutX)+sum(scoreOutY);
elseif yflag==0
    score=sum(scoreOutX);
end

end