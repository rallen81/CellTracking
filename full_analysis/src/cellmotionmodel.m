function [n,xout,A,controldata,coord,state,cellstepstotal]=cellmotionmodel(p_best,stateswitch)
%close all;
% fix alpha, to avoid degeneracy


if stateswitch==0
alpha=p_best(1);
end

if stateswitch==1
alpha=1;
end

if stateswitch==2
alpha=0;
end
mur1=p_best(2);

sigmar1=p_best(3);

sigmatheta1=p_best(4);

mur2=p_best(5);

sigmar2=p_best(6);

sigmatheta2=p_best(7);







t=1;
n=t;
run=100000;
dt=1;
coord=zeros(2,n);
state=zeros(1,n);
diff=zeros(2,n);
A=zeros(2,run,n);

thetazero=0;
    
for m=1:run
   
    thetazero2=thetazero;%random('unif',0,b);
 
for i=1:n
 
    if m/run<alpha
        mur=mur1;
        sigmar=sigmar1;
        sigmatheta=sigmatheta1;
        state(m)=1;
    else if m/run>=alpha
        mur=mur2;
        sigmar=sigmar2;
        sigmatheta=sigmatheta2;
        end
        state(m)=2;
    end
    
Vr=random('Normal',mur,sigmar);

Vtheta=thetazero2+random('Normal',0,sigmatheta);
thetazero2=Vtheta;

theta=Vtheta;    

if i~=1
coord(:,i)=[coord(1,i-1)+dt*Vr*cos(theta),coord(2,i-1)+dt*Vr*sin(theta)];
diff(1,i)=coord(1,i)-coord(1,i-1);
diff(2,i)=coord(2,i)-coord(2,i-1);
else
coord(:,i)=[dt*Vr*cos(theta),dt*Vr*sin(theta)];
diff(1,i)=coord(1,i);
diff(2,i)=coord(2,i);
end

%diff(i*m,1)=coord(1,i);%-coord(1,i-1,m);
%diff(i*m,2)=coord(2,i);%-coord(2,i-1,m);
%thetastore(i)=theta;

end
%theta=random('unif',0,b);
a=(m-1)*n+1;
b=m*n;;
size(A);
%A(:,a:b)=diff;
A(:,m,:)=diff;

%plot(coord(1,:),coord(2,:));figure(gcf);
%A(1,m)=Vtheta;
add(:,2:3)=coord';
add(:,1)=1:1:n;
add(2:n,4:5)=coord(:,2:n)'-coord(:,1:n-1)';
add(:,6)=m;
%add(:,7)=thetastore;

if m==1
out=add;
else
out=[out;add];
end


end

   controldata=out; 
    figure(2);

[n,xout]=scalehist(A(1,:,t),60);
 %n=n/(sum(n)*(xout(2)-xout(1)));
 %hold on;
 %bar(xout,n);figure(gcf);

 cellstepstotal=A';
 cellstepstotal(:,3)=1;
 cellstepstotal(:,4)=1:run;
end