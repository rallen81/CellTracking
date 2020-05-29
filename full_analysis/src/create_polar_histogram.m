function [h,dtheta] = create_polar_histogram(dataset,savestr)

load(dataset)
controldata(:,1)=round(controldata(:,1));
controldata(:,6)=round(controldata(:,6));

figure
hold on
I=abs(controldata(:,4))<50;
%histogram(controldata(I,4),50)

scale_controldata=(controldata(I,4)-mean(controldata(I,4)))/std(controldata(I,4));

[hx,px]=kstest(scale_controldata)

[f,x_values]=ecdf(scale_controldata);
F=plot(x_values,f);
set(F,'LineWidth',2);

G=plot(-10:0.1:10,normcdf(-10:0.1:10,0,1),'r-');
set(G,'LineWidth',2);


I=abs(controldata(:,5))<50;
%histogram(controldata(I,5),50)

scale_controldata=(controldata(I,5)-mean(controldata(I,5)))/std(controldata(I,5));

[hy,py]=kstest(scale_controldata)

[f,x_values]=ecdf(scale_controldata);
F2=plot(x_values,f);
set(F2,'LineWidth',2);

legend([F F2 G],'Scaled Empirical CDF for x','Scaled Empirical CDF for y', 'Standard Normal CDF', 'Location', 'SE')
xlabel('Distance (Normalized)');
ylabel('CDF')
xlim([-10,10])
set(gca,'FontSize',20)
% calculate dr
S=size(controldata);
N=S(1);
k=1;
i=2;
skiplist=[];
threshold=0.1^2;
while i<=N

if controldata(i,1)==1 || controldata(i,1)<controldata(i-1,1)
i=i+1;
kstart=k;
end

%%
if i<N
   
if (controldata(i,1)==controldata(i-1,1)+1)
    
if (controldata(i,4).^2+controldata(i,5).^2).^0.5>threshold
dr(k,1)=controldata(i,6);
dr(k,2)=controldata(i,1);
dr(k,3)=(controldata(i,4).^2+controldata(i,5).^2).^0.5;

   
x1(k,1)=controldata(i,6);
x1(k,2)=controldata(i,1);
x1(k,3)=controldata(i,2);

y1(k,1)=controldata(i,6);
y1(k,2)=controldata(i,1);
y1(k,3)=controldata(i,3);


theta(k,1)=controldata(i,6);
theta(k,2)=controldata(i,1);
theta(k,3)=atan2(controldata(i,5),controldata(i,4));
calctheta(k)=1;


else
dr(k,1)=controldata(i,6);
dr(k,2)=controldata(i,1);
dr(k,3)=(controldata(i,4).^2+controldata(i,5).^2).^0.5;


theta(k,1)=controldata(i,6);
theta(k,2)=controldata(i,1);   

x1(k,1)=controldata(i,6);
x1(k,2)=controldata(i,1);
x1(k,3)=controldata(i,2);

y1(k,1)=controldata(i,6);
y1(k,2)=controldata(i,1);
y1(k,3)=controldata(i,3);


if k>1
    if calctheta(k-1)==1
        theta(k,3)=theta(k-1,3);
        calctheta(k)=-1;
    elseif k>2
            if calctheta(k-2)==1
                theta(k,3)=theta(k-2,3);
                calctheta(k)=-2;
            else
                theta(k,3)=2*pi*rand;
                calctheta(k)=0;
            end
    else
     theta(k,3)=2*pi*rand;
     calctheta(k)=0;
    end
    
else
     theta(k,3)=2*pi*rand;
     calctheta(k)=0;
end
end

k=k+1;
i=i+1;

elseif controldata(i,1)>(controldata(i-1,1)+1)
  
m=i;

    while controldata(m,6)==controldata(i,6)
        m=m+1;
    end
 i=m; 
 k=kstart;

end
end
if i==N
i=i+1;
end
%
end

S2=size(dr);
n=max(S2);
count1=1;
count2=2;

%makepolarcoordinates
k=1;
k1=1;


Ntheta=max(size(theta));


for j=2:Ntheta
    
    
    if theta(j,2)>theta(j-1,2)
    dtheta(k)=theta(j,3)-theta(j-1,3);


while dtheta(k)>pi
   dtheta(k)=dtheta(k)-2*pi ;
end

while dtheta(k)<-pi
   dtheta(k)=dtheta(k)+2*pi ;
end

 k=k+1;
    end
end


newN=max(size(x1));
%A=rotate(x1KD(2:newN,3)-x1KD(1:newN-1,3),y1KD(2:newN,3)-y1KD(1:newN-1,3),-thetaKD(2:(newN),3));
%B=rotate(x1KD(3:newN,3)-x1KD(1:newN-2,3),y1KD(3:newN,3)-y1KD(1:newN-2,3),-thetaKD(2:(newN-1),3));
%[diffx1(:,3:4)]=B-A(1:newN-2,:);

%diffx1(1:(newN-1),1)=y1KD(2:newN,1);
%diffx1(1:(newN-1),2)=y1KD(2:newN,2);

cellnumber=max(controldata(:,6));
i=1;
for i=1:newN
    
  
    x1new(x1(i,2)-1,x1(i,1),1)=x1(i,3);
    y1new(y1(i,2)-1,y1(i,1),1)=y1(i,3);
    thetanew(theta(i,2)-1,theta(i,1))=theta(i,3);
    i=i+1;
  
end


%REMOVE ZERO COLUMNS
%x1new( :, ~all(x1new,1) ) = [];  %columns
%y1new( :, ~all(y1new,1) ) = [];  %columns

nCell2d=size(x1new);
nCell=nCell2d(1);


%% change to size of x1new
for i=1:nCell2d(2)
  
    Anew=rotate(x1new(2:nCell,i)-x1new(1:nCell-1,i),y1new(2:nCell,i)-y1new(1:nCell-1,i),-thetanew(2:(nCell),i));
    Bnew=rotate(x1new(3:nCell,i)-x1new(1:nCell-2,i),y1new(3:nCell,i)-y1new(1:nCell-2,i),-(thetanew(2:(nCell-1),i)));
   
    cellsteps(:,1:2)=Bnew-Anew(1:nCell-2,:);
    cellsteps(:,3)=i;
     cellsteps(:,4)=3:nCell;
     
    
   % [mur,sigmar,sigmatheta,score,murcurrent,sigmarcurrent,sigmathetacurrent,scorecurrent] = fitdata(cellsteps,1);
   % scorestore(i)=scorecurrent;
   % mustore(i)=murcurrent;
   % sigmacurrentstore(i)=sigmarcurrent;
   % sigmathetacurrentstore(i)=sigmathetacurrent;
   if i~=1
    cellstepstotal=[cellstepstotal; cellsteps];
   else
       cellstepstotal=cellsteps;
   end
   
   
end


%numb=35;%30;

%numoffits=8; %80;




figure
h=polarhistogram(dtheta,'Normalization','probability');
ax=gca;
ax.ThetaAxisUnits='radians';
ax.FontSize=24;
ax.LineWidth=2;
thetaticks(0:pi/4:2*pi-pi/4);


cd figures
savefig(savestr)
print(savestr,'-deps')

cd ..
end

