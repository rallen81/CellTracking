function prepare_fit_tracks(loadstr,savestr,numb,numoffits)

% this function firstly prepares the data for fitting by transforming triplets in
% the X-Y coordinate data by rigid body rotation such that Y2-Y1 = 0, i.e.
% perpendicular to x-axis for the first pair of coordinates. 

% loadstr:
%string of data load, manipulate and fit  must contain n x 6 
% variable called controldata with:
%     1st column step number (1 to m integer)
%     2nd & 3rd column X-Y coordinates
%     4th & 5th column Delta X and Y (not used, I believe)
%     6th column, cell number (1 to k integer)

% savestr: 
% string determining filename to save outputs too.

% numb
% number of bins to sort the data into, select such that histogram looks
% reasonable.

% numberoffits
% how many simulated annealing runs to perform



stream0 = RandStream('mt19937ar','Seed',sum(100*clock));
RandStream.setGlobalStream(stream0);


load(loadstr)


%% parse the data
controldata(:,1)=round(controldata(:,1));
controldata(:,6)=round(controldata(:,6));
 

% calculate dr
S=size(controldata);
N=S(1);
k=1;
i=2;

threshold=0.1^2;

%% go through the data line by line, store calculations and seperate out by cell
while i<=N

if controldata(i,1)==1 || controldata(i,1)<controldata(i-1,1)
i=i+1;
kstart=k;
end


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



k=1;
Ntheta=max(size(theta));

%% calculate dtheta
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

cellnumber=max(controldata(:,6));
i=1;
for i=1:newN
    
  
    x1new(x1(i,2)-1,x1(i,1),1)=x1(i,3);
    y1new(y1(i,2)-1,y1(i,1),1)=y1(i,3);
    thetanew(theta(i,2)-1,theta(i,1))=theta(i,3);
    i=i+1;
  
end


nCell2d=size(x1new);
nCell=nCell2d(1);

%% Do rotation
for i=1:nCell2d(2)
  
    Anew=rotate(x1new(2:nCell,i)-x1new(1:nCell-1,i),y1new(2:nCell,i)-y1new(1:nCell-1,i),-thetanew(2:(nCell),i));
    Bnew=rotate(x1new(3:nCell,i)-x1new(1:nCell-2,i),y1new(3:nCell,i)-y1new(1:nCell-2,i),-(thetanew(2:(nCell-1),i)));
   
    cellsteps(:,1:2)=Bnew-Anew(1:nCell-2,:);
    cellsteps(:,3)=i;
     cellsteps(:,4)=3:nCell;
     
   
   if i~=1
    cellstepstotal=[cellstepstotal; cellsteps];
   else
       cellstepstotal=cellsteps;
   end
   
   
end

%% do fitting, store last 500 iterations
alphastore=zeros(500,numoffits); 
murstore=zeros(500,numoffits);
sigmarstore=zeros(500,numoffits);
sigmathetastore=zeros(500,numoffits);

murstore2=zeros(500,numoffits);
sigmarstore2=zeros(500,numoffits);
sigmathetastore2=zeros(500,numoffits);
scorestore=zeros(500,numoffits);


parfor j=1:numoffits
    j

[mur,sigmar,sigmatheta,mur2,sigmar2,sigmatheta2,alpha,~,score,~,~,~,~,~,~,~,fitteddata,out2,xout,yout] = fitTdata(cellstepstotal,numb,4000,1);

maxA=numel(mur);

alphastore(:,j)=alpha((maxA-500+1):maxA);
murstore(:,j)=mur((maxA-500+1):maxA);
sigmarstore(:,j)=sigmar((maxA-500+1):maxA);
sigmathetastore(:,j)=sigmatheta((maxA-500+1):maxA);

murstore2(:,j)=mur2((maxA-500+1):maxA);
sigmarstore2(:,j)=sigmar2((maxA-500+1):maxA);
sigmathetastore2(:,j)=sigmatheta2((maxA-500+1):maxA);  
scorestore(:,j)=score((maxA-500+1):maxA);  


end
[mur,sigmar,sigmatheta,mur2,sigmar2,sigmatheta2,alpha,~,score,~,~,~,~,~,~,~,fitteddata,out2] = fitTdata(cellstepstotal,numb,10,1);
alphastore=alphastore(:);
murstore=murstore(:);
sigmarstore=sigmarstore(:);
sigmathetastore=sigmathetastore(:);
murstore2=murstore2(:);
sigmarstore2=sigmarstore2(:);
sigmathetastore2=sigmathetastore2(:);
scorestore=scorestore(:);

close all

%% save results
cd results
save(savestr,'murstore','murstore2','sigmathetastore','sigmathetastore2','sigmarstore','sigmarstore2','alphastore','scorestore','fitteddata','numb','cellstepstotal')
cd ..
