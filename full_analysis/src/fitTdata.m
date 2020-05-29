function [mur,sigmar,sigmatheta,mur2,sigmar2,sigmatheta2,alpha,alphamin,score,murmin,sigmarmin,sigmathetamin,murmin2,sigmarmin2,sigmathetamin2,scoremin,cellnewx,out,xout,yout] = fitTdata(cellsteps,n,search,scale)

% This function performs a custom simulated annealing routine to find best
% estimates for the model parameters.

% cellsteps is the transformed data from prepare_fit_tracks
% n is the number of bins to sort the data into, select such that histogram looks
% reasonable.
% search is the number of iterations of the SA
% scale increases the threshold for throwing data out. (helpful to throw
% out large outliers)


%% Preliminaries and Remove Outliers
% initial guess 
higher=[20 20 3.5   20 20 3.5   1];
lower= [0 1 0.1   0  1  0.1 0.01];

FittingParams = lower + (higher-lower).*rand(1,7);
Fittingcurrent=FittingParams;

mur=zeros(search,1);
score=zeros(search,1);
mur2=mur;
sigmar=mur;
sigmatheta=mur;
alpha=mur;
sigmar2=mur;
sigmatheta2=mur;

mur(1)=FittingParams(1);
sigmar(1)=FittingParams(2);
thetazero=0;
sigmatheta(1)=FittingParams(3);

mur2(1)=FittingParams(4);
sigmar2(1)=FittingParams(5);
thetazero2=0;
sigmatheta2(1)=FittingParams(6);
alpha(1)=FittingParams(7);

i=1;

murcurrent=mur(1);
sigmarcurrent=sigmar(1);
sigmathetacurrent=sigmatheta(1);

murcurrent2=mur2(1);
sigmarcurrent2=sigmar2(1);
sigmathetacurrent2=sigmatheta2(1);


alphacurrent=alpha(1);
count=1;

for k=1:max(size(cellsteps))

    if(abs(cellsteps(k,1))<20*scale)
        if (abs(cellsteps(k,2))<20*scale) %10, %20, or 30 used
        cellnewx(count,:)=cellsteps(k,:);
        count=count+1;
        end
    end

end

%% score the initial guess
p(1)=alpha(1);
p(2)=mur(1);
p(3)=sigmar(1);
p(4)=sigmatheta(1);
p(5)=mur2(1);
p(6)=sigmar2(1);
p(7)=sigmatheta2(1);
int_acc=1000;
[nx,xout]=scalehist(cellnewx(:,1),n,[min(cellnewx(:,1)),max(cellnewx(:,1))]);
[ny,yout]=scalehist(cellnewx(:,2),n,[min(cellnewx(:,2)),max(cellnewx(:,2))]);
[scorecurrent,xout,~]=score_SA(p,xout,nx,yout,ny,int_acc,0);

score(1)=scorecurrent;

i=1;

% inital noise to guess new parameters
sd=0.15;
T=1;

%% perform the SA
for i=2:search
    if mod(i,100)==0
    i
    end
    if T<3750
    T=0.996*T;
    end

       
     if i>search*4/5
     int_acc=500;
     else
     int_acc=50;
     end
        
        
  mur(i)=  abs(murcurrent+random('norm',0,sd));
  sigmar(i)=abs(sigmarcurrent+random('norm',0,sd));
  sigmatheta(i)=abs(sigmathetacurrent +random('norm',0,sd));
  
mur2(i)=  murcurrent2+random('norm',0,sd);
sigmar2(i)=abs(sigmarcurrent2+random('norm',0,sd));
sigmatheta2(i)=abs(sigmathetacurrent2 +random('norm',0,sd));
 alpha(i)=abs(alphacurrent+random('norm',0,sd/5));


if mur(i)>higher(1)
    mur(i)=higher(1);
end

if mur(i)<lower(1)
    mur(i)=lower(1);
end


if sigmar(i)>higher(2)
    sigmar(i)=higher(2);
end

if  sigmar(i)<lower(2)
    sigmar(i)=lower(2);
end

if sigmatheta(i)>higher(3)
    sigmatheta(i)=higher(3);
end

if  sigmatheta(i)<lower(3)
    sigmatheta(i)=lower(3);
end


if mur2(i)>higher(4)
    mur2(i)=higher(4);
end

if mur2(i)<lower(4)
    mur2(i)=lower(4);
end


if sigmar2(i)>higher(5)
    sigmar2(i)=higher(5);
end

if  sigmar2(i)<lower(5)
    sigmar2(i)=lower(5);
end

if sigmatheta2(i)>higher(6)
    sigmatheta2(i)=higher(6);
end

if  sigmatheta2(i)<lower(6)
    sigmatheta2(i)=lower(6);
end


if alpha(i)>higher(7)
    alpha(i)=higher(7);
end

if  alpha(i)<lower(7)
    alpha(i)=lower(7);
end

  

  if sigmatheta(i)>3.5
      sigmatheta(i)=3.5;
  end
  
  if sigmatheta2(i)>3.5
      sigmatheta2(i)=3.5;
  end
  
  
  if sigmatheta(i)<0.2
      sigmatheta(i)=0.2+sigmatheta(i);
  end
  
  if sigmatheta2(i)<0.2
      sigmatheta2(i)=0.2+sigmatheta(i);
  end
  

  if sigmar(i)<0.2
      sigmar(i)=0.2+sigmar(i);
  end
  
  if sigmar2(i)<0.2
      sigmar2(i)=0.2+sigmar2(i);
  end
  
 
p(1)=alpha(i);
p(2)=mur(i);
p(3)=sigmar(i);
p(4)=sigmatheta(i);
p(5)=mur2(i);
p(6)=sigmar2(i);
p(7)=sigmatheta2(i);

[score(i),~,~]=score_SA(p,xout,nx,yout,ny,int_acc,0);

if score(i)<scorecurrent
   
   murcurrent=mur(i);
   sigmarcurrent=sigmar(i);
   sigmathetacurrent=sigmatheta(i);
   scorecurrent=score(i);
   out(i)=scorecurrent;

   murcurrent2=mur2(i);
   sigmarcurrent2=sigmar2(i);
   sigmathetacurrent2=sigmatheta2(i);
   alphacurrent=alpha(i);
   
  Fittingcurrent=[murcurrent,sigmarcurrent,sigmathetacurrent,murcurrent2,sigmarcurrent2,sigmathetacurrent2,alphacurrent];
else
   test=exp((scorecurrent-score(i))/(0.1*T));
    
   r=unifrnd(0,1);
   if r<test
          murcurrent=mur(i);
   sigmarcurrent=sigmar(i);
   sigmathetacurrent=sigmatheta(i);
   scorecurrent=score(i);
   out(i)=scorecurrent;
   
   
   murcurrent2=mur2(i);
   sigmarcurrent2=sigmar2(i);
   sigmathetacurrent2=sigmatheta2(i);
   alphacurrent=alpha(i);
   
   end    
end
 out(i)=scorecurrent;

end


scoremin=min(score);
minvalue=find(score==scoremin);
sigmarmin=sigmar(minvalue(1));
sigmathetamin=sigmatheta(minvalue(1));
murmin=mur(minvalue(1));
sigmarmin2=sigmar2(minvalue(1));
sigmathetamin2=sigmatheta2(minvalue(1));
murmin2=mur2(minvalue(1));
alphamin=alpha(minvalue(1));


end

