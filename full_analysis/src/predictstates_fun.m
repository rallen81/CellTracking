function [plotstate,state,state_forced]=predictstates_fun(p_best,cellstepstotal)


% predict the states based on bayes theorem, best parameter guess p_best,
% and the data (cellstepstotal).



alpha=p_best(1);
mur1=p_best(2);
sigmar1=p_best(3);
sigmatheta1=p_best(4);
mur2=p_best(5);
sigmar2=p_best(6);
sigmatheta2=p_best(7);


%%%%%%%%%%%%%%%%%%%%%%%
%%%process to make tracks all the same length

statesetsize=1; %should be odd
x=cellstepstotal(:,1);
M=max(size(x));

deltaX=0.01;
xin=-30:deltaX:30;

[Xf,x]=mysprocess(sigmar1,mur1,sigmatheta1,0,1,xin,1000);
[Xf2,x2]=mysprocess(sigmar2,mur2,sigmatheta2,0,1,xin,1000);

x=cellstepstotal(:,1);

index=1;
jlower=floor(statesetsize/2);
tracklength=max(cellstepstotal(:,4))-min(cellstepstotal(:,4))+1;
cellnumber=max(cellstepstotal(:,3));
ind=1;

for c=1:cellnumber
Logical_by_cell=cellstepstotal(:,4)==c;
TrackLength(c)=sum(Logical_by_cell);
min_track_length=min(TrackLength(c));


end

for cellstepstotal=1:cellnumber
tracklength;
for m=1+jlower:tracklength-jlower
i=(cellstepstotal-1)*tracklength+m;

%%%statesetsize=3
  %iout(c,m-jlower)=cellstepstotal(i,4)+(c-1)*(tracklength+2*jlower+1);
  
 %iout(c,m-jlower)=cellstepstotal(i,4)+(c-1)*(tracklength+3);
  
if index==0
    index=1;
end

for j=1:statesetsize
    
    

    
index=find(xin<(x(i-jlower+j-1)+deltaX/2) & xin>(x(i-jlower+j-1)-deltaX/2));
x(i-jlower+j-1);
  
if numel(index)==0

    
    [Xff,x22]=mysprocess(sigmar1,mur1,sigmatheta1,0,1,x(i-jlower+j-1),1000);
    [Xff2,x22]=mysprocess(sigmar2,mur2,sigmatheta2,0,1,x(i-jlower+j-1),1000);
    
    pdf1(j)=Xff;
    pdf2(j)=Xff2;
    
    
    
else

  pdf1(j)=Xf(index(1));
  pdf2(j)=Xf2(index(1));
end 
end

likelihoodstate1=prod(pdf1);
likelihoodstate2=prod(pdf2);

PState1givenx=likelihoodstate1*alpha/(likelihoodstate2*(1-alpha)+likelihoodstate1*alpha);
PState2givenx=likelihoodstate2*(1-alpha)/(likelihoodstate2*(1-alpha)+likelihoodstate1*alpha);


if PState1givenx>PState2givenx
  
    state(i)=1;
    pvalue(i)=PState2givenx;
   
else
  
    state(i)=2;
    pvalue(i)=PState1givenx;
end

    if (likelihoodstate2*(1-alpha)+likelihoodstate1*alpha)==0
    pvalue(i)=0;
    state(i)=2;
    end
    

   
   if state(i)==1
   plotstate1(i)=1-pvalue(i);
   plotstate2(i)=0;
   plotstate(cellstepstotal,m-jlower)=1+pvalue(i);
   
 
   end
   
   if state(i)==2
   plotstate2(i)=1-pvalue(i);
   plotstate1(i)=0;
   plotstate(cellstepstotal,m-jlower)=2-pvalue(i);
   
   end
   
   
  
end
 thresh=prctile(plotstate(:),alpha*100);
 
state_forced_idx=plotstate(:)<thresh;
state_forced=state;
state_forced(state_forced_idx)=1;
state_forced(~state_forced_idx)=2;
end
%cline(randn*0.15+controldata(iout(cell,1):iout(cell,tracklength-2*jlower),2)-controldata(iout(cell,1),2),randn*0.15+controldata(iout(cell,1):iout(cell,tracklength-2*jlower),3)-controldata(iout(cell,1),3),round(plotstate(cell,1:tracklength-2*jlower)))

    %X=controldata(iout(cell,1):iout(cell,tracklength-2*jlower),2);
    %Y=controldata(iout(cell,1):iout(cell,tracklength-2*jlower),3);
    
    %Xstate1=X(round(plotstate(cell,1:tracklength-2*jlower))==1);
    %Ystate1=Y(round(plotstate(cell,1:tracklength-2*jlower))==1);
    
    
    
    
    
    
  