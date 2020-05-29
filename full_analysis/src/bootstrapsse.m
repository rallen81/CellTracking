function [ CI_5,CI_95,SSE,SSE2,nx,xout] = bootstrapsse(data,numb)
%data is an n x 4 matrix,
% column 1 is the data we are interested in
% column 2 is additional data we might be interested in
% column 3 is an observation set (here a cell track) 
% column 4 is a time index (not used)

[nx,xout]=scalehist(data(:,1),numb,[min(data(:,1)),max(data(:,1))]);
%[score,xout,Xf]=score_SA(p,xout,nx,1000);
for i=1:4000
    
y1=1:max(data(:,3));
y2=randsample(max(data(:,3)),max(data(:,3)),1); %nchoosen with replacement

data1=resampled_data(y1); %original data
data2=resampled_data(y2); %random sample

data3=randsample(data(:,1),numel(data(:,1)),1);
[nx2,xout]=hist(data3,xout);
nx2=nx2./(sum(nx2)*(xout(2)-xout(1)));
SSE(i)=sum((data1-data2).^2);
SSE2(i)=sum((data1-nx2).^2);
end


CI_5=prctile(SSE,5);
CI_95=prctile(SSE,95);

    function data2=resampled_data(y)
        for j=1:numel(y)
        dataI=data(:,3)==y(j);
        [ntemp(j,:),xouttemp]=hist(data(dataI,1),xout);
        end
        data2=sum(ntemp);
        
        data2=data2./(sum(data2)*(xout(2)-xout(1)));
    end

    
end

