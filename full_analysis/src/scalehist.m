function [n,xout] = scalehist(A,num,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if isempty(varargin)
[n,xout]=hist(A,num);
 n=n/(sum(n)*(xout(2)-xout(1)));
 bar(xout,n);figure(gcf);
else

    if numel(varargin)==1
    
    Bounds=varargin{1};
 Bins=linspace(Bounds(1),Bounds(2),num);
    if Bounds(1)<0
   
    
    [~,shiftI]=min(abs(Bins));
    
    %make sure there is a bin centered on 0
    Bins=Bins-Bins(shiftI);
    
    end
    figure
    [n,xout]=hist(A,Bins);
 n=n/(sum(n)*(xout(2)-xout(1)));
 bar(xout,n);figure(gcf);
        
    end
    
end

end

