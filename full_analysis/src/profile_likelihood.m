function [p_best_opt]=profile_likelihood(p_best,cellstepstotal,num_of_checks,fitteddata,numb,flag_profile,yflag,varargin)

% This function sweeps through parameter space by fixing parameters at a
% certain value and then optimizing to find the best fit around the other
% parameters. See paper supplementary materials

% p_best, best estimate of the global minimum.
% fitteddata, data we are fitting
% numb - number of bins,
% switch for profile, 1 perform profile, 0 just return refined best fit.
% use y-data also
%optional arguments
% savestr save name for figures/data
% scale - scale factor for deciding range to search for parameters in

close all

if isempty(varargin)==0
scale=varargin{2};
else
scale=1;
end

[nx,xout]=scalehist(fitteddata(:,1),numb,[min(fitteddata(:,1)),max(fitteddata(:,1))]);
[ny,yout]=scalehist(fitteddata(:,2),numb,[min(fitteddata(:,2)),max(fitteddata(:,2))]);
%p_best(1)=bestalpha;

%% ranges to search in

% fix alpha, to avoid degeneracy

upper_p(1)=p_best(1)*1.0001; %
lower_p(1)=p_best(1)*0.9999;

%p_best(2)=bestmur;
upper_p(2)=max(xout)*scale;
lower_p(2)=0;

%p_best(3)=bestsigmar;
upper_p(3)=max(xout)*scale;
lower_p(3)=0.01; %measurement accuracy (essentially)

%p_best(4)=bestsigmatheta;
upper_p(4)=3.5;
lower_p(4)=0.01;

%p_best(5)=bestmur2;
upper_p(5)=max(xout)*scale;
lower_p(5)=0;

%p_best(6)=bestsigmar2;
upper_p(6)=max(xout)*scale;
lower_p(6)=0.01; %measurement accuracy (essentially)

%p_best(7)=bestsigmatheta2;
upper_p(7)=3.5;
lower_p(7)=0.01;


%% re-optimize to find truly best
options = optimoptions('fmincon');
options.Display='iter';
options.UseParallel=1;

best=score_SA(p_best,xout,nx,yout,ny,1000,yflag)
p_store=zeros(7,num_of_checks+1,7);
val_store=zeros(7,num_of_checks+1);
p_singleton=zeros(7,num_of_checks+1);

f=@(p)score_SA(p,xout,nx,yout,ny,1000,yflag);
p_best
lower_p;
upper_p;
[pout,fval]=fmincon(f,p_best,[],[],[],[],lower_p,upper_p,[],options);
fval
p_best_opt=pout;

%% run through each parameter fixing at every point in range and re_optimizing.
if flag_profile==1
parfor k=1:7
    upper_p_temp=upper_p;
    lower_p_temp=lower_p;
    
    p_temp=p_best_opt;
    
   
    for j=1:num_of_checks+1
        p_temp(k)=lower_p(k)+(upper_p(k)-lower_p(k))*(j-1)/num_of_checks;
        if k==1
            p_temp(k)=(j-1)/num_of_checks;
            upper_p_temp=p_temp*1.001;
            lower_p_temp=p_temp*0.999;
        end
        j
        f=@(p)score_SA(p,xout,nx,yout,ny,1000,yflag);
        upper_p_temp(k)=p_temp(k)*1.005;
        lower_p_temp(k)=p_temp(k)*0.995;
        p_temp
        
        [pout,fval]=fmincon(f,p_temp,[],[],[],[],lower_p_temp,upper_p_temp,[],options);
        
       fval
        p_store(k,j,:)=pout;
        val_store(k,j)=fval;
        p_singleton(k,j)=p_temp(k);
    end
    

end
end
%remove-outliers

for k=1:7
for j=2:num_of_checks
if val_store(k,j)>1.2*val_store(k,j-1) && val_store(k,j)>1.2*val_store(k,j+1)
val_store(k,j)=0.5*(val_store(k,j-1)+val_store(k,j+1));
end
end




end

%% predict the states based on data
[plotstate,state]=predictstates_fun(p_best,cellstepstotal);

%% save results
if isempty(varargin)==0
cd results
savestr=varargin{1};
save(savestr);
cd ..
end
end