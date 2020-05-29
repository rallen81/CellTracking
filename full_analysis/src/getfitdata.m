function p_best=getfitdata(loadstr2,definedsplit,variablesplit,unit_conv)

%this function loads the data saved from prepare_fit_tracks and identifies
%and plots the best fit. It also can be used to look at this histograms of all fits,
%or the best fits. 


%loadstr:
    % file name of the data saved from prepare_fit_tracks
% definedsplit:
    % Values less than definiedsplit will be state 1, otherwise
    % state 2
%variablesplit:
    % switch to define which parameter to split the states on:
        % 1 = mur
        % 2 = sigmar
        % 3 = sigmatheta
%unit_conv:
    % conversion factor to change units to microns for plotting
    
load(loadstr2)

%% population fits
murcombined=[murstore; murstore2];
sigmathetacombined=[sigmathetastore; sigmathetastore2];
sigmarcombined=[sigmarstore; sigmarstore2];
alphacombined=[alphastore; (1-alphastore)];
scorecombined=[scorestore; scorestore];

[n,yout]=scalehist(fitteddata(:,2),numb,[min(fitteddata(:,2)),max(fitteddata(:,2))]);
[n,xout]=scalehist(fitteddata(:,1),numb,[min(fitteddata(:,1)),max(fitteddata(:,1))]);
xFitted=xout;
yFitted=yout;
A=zeros(2,5);
B=zeros(2,5);
C=zeros(2,5);
D=zeros(2,5);

[CI_5,CI_95,SSE,SSE2] = bootstrapsse(fitteddata,numb);
top5=prctile(SSE,99);



bestfits=scorecombined<=top5;
murbest=murcombined(bestfits);
scorecombined=scorecombined(bestfits);

sigmarbest=sigmarcombined(bestfits);
sigmathetabest=sigmathetacombined(bestfits);
alphabest=alphacombined(bestfits);




if variablesplit==1
hist(murbest,500);figure(gcf); %velocity
state1=murbest<definedsplit & murbest>0;
state2=murbest>=definedsplit;
end

if variablesplit==2
hist(sigmarbest,500);figure(gcf); %variation in velocity
state1=sigmarbest<definedsplit ;
state2=sigmarbest>=definedsplit;
end

if variablesplit==3
hist(sigmathetabest,500);figure(gcf); %persistence 
state1=sigmathetabest<definedsplit ;
state2=sigmathetabest>=definedsplit;
end
%%%%%%%%%%%

mur1=murbest(state1);
mur2=murbest(state2);
sigmar1=sigmarbest(state1);
sigmar2=sigmarbest(state2);
sigmatheta1=sigmathetabest(state1);
sigmatheta2=sigmathetabest(state2);


medmur1=median(murbest(state1));
meanmur1=mean(murbest(state1));

medmur2=median(murbest(state2));
meanmur2=mean(murbest(state2));

stdmur1=std(murbest(state1));
stdmur2=std(murbest(state2));

fifththmur1=prctile(murbest(state1),5);
ninefifthmur1=prctile(murbest(state1),95);

fifththmur2=prctile(murbest(state2),5);
ninefifthmur2=prctile(murbest(state2),95);


A(1,1)=medmur1;
A(1,2)=meanmur1;
A(1,3)=stdmur1;
A(1,4)=fifththmur1;
A(1,5)=ninefifthmur1;

A(2,1)=medmur2;
A(2,2)=meanmur2;
A(2,3)=stdmur2;
A(2,4)=fifththmur2;
A(2,5)=ninefifthmur2;
A=A';

medsigmar1=median(sigmarbest(state1));
meansigmar1=mean(sigmarbest(state1));

medsigmar2=median(sigmarbest(state2));
meansigmar2=mean(sigmarbest(state2));

stdsigmar1=std(sigmarbest(state1));
stdsigmar2=std(sigmarbest(state2));

fifththsigmar1=prctile(sigmarbest(state1),5);
ninefifthsigmar1=prctile(sigmarbest(state1),95);

fifththsigmar2=prctile(sigmarbest(state2),5);
ninefifthsigmar2=prctile(sigmarbest(state2),95);

B(1,1)=medsigmar1;
B(1,2)=meansigmar1;
B(1,3)=stdsigmar1;
B(1,4)=fifththsigmar1;
B(1,5)=ninefifthsigmar1;


B(2,1)=medsigmar2;
B(2,2)=meansigmar2;
B(2,3)=stdsigmar2;
B(2,4)=fifththsigmar2;
B(2,5)=ninefifthsigmar2;
B=B';
medsigmatheta1=median(sigmathetabest(state1));
meansigmatheta1=mean(sigmathetabest(state1));

medsigmatheta2=median(sigmathetabest(state2));
meansigmatheta2=mean(sigmathetabest(state2));

stdsigmatheta1=std(sigmathetabest(state1));
stdsigmatheta2=std(sigmathetabest(state2));

fifththsigmatheta1=prctile(sigmathetabest(state1),5);
ninefifthsigmatheta1=prctile(sigmathetabest(state1),95);

fifththsigmatheta2=prctile(sigmathetabest(state2),5);
ninefifthsigmatheta2=prctile(sigmathetabest(state2),95);

C(1,1)=medsigmatheta1;
C(1,2)=meansigmatheta1;
C(1,3)=stdsigmatheta1;
C(1,4)=fifththsigmatheta1;
C(1,5)=ninefifthsigmatheta1;


C(2,1)=medsigmatheta2;
C(2,2)=meansigmatheta2;
C(2,3)=stdsigmatheta2;
C(2,4)=fifththsigmatheta2;
C(2,5)=ninefifthsigmatheta2;
C=C';




medalpha=median(alphabest(state1));
meanalpha=mean(alphabest(state1));
stdalpha=std(alphabest(state1));

fifthtalpha=prctile(alphabest(state1),5);
ninefifthalpha=prctile(alphabest(state1),95);



D(1,1)=medalpha;
D(1,2)=meanalpha;
D(1,3)=stdalpha;
D(1,4)=fifthtalpha;
D(1,5)=ninefifthalpha;
D=D';

%% best fit
figure(10)
[n,xout]=scalehist(fitteddata(:,1),xFitted);


xmin=min(xout);
xmax=max(xout);
xplot=xmin:0.05:xmax;

minscore=min(scorestore);
a=find(scorestore==minscore);
a;
a=a(1);



bestsigmar=sigmarstore(a);
bestmur=murstore(a);
bestsigmatheta=sigmathetastore(a);

bestsigmar2=sigmarstore2(a);
bestmur2=murstore2(a);
bestsigmatheta2=sigmathetastore2(a);

bestalpha=alphastore(a);

p_best(1)=bestalpha;
p_best(2)=bestmur;
p_best(3)=bestsigmar;
p_best(4)=bestsigmatheta;
p_best(5)=bestmur2;
p_best(6)=bestsigmar2;
p_best(7)=bestsigmatheta2;

%% finalize the best fit with local optimization
[p_best]=profile_likelihood(p_best,cellstepstotal,80,fitteddata,numb,0,0);

%% plot the best fit
[Xf,u]=mypdfmaster(p_best(3),p_best(2),p_best(4),0,xplot);
[Xf2,u2]=mypdfmaster(p_best(6),p_best(5),p_best(7),0,xplot);
Xfinal=Xf*p_best(1)+Xf2*(1-p_best(1));

figure
[n,xout]=scalehist(fitteddata(:,1),xFitted);
hold off

figuremain=mainfitfigure(xout*unit_conv,n,xplot*unit_conv,Xfinal);
figurein_1=insetfitfigure2(xplot*unit_conv,Xf,1);
figurein_2=insetfitfigure2(xplot*unit_conv,Xf2,2);
[h1,~,inseth]=inset2(figuremain, figurein_1,figurein_2);

pos = [1000         828         625         503];
set(gcf,'position',pos);
cwd=pwd;


%% save the figures and data
cd figures
strcat(loadstr2,'_DeltaXfit_2insets');
savefig(gcf,strcat(loadstr2(1,9:end-4),'_DeltaXfit_2insets'));
cd(cwd);

[Yf,u]=mypdfmaster(p_best(3),p_best(2),p_best(4),pi/2,xplot);
[Yf2,u2]=mypdfmaster(p_best(6),p_best(5),p_best(7),pi/2,xplot);
Yfinal=Yf*p_best(1)+Yf2*(1-p_best(1));




%hold off
figure
[n,xout]=scalehist(fitteddata(:,2),xFitted);
%hold on;

figuremain2=mainfitfigure(xout*unit_conv,n,xplot*unit_conv,Yfinal);
xlabel('Y (\mum)')

%figurein2=insetfitfigure(xplot*unit_conv,p_best(1)*Yf,xplot*unit_conv,(1-p_best(1))*Yf2);
%[h2,~,inseth2]=inset(figuremain2, figurein2);

figurein_12=insetfitfigure2(xplot*unit_conv,Yf,1);
figurein_22=insetfitfigure2(xplot*unit_conv,Yf2,2);
[h2,~,inseth2]=inset2(figuremain2, figurein_12,figurein_22);


pos = [1000         828         625         503];
set(gcf,'position',pos);
hold off;


cwd=pwd;
cd figures
strcat(loadstr2,'_DeltaYfit_2insets');
savefig(gcf,strcat(loadstr2(1,9:end-4),'_DeltaYfit_2insets'));
cd(cwd);


save(loadstr2,'murstore','murstore2','sigmathetastore','sigmathetastore2','sigmarstore','sigmarstore2','alphastore','scorestore','fitteddata','numb','p_best','cellstepstotal')
 %savestr=strcat('analysis_',loadstr2);
% save(savestr)
close all
end

