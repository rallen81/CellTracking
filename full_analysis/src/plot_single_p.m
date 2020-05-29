function [figuremain,X1,Y1] = plot_single_p(p_best_1,p_best,numb,fitteddata,unit_conv,state_select)

[n,xout]=scalehist(fitteddata(:,1),numb,[min(fitteddata(:,1)),max(fitteddata(:,1))]);
xFitted=xout;
xmin=min(xout);
xmax=max(xout);
xplot=xmin:0.05:xmax;

figure
[n,xout]=scalehist(fitteddata(:,1),xFitted);

if state_select == 1
[Xf,u]=mypdfmaster(p_best_1(3),p_best_1(2),p_best_1(4),0,xplot);
Xfinal=Xf;
end

if state_select == 2
[Xf2,u2]=mypdfmaster(p_best_1(6),p_best_1(5),p_best_1(7),0,xplot);
Xfinal = Xf2;
end


[figuremain2,X1,Y1]=mainfitfigure(xout*unit_conv,n,xplot*unit_conv,Xfinal);
%figurein_1=insetfitfigure2(xplot*unit_conv,Xf,1);
%figurein_2=insetfitfigure2(xplot*unit_conv,Xf2,2);
%[h1,~,inseth]=inset2(figuremain, figurein_1,figurein_2);




%% plot the best fit
[Xf,u]=mypdfmaster(p_best(3),p_best(2),p_best(4),0,xplot);
[Xf2,u2]=mypdfmaster(p_best(6),p_best(5),p_best(7),0,xplot);
Xfinal=Xf*p_best(1)+Xf2*(1-p_best(1));

figure
[n,xout]=scalehist(fitteddata(:,1),xFitted);
hold off

figuremain=mainfitfigure(xout*unit_conv,n,xplot*unit_conv,Xfinal);
hold on

plot(X1,Y1,'--k','Linewidth',2)
figurein_1=insetfitfigure2(xplot*unit_conv,Xf,1);
figurein_2=insetfitfigure2(xplot*unit_conv,Xf2,2);
[h1,~,inseth]=inset2(figuremain, figurein_1,figurein_2);
insetfigure=h1;
pos = [1000         828         625         503];
set(gcf,'position',pos);
cwd=pwd;









end

