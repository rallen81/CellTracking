function [cellstepstotal4,state1bycell4,state1forced4] = plot_state_predict(MEFcontrol, MEFKD,Helacontrol,HelaRac1,unit_conv1,unit_conv2)
%% This function is primarily for plotting tracks with the states predictions, arguments are the structures unpacked as below.

%% unpack structures
cellstepstotal1=MEFcontrol.cellstepstotal;
plotstate1=MEFcontrol.plotstate;
alpha1=MEFcontrol.p_best(1);

cellstepstotal2=MEFKD.cellstepstotal;
plotstate2=MEFKD.plotstate;
alpha2=MEFKD.p_best(1);

cellstepstotal3=Helacontrol.cellstepstotal;
plotstate3=Helacontrol.plotstate;
alpha3=Helacontrol.p_best(1);


cellstepstotal4=HelaRac1.cellstepstotal;
plotstate4=HelaRac1.plotstate;
alpha4=HelaRac1.p_best(1);


%% use place_predict function to parse state prediction
[cellstepstotal1,state1bycell1,state1forced1,thresh1]=place_predict(cellstepstotal1,plotstate1,alpha1);
[cellstepstotal2,state1bycell2,state1forced2,thresh2]=place_predict(cellstepstotal2,plotstate2,alpha2);
[cellstepstotal3,state1bycell3,state1forced3,thresh3]=place_predict(cellstepstotal3,plotstate3,alpha3);
[cellstepstotal4,state1bycell4,state1forced4,thresh4]=place_predict(cellstepstotal4,plotstate4,alpha4);


%% plot the first box plot
G=[zeros(size(state1bycell1)) ones(size(state1bycell2)) 2*ones(size(state1bycell3))];

bh=boxplot([state1bycell1, state1bycell2, state1bycell3],G,'labels',{'MEF Control','MEF RhoG KD','HeLa Control'});

for k=1:numel(bh(:)) % <- # graphics handles/x
     set(bh(k),'linewidth',3);
end
axes1=gca;
set(axes1,'FontSize',16,'TickLabelInterpreter','none','XGrid','on','XTick',...
    [1 2 3],'YGrid','on');
ylabel('% Time in State One')

pos = [1000         828         625         503];
set(gcf,'position',pos);

cwd=pwd;
cd figures
savefig(gcf,'BoxPlot_unforced');
cd(cwd);

%% plot the second box plot

bh2=boxplot([state1forced1, state1forced2, state1forced3],G,'labels',{'MEF Control','MEF RhoG KD','HeLa Control'});

for k=1:numel(bh2(:)) % <- # graphics handles/x
     set(bh2(k),'linewidth',3);
end
axes2=gca;
set(axes2,'FontSize',24,'TickLabelInterpreter','none','XGrid','on','XTick',...
    [1 2 3],'YGrid','on');
ylabel('% Time in State One')

pos = [1000         828         625         503];
set(gcf,'position',pos);
cwd=pwd;
cd figures
savefig(gcf,'BoxPlot_forced');
cd(cwd);





%% plot MEF control track

h_combined=figure;
flag_forced=1;
[~,cell_plot]=min(abs(50-state1forced1));
cellsteps1=cellstepstotal1(cellstepstotal1(:,3)==cell_plot,:);
state1=plotstate1(cell_plot(1),:);
[h1] = plot_track_with_states(cellsteps1,state1,thresh1,flag_forced,unit_conv1,[0,0],h_combined);
cwd=pwd;
cd figures
savefig(h1,strcat('MEFcontrol_track',num2str(cell_plot)));
cd(cwd);

figure(h_combined);

%% plot MEF KD track
[~,cell_plot]=min(abs(50-state1forced2));
cellsteps2=cellstepstotal2(cellstepstotal2(:,3)==cell_plot,:);
state2=plotstate2(cell_plot,:);
[h2] = plot_track_with_states(cellsteps2,state2,thresh2,flag_forced,unit_conv1,[0,20],h_combined);
cwd=pwd;
cd figures
savefig(h2,strcat('MEFKD_track',num2str(cell_plot)));
cd(cwd);
%% plot HeLa Control track
[~,cell_plot]=min(abs(50-state1forced3));
cellsteps3=cellstepstotal3(cellstepstotal3(:,3)==cell_plot,:);
state3=plotstate3(cell_plot,:);
%

[h3] = plot_track_with_states(cellsteps3,state3,thresh3,flag_forced,unit_conv2,[-113.18,-20],h_combined);

%% plot HeLa Control track
[~,cell_plot]=min(abs(15-state1forced4));
cellsteps4=cellstepstotal4(cellstepstotal4(:,3)==cell_plot,:);
state4=plotstate4(cell_plot,:);
%
cell_plot
size(cellsteps4)
size(plotstate4)
size(state4)
size(cellstepstotal4)
h_sep=figure;
[h_sep2] = plot_track_with_states(cellsteps4,state4,thresh4,flag_forced,unit_conv2,[-113.18,-20],h_sep);

figure(h_sep2)
xlabel('X (\mum)')
ylabel('Y (\mum)')
xlim([-10 150])
ylim([-60 100])
axis square
set(gca,'XTick',[0:20:140])
set(gca,'YTick',[-60:20:100])
cwd=pwd;
cd figures
savestr=strcat('HeLaRac1_track',num2str(cell_plot));
savefig(h3,strcat('HeLaControl_track',num2str(cell_plot)));
savefig(h_sep2,savestr);
cd(cwd);




figure(h_combined)
text(0,7,'MEF Control','FontSize',14,'FontWeight','bold');
text(30,-30,'MEF RhoG KD','FontSize',14,'FontWeight','bold');
text(65,55,'HeLa Control','FontSize',14,'FontWeight','bold');
set(h_combined,'position',[937   655   823   667]);
xlabel('X (\mum)')
ylabel('Y (\mum)')
xlim([-10 150])
ylim([-60 100])
axis square
set(gca,'XTick',[0:20:140])
set(gca,'YTick',[-60:20:100])
cwd=pwd;
cd figures
strcat('HeLaControl_track',num2str(cell_plot));
savefig(h3,strcat('HeLaControl_track',num2str(cell_plot)));
savefig(h_combined,'Combined_track');
save('plotstatepredict_out.mat');
cd(cwd);

close all
%% place_predict subfunction
function [cellstepstotalout,state1bycell,state1_forced,thresh]= place_predict(cellstepstotal,plotstate,alpha)
 cellstepstotalout=cellstepstotal;
 thresh=prctile(plotstate(:),alpha*100);
    for i=1:numel(plotstate(:,1))

        idx=cellstepstotal(:,3)==3;
        
           place=plotstate(i,:)';
           cellstepstotalout(idx,5)=place;
        
           state1bycell(i)=100*sum(round(plotstate(i,:))==1)/numel(plotstate(i,:));
           state1_forced(i)=100*sum(plotstate(i,:)<thresh)/numel(plotstate(i,:));
    end
    
  
    
end


figure
scatter(state1bycell4,state1forced4)


end

