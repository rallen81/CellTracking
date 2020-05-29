%load the data- specifically my state prediction
%load('CI_HeLaWTRac1Fret_final')

%str='foci_save';
%count_foci(plotstate,cellstepstotal,str);

%all foci analyzed is my final analysis, uncomment out one of these to load any
%repeated analysis.
load('all foci analyzed')
%load('foci_save')


include=[1 2 3 4 5 6];
plotstate2=plotstate(include,:);
FociNumber_auto=FociNumber_auto(include,:);

A=round(plotstate2);

state1=A==1;
state2=A==2;



figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.11 0.347700693756194 0.815]);
scatter(round(plotstate2(:))+0.1*rand(numel(plotstate2(:)),1)-0.05,FociNumber_auto(:))
hold on
errorbar([1 2],[mean(FociNumber_auto(state1)) mean(FociNumber_auto(state2))],[std(FociNumber_auto(state1)) std(FociNumber_auto(state2))],'s','MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',3 )

% Create ylabel
ylabel('Number of Foci Identified in Frame');

% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0 7]);
% Set the remaining axes properties
set(axes1,'FontSize',16,'XGrid','on','XTick',[1 2],'XTickLabel',...
    {'State One','State Two'},'YGrid','on');
ylim([0 7])

cwd=pwd;
cd figures
savefig(gcf,'bootstrap_state_scatter');
cd(cwd);


%hold on
%histogram(FociNumber_auto(state1),'Normalization','probability') ;
%hold on
%histogram(FociNumber_auto(state2),'Normalization','probability') ;
[fakeD,h] = bootstrap_state_data(FociNumber_auto(:),state1(:),state2(:));

