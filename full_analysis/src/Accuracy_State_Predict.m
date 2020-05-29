load('CI_sheet1.mat')
[n,xout,A,controldata,coord,state_simulated,cellstepstotal_simulated]=cellmotionmodel(p_best);

[plotstate,state_predicted,state_forced]=predictstates_fun(p_best,cellstepstotal_simulated);
%state_predicted=state_forced;
acc=state_simulated-state_predicted;
correct=acc==0;
accuracy_MEF_control=sum(correct)/numel(acc)

state1=state_simulated==1;
acc=state_simulated(state1)-state_predicted(state1);
correct=acc==0;
accuracy_MEF_control_s1=sum(correct)/numel(acc)

acc=state_simulated(~state1)-state_predicted(~state1);
correct=acc==0;
accuracy_MEF_control_s2=sum(correct)/numel(acc)


clear p_best;

load('CI_sheet2.mat')


[n,xout,A,controldata,coord,state_simulated,cellstepstotal_simulated]=cellmotionmodel(p_best);

[plotstate,state_predicted,state_forced]=predictstates_fun(p_best,cellstepstotal_simulated);
%state_predicted=state_forced;


acc=state_simulated-state_predicted;
correct=acc==0;

accuracy_MEF_RhoGKD=sum(correct)/numel(acc)


state1=state_simulated==1;
acc=state_simulated(state1)-state_predicted(state1);
correct=acc==0;
accuracy_MEF_RhoGKD_s1=sum(correct)/numel(acc)

acc=state_simulated(~state1)-state_predicted(~state1);
correct=acc==0;
accuracy_MEF_RhoGKD_s2=sum(correct)/numel(acc)


clear p_best;

load('CI_RhoG_CW_WT.mat')

[n,xout,A,controldata,coord,state_simulated,cellstepstotal_simulated]=cellmotionmodel(p_best);

[plotstate,state_predicted,state_forced]=predictstates_fun(p_best,cellstepstotal_simulated);
%state_predicted=state_forced;

acc=state_simulated-state_predicted;
correct=acc==0;

accuracy_HeLaWT=sum(correct)/numel(acc)

state1=state_simulated==1;
acc=state_simulated(state1)-state_predicted(state1);
correct=acc==0;
accuracy_HeLaWT_s1=sum(correct)/numel(acc)

acc=state_simulated(~state1)-state_predicted(~state1);
correct=acc==0;
accuracy_HeLaWT_s2=sum(correct)/numel(acc)


clear p_best;


load('CI_HeLaWTRac1Fret_final.mat')

[n,xout,A,controldata,coord,state_simulated,cellstepstotal_simulated]=cellmotionmodel(p_best);

[plotstate,state_predicted,state_forced]=predictstates_fun(p_best,cellstepstotal_simulated);

%state_predicted=state_forced;

acc=state_simulated-state_predicted;
correct=acc==0;

accuracy_HeLaRac1Fret=sum(correct)/numel(acc)

state1=state_simulated==1;
acc=state_simulated(state1)-state_predicted(state1);
correct=acc==0;
accuracy_HeLaRac1Fret_s1=sum(correct)/numel(acc)

acc=state_simulated(~state1)-state_predicted(~state1);
correct=acc==0;
accuracy_HeLaRac1Fret_s2=sum(correct)/numel(acc)

clear p_best;
