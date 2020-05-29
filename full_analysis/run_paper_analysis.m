function run_paper_analysis(analyze1,analyze2,short_flag)

%Main function to reproduce code results. First time run will need to use
%analyze1=1, analyze2=1 OR copy the paper results into the results folder.

%Code will take about 24-48 hours to complete in total run on a 2017
%Macbook Pro using parallel processing (things should scale 

if short_flag~=1

%number of bins, larger data-sets with larger variance need bigger n    
n1=35;
n2=120;

% Number of SA runs. Can get away with fewer, but want to try and get close
% to global minimum...
runs=100;

% How many points to check in the confidence interval profiling.
profile_n=40;


else

n1=10;
n2=20;
runs=4;
profile_n=3;

end

% add the source code path
addpath('src')
addpath('data')

if analyze1==1

%note pixel units
prepare_fit_tracks('mef_control','mef_control_analyzed',n1,runs);

[h] = create_polar_histogram('mef_control','MEF_control_polar');

%note pixel units
prepare_fit_tracks('mef_rhog_kd.mat','mef_rhog_kd_analyzed',n1,runs);
[h] = create_polar_histogram('mef_rhog_kd.mat','MEF_RhoGKD_polar');

%note micron units
prepare_fit_tracks('HeLa_control_um','HeLa_control_um_analyzed',n2,runs);
[h] = create_polar_histogram('HeLa_control_um.mat','HelaWT_polar');


%note micron units
prepare_fit_tracks('HeLa_rac1_biosensor','HeLa_rac1_biosensor_analyzed',n1,runs);
[h] = create_polar_histogram('HeLa_rac1_biosensor.mat','HeLa_rac1_biosensor_polar');
end


unit_conv=0.465;



%split on top 5% mur distribution
p_best1=getfitdata('results/mef_control_analyzed.mat',1,3,unit_conv);
%split on top 5% sigmatheta distribution
p_best2=getfitdata('results/mef_rhog_kd_analyzed.mat',4,1,unit_conv);
p_best3=getfitdata('results/HeLa_control_um_analyzed.mat',1,3,1);
p_best4=getfitdata('results/HeLa_rac1_biosensor_analyzed.mat',1.8,1,1);

%analyze2=0;
if analyze2==1
    
clear p_best cellstepstotal fitteddata numb
load('results/mef_control_analyzed.mat')
[p_best_opt1]=profile_likelihood(p_best,cellstepstotal,profile_n,fitteddata,numb,1,0,'CI_mef_control.mat',1);

clear p_best cellstepstotal fitteddata numb
load('results/mef_rhog_kd_analyzed.mat')
[p_best_opt2]=profile_likelihood(p_best,cellstepstotal,profile_n,fitteddata,numb,1,0,'CI_mef_rhog_kd.mat',1);

clear p_best cellstepstotal fitteddata numb
load('results/HeLa_control_um_analyzed')
[p_best_opt3]=profile_likelihood(p_best,cellstepstotal,profile_n,fitteddata,numb,1,0,'CI_HeLa_control.mat',1);

clear p_best cellstepstotal fitteddata numb
load('results/HeLa_rac1_biosensor_analyzed')
[p_best_opt4]=profile_likelihood(p_best,cellstepstotal,profile_n,fitteddata,numb,1,0,'CI_HeLa_rac1_biosensor.mat',0.5);
end

clear fitteddata numb p_best p_singleton val_store
unit_conv=0.465;
load('results/CI_mef_control.mat')
[U_CI_1,L_CI_1,U_CI_2,L_CI_2] = Identify_CI(fitteddata,numb,p_best,p_singleton,val_store,unit_conv, 'MEF Control');

clear fitteddata numb p_best p_singleton val_store
unit_conv=0.465;
load('results/CI_mef_rhog_kd.mat')
[U_CI_1,L_CI_1,U_CI_2,L_CI_2] = Identify_CI(fitteddata,numb,p_best,p_singleton,val_store,unit_conv,'MEF RhoG KD');

clear fitteddata numb p_best p_singleton val_store
load('results/CI_HeLa_control.mat')
[U_CI_1,L_CI_1,U_CI_2,L_CI_2] = Identify_CI(fitteddata,numb,p_best,p_singleton,val_store,1, 'HeLa Control');

clear fitteddata numb p_best p_singleton val_store
load('results/CI_HeLa_rac1_biosensor.mat')
CI=66;
thresh_version=1;
[U_CI_1,L_CI_1,U_CI_2,L_CI_2] = Identify_CI(fitteddata,numb,p_best,p_singleton,val_store,1, 'HeLa Control + Rac1 Biosensor',CI,thresh_version);


%clear all
unit_conv=0.465;


MEFcontrol=load('results/CI_mef_control.mat','cellstepstotal','plotstate','p_best');
MEFKD=load('results/CI_mef_rhog_kd.mat','cellstepstotal','plotstate','p_best');
HeLaControl=load('results/CI_HeLa_control.mat','cellstepstotal','plotstate','p_best');
HeLaRac1=load('results/CI_HeLa_rac1_biosensor.mat','cellstepstotal','plotstate','p_best');

[cellstepstotal1,state1bycell1,state1forced1] = plot_state_predict(MEFcontrol, MEFKD,HeLaControl,HeLaRac1,unit_conv,1);

%Run the foci analysis (pre-loads image analysis componenr from count_foci.m) 
Calculate_Plot_Foci
end






