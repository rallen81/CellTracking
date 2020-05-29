function [fakeD,h] = bootstrap_state_data(Data,state1_logical,state2_logical)


proportionstate1=sum(state1_logical)/numel(state1_logical);

realD=mean(Data(state1_logical))-mean(Data(state2_logical));

state1_fake=logical(size(state1_logical));

n=50000;
fakeD=zeros(n,1);
for i=1:n
   
    for j=1:numel(Data)
    p=rand;
    
    if p<proportionstate1
    
        state1_fake(j)=1;
    else 
        state1_fake(j)=0;
        
    end
    
    end
    
    fakeD(i)=mean(Data(state1_fake))-mean(Data(~state1_fake));
end

figure
h=histogram(fakeD,25)
hold on

h=plot([realD realD],[0 max(h.Values)],'--','LineWidth',2)
set(gca,'FontSize',16);
xlabel('\Delta mean foci')
ylabel('Frequency Density')

cwd=pwd;
cd figures
savefig(gcf,'bootstrap_state');
cd(cwd);

end

