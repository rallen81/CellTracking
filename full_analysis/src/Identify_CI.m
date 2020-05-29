function [U_CI_1,L_CI_1,U_CI_2,L_CI_2] = Identify_CI(fitteddata,numb,p_best,p_singleton,val_store,unit_conv,dataset,varargin)

% This function takes the variables generated from profile_likelihood and
% calcuates the CI by checking the intercept between some threshold (by
% bootstrapping) and the SSE of the fit along each parameter axis

[CI_5,CI_95,SSE,SSE2] = bootstrapsse(fitteddata,numb);
str={'\alpha', '\mu_1 (\mum)','\sigma_{r_1} (\mum)','\sigma_{{\theta}_1} (radians)','\mu_2 (\mum)','\sigma_{r_2} (\mum)','\sigma_{{\theta}_2} (radians)'};

if isempty(varargin)==0
thresh1=prctile(SSE,varargin{1}); %resample cells
thresh2=prctile(SSE2,varargin{1}); %resample steps
else
thresh1=prctile(SSE,95);
thresh2=prctile(SSE2,95);
end

if isempty(varargin)==0
if varargin{2}==1;
    thresh1=thresh1;
else
    thresh1=thresh2;
end
end
thresh1
figure

%% do the interpolation if possible, and make subplots
for i=1:numel(p_best)
    L_CI_1(i)=0;
    L_CI_2(i)=0;
    U_CI_1(i)=max(p_singleton(i,:));
    U_CI_2(i)=max(p_singleton(i,:));
   
    [curve,gof]=fit(p_singleton(i,:)',val_store(i,:)','pchipinterp');
    
    %%upperCI
    
    x0=p_best(i)*1.1;
    
    obj1=@(x)curve(x)-thresh1;
    obj2=@(x)curve(x)-thresh2;
    xplot=min(p_singleton(i,:)):0.01:max(p_singleton(i,:));
    yplot=curve(xplot);
    
   
    
    
    
    Dist1=1e256;
    Dist2=1e256;
    Dist3=1e256;
    Dist4=1e256;
    
    for k=1:numel(p_singleton)
        
        if p_singleton(k)<p_best(i)
            
           
    tempL=fzero(obj1,p_singleton(k));
    temp=abs(tempL-p_best(i));
    
    tempL2=fzero(obj2,p_singleton(k));
    temp2=abs(tempL2-p_best(i));
    
    if temp<Dist1 && tempL<p_best(i)
        if tempL<0;
            tempL=0;
        end
    L_CI_1(i)=tempL;
    Dist1=temp;
    end
    
    if temp2<Dist2 && tempL2<p_best(i)
        if tempL2<0;
            tempL2=0;
        end
    L_CI_2(i)=tempL2;
    Dist1=temp2;
    end
    
        end
        
        
        if p_singleton(k)>=p_best(i)
            tempL3=fzero(obj1,p_singleton(k));
            temp3=abs(tempL3-p_best(i));
   
    tempL4=fzero(obj2,p_singleton(k));
    temp4=abs(tempL4-p_best(i));
    
    if temp3<Dist3 && tempL3>p_best(i)
    U_CI_1(i)=tempL3;
    Dist3=temp3;
    end
    
    if temp4<Dist4 && tempL4>p_best(i)
    U_CI_2(i)=tempL4;
    Dist4=temp4;
    end
    
        end
    end
    
     subplot(7,1,i);
     if i==2 || i==3 || i==5 || i==6
     unit_conv2=unit_conv;
     else
         unit_conv2=1;
     end
     plot(xplot*unit_conv2,yplot,'Linewidth',3);
      hold on
     plot(xplot*unit_conv2,repmat(thresh1,size(xplot)),':r','Linewidth',1);
   %  plot(xplot,repmat(thresh2,size(xplot)),':k','Linewidth',1);
    
     yThresh=[0 max(val_store(i,:))];
  
     plot([max(L_CI_1(i),0) max(L_CI_1(i),0)]*unit_conv2,yThresh,'--r','Linewidth',2)
   %  plot([max(L_CI_2(i),0) max(L_CI_2(i),0)],yThresh,'--k','Linewidth',2)  
     
     plot([U_CI_1(i) U_CI_1(i)]*unit_conv2,yThresh,'--r','Linewidth',2)
  %   plot([U_CI_2(i) U_CI_2(i)],yThresh,'--k','Linewidth',2)  
     
     ylim([0,max(thresh1)*1.5]);
     xlim([0,p_singleton(i,end)*unit_conv2]);
     xlabel(str(i))
     ylabel('SSE')
     set(gca,'FontSize',16);
   %  if i==1
   %      title(dataset)
   %  end
end
   
%% clean up figure and save
figure(gcf)

ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 
1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');

text(0.5, 1,dataset,'HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',24,'FontWeight','bold');


set(gcf,'units','normalized','outerposition',[0.3000    0.0917    0.1113    0.7799])
cwd=pwd;
cd figures
strcat(dataset,'_p_summary');
savefig(gcf,strcat(dataset,'_p_summary'));
cd(cwd);
%title(dataset)
close all
end

