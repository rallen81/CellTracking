function [h] = plot_track_with_states(cellsteps,state,thresh,flag_forced,unit_conv,start_pt,h_other)


h=figure;
hold on

x_current=0;start_pt(1);
y_current=0;start_pt(2);

for i=1:numel(state)


x_new=x_current+cellsteps(i,1)*unit_conv;
y_new=y_current+cellsteps(i,2)*unit_conv;

if cellsteps(i,1)==0 && cellsteps(i,2)==0

    x_new=x_current+0.4*(randn-0.5);
    y_new=y_current+0.4*(randn-0.5);
end

if flag_forced==0
if round(state(i))==1
    col='b';
elseif round(state(i))==2
    col='r';
end
end

if flag_forced==1
if state(i)<thresh
    col='b';
elseif state(i)>=thresh
    col='r';
end
end

figure(h_other)
hold on
plot([x_current,x_new]-start_pt(1),[y_current,y_new]-start_pt(2),col,'LineWidth',2);
set(gca,'FontSize',24,'XGrid','on','YGrid','on')
axis square

figure(h)
hold on
plot([x_current,x_new],[y_current,y_new],col,'LineWidth',2);
set(gca,'FontSize',24,'XGrid','on','YGrid','on')
axis square

x_current=x_new;
y_current=y_new;
end
end

