function  [blurred,cell]=count_foci(~,cellstepstotal,str)
%semi_automated

N = input('How Many Cells Are We Analyzing? ');

dir=pwd;

%loop through the cells
for k=1:N

    folder_name = uigetdir;
    cd(folder_name);

    
Namelist=importdata('Namelist.txt');

Idx=cellstepstotal(:,3)==k;
cellsteps=cellstepstotal(Idx,:);

%loop through the frames (starting at frame 3 which is the first I can
%analyze. Only look at the first ~100 frames because Chris thought the
%cells were sick after this point


for i=cellsteps(1,4):cellsteps(end,4)

%initialize the foci count    
count=0;

% load the frame names
Frame=char(Namelist(i));
% load the frame image
I=imread(Frame);
% set a level for cell mask
level=graythresh(I);%for uint16 background subtracted
BWcell=im2bw(I,level);

%identify the cell and blur it to make nice foci.
BW2=imfill(BWcell, 'holes');
BWcell = +BWcell;
h=fspecial('gaussian', [5 5], 2);
blurred=  imfilter(I, h,'replicate');

%identify the boundaries in the image and find the biggest cell to analyze
    [B,L,N,A]=bwboundaries(BWcell);%,'noholes');
    sizeobject=zeros(N,1);
    
    for i2=1:N
       object=find(L==i2); 
       sizeobject(i2)=numel(object);
    end    
    
   
    cellnumber=find(sizeobject==max(sizeobject));
    
%mask everything apart from the cell
    mask=find(L~=cellnumber);
    cell=L==cellnumber;
    blurred(mask)=0;
    I(mask)=0;
    BW2(mask)=0;
    D=bwdist(~BW2);
    
% plot the cell

    figure(1);
    subplot(1,3,1)
    imagesc(blurred)
    axis square
    
% calculate the average intensity    
average=mean(blurred(cell));

% some intensity levels for plotting
level1=average*1.2;
level2=average*1.3;
level3=average*1.6;


% calculate some properties of the cell (not really used ultimatel)
statsTotal = regionprops(BW2,I,'Area','Eccentricity','Centroid','WeightedCentroid');
FakeR=((statsTotal.Area)^0.5)/pi;
polarity_temp=(norm(statsTotal.Centroid-statsTotal.WeightedCentroid))/FakeR;


% make a thresholded image based on level1
BW_1=im2bw(blurred,level1/2^16);
BW2_1=imfill(BW_1, 'holes');
BW2_1 = +BW2_1;
[B_1,L_1,N_1,A_1]=bwboundaries(BW2_1,'noholes');

% make a thresholded image based on level2
BW_2=im2bw(blurred,level2/2^16);
BW2_2=imfill(BW_2, 'holes');
BW_2 = +BW_2;
[B_2,L_2,N_2,A_2]=bwboundaries(BW2_2,'noholes');

% make a thresholded image based on level3
BW_3=im2bw(blurred,level3/2^16);
BW2_3=imfill(BW_3, 'holes');
BW2_3 = +BW2_3;
[B_3,L_3,N_3,A_3]=bwboundaries(BW2_3,'noholes');
    
% plot the image again to overlay foi boundary
    figure(1)
    
    subplot(1,3,2)
    imagesc(blurred)
    axis square
    
    hold on
    
    % get the foci stats
      stats = regionprops(L_3,'centroid','Area');
      
%{      
for k2 = 1:length(B_1)
   boundary_1 = B_1{k2};
  % if stats(k2).MajorAxisLength> 5  && (stats(k2).MajorAxisLength./stats(k2).MinorAxisLength)<2  %only consider level2 regions bigger than 5 x 2.5 pixels, and not too eliposidal
%   plot(boundary_1(:,2), boundary_1(:,1), 'w', 'LineWidth', 2)
 %  end
end

for k2 = 1:length(B_2)
   boundary_2 = B_2{k2};
  %  if stats(k2).Area> 20%only consider level2 regions bigger than 10 x 5 pixels, and not too eliposidal
  % plot(boundary_2(:,2), boundary_2(:,1), 'r', 'LineWidth', 2)
  %  end
end
%}
      
      
for k2 = 1:length(B_3)
    m=1;
   boundary_3 = B_3{k2};
  
   % calculate the minimal distance between foci and edge.
   mask=L_3~=k2;
   Temp=BW_3;
   Temp(mask)=0;
   Edge_Distance=min(D(logical(Temp)));
    
   % count and plot foci that are greater than 100 pixels, close to the
   % edge and at level3 average intensity
   if stats(k2).Area> 100 && Edge_Distance(1)<5
   plot(boundary_3(:,2), boundary_3(:,1), 'b', 'LineWidth', 2)
   %figure(2)
   % imagesc(blurred)
  %  hold on
 %   axis square
   plot(boundary_3(:,2), boundary_3(:,1), 'b', 'LineWidth', 2)
   %hold off
   count=count+1;
    end
end
hold off

figure(gcf)

%pos=[-1634         518        1661         825];
%set(gcf, 'Position',pos);
count

FociNumber_auto(k,i-cellsteps(1,4)+1)=count;%input('How Many Foci Are There? ');
Eccentricity(k,i-cellsteps(1,4)+1)=statsTotal.Eccentricity;
Polarity(k,i-cellsteps(1,4)+1)=polarity_temp;

%uncomment to do manually, and get the figure

FociNumber_manual(k,i-cellsteps(1,4)+1)=input('How Many Foci Are There? ');
%FociSwitch(k,i-cellsteps(1,4)+1)=input('were the foci stable? ');

%end
subplot(1,3,3)
imagesc(blurred)
axis square

%close all


end
save('temp_foci')
end

cd(dir)
save(str);

end
