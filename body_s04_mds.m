clear all
close all

addpath('./external/bodyspm')
load('output/body/bspm_ttest.mat');
mask=uint8(imread('bodySPM_base3.png'));
in_mask=find(mask>128); % list of pixels inside the mask
base=uint8(imread('bodySPM_base2.png'));
base2=base(10:531,33:203,:); % single image base

load_labels % makes variable labels

mean_data=zeros(length(in_mask),27); % will contain effect sizes vectors
for n=1:27
	temp=bspm.ttest.es(:,:,n);    % uses the effect sizes
	mean_data(:,n)=temp(in_mask); % only pixel in masks
end
mean_data=mean_data'; 

%% Multidimensional scaling

disp('CMDscale')
body_mds.cmd=cmdscale(squareform(pdist(mean_data)),2); % Euclidean distance

%plot
mds_out=body_mds.cmd;

for s=1:27
    figure(2000)
    hold on
    plot(mds_out(s,1),mds_out(s,2),'o','MarkerSize',5)
    al='left';
    th=text(mds_out(s,1),mds_out(s,2),labels_en{s},'Rotation',0,'FontSize',12,'HorizontalAlignment',al);
end

save(['output/body/body_mds.mat'],'body_mds','mean_data','labels');

