clear all
close all
addpath(genpath('external/'))
addpath(genpath('external/colorbrewer'))

% data from the similarity task
simdata=load('output/sim/sim_mds.mat');   % Euclidean distance to begin with. Uses the mean distance.

% data from the dimension rating task
minddata=load('output/mind/mind_mds.mat'); %  non-transformed data

% data from the body painting task
bodydata=load('output/body/body_mds.mat'); % uses the euclidean distance from effect sizes (mean_data is effect sizes)

reorder = false

alldimlabels={
'Feeling similarity'
'Bodily saliency' % How much feels in body question
'Mental saliency' % How much involves the mind
'Valence' % emotional valence question
'Controllability' % how much you can control question (used to be agency)
'Last experienced' % last time question
'Touch' % how strongly associate with touch
'Bodily maps'
}


%% Preparing the distance matrices

ids=find(triu(ones(27),1)); % top triangle

Y_dm=simdata.mean_data_D;       % similarity data is already in a distance matrix form
all_dm(:,:,1)=Y_dm;             % add to all_dm
Y=zscore(Y_dm(ids));            % convert distances to zscores 
X1=[];
% X1 mind
for mr=1:6 % 
    X1_dm(:,:,mr)=squareform(pdist(minddata.mean_data(:,mr)));     % means of all dimensions, euclidean distances
    temp=X1_dm(:,:,mr);
    X1=[X1 zscore(temp(ids))]; % to zscores
end
all_dm(:,:,2:7)=X1_dm;

% Add distances of the body maps, we use euclidean distance
X2_dm_1=squareform(pdist(bodydata.mean_data));
X2_dm=X2_dm_1;
X2=[zscore(X2_dm_1(ids))]; 

all_dm(:,:,8)=X2_dm_1;

X=[X1 X2];

%% Mantel between each pair of distance matrices

tempdata=[Y X];  % sim zscores as a first column

if reorder == true
    tempdata = [tempdata(:,1) tempdata(:,3:4) tempdata(:,2) tempdata(:,7:8) tempdata(:,5:6)] % reordering to make the figure easier to read
end

rda_r=zeros(size(tempdata,2)); 
rda_p=zeros(size(tempdata,2));
% run mantel to test significance

for d1=1:size(tempdata,2)
        
    for d2=(d1+1):size(tempdata,2)
        idsOI=(1:351)';
        SS=27;
        mat1=zeros(SS);
        mat1(find(triu(ones(SS),1)))=tempdata(idsOI,d1);
        mat1=mat1+mat1';   % fills the lower triange
            
        mat2=zeros(SS);
        mat2(find(triu(ones(SS),1)))=tempdata(idsOI,d2);
        mat2=mat2+mat2';   
        rng(0); % same permutations for each pair of dimensions
        [rtemp ptemp]=bramila_mantel(mat1,mat2,5000,'spearman')
        rda_r(d1,d2)=rtemp;     % input the result of the mantel test
        rda_p(d1,d2)=ptemp;
    end
end

rda_r=rda_r+rda_r'+eye(size(rda_r));  

q=mafdr(rda_p(find(triu(ones(size(rda_p)),1))),'BHFDR','true');
mask_p=zeros(size(rda_p));
mask_p(find(triu(ones(size(rda_p)),1)))=q<0.05;

mask_p(find(triu(ones(size(rda_p)),1)))=mask_p(find(triu(ones(size(rda_p)),1)))+(q<0.005); 



%% Make the RDA figure
if reorder == true
    alldimlabels = [alldimlabels(1) alldimlabels(3) alldimlabels(4) alldimlabels(2) alldimlabels(7) alldimlabels(8) alldimlabels(5) alldimlabels(6)] % for reordering the labels
end

tickfontsize = 18;

figure(1)
imagesc(rda_r,[min(min(rda_r)) 1]) 
hold on
% add a plot of stars
for d1=1:size(mask_p,2)
    for d2=(d1+1):size(mask_p,2)
        if(mask_p(d1,d2)==1) text(d1,d2,'*','FontSize',20,'HorizontalAlignment','center'); end
        if(mask_p(d1,d2)==2) text(d1,d2,'**','FontSize',20,'HorizontalAlignment','center'); end
        if(mask_p(d1,d2)==1) text(d2,d1,'*','FontSize',20,'HorizontalAlignment','center'); end
        if(mask_p(d1,d2)==2) text(d2,d1,'**','FontSize',20,'HorizontalAlignment','center'); end
    end
end

map = brewermap(16, 'YlOrRd')
%map=cbrewer('seq','YlOrRd',16);
map=[0.5020/2         0    0.1490/2;flipud(map);1 1 1;];
%map= map + abs(min(map,[], 'all')) 
%map = map ./ max(map(:));
colormap(map) % change colors
hcb=colorbar
ylabel(hcb,'Similarity (Spearman''s correlation)')
set(gca,'YTick',1:length(alldimlabels))
set(gca,'YTickLabel',alldimlabels,'FontSize', tickfontsize)
set(gca,'XTick',1:length(alldimlabels))
set(gca,'XTickLabel',alldimlabels,'XTickLabelRotation',90)
axis square
set(gcf,'Color','white')
set(gcf,'Units', 'Pixels', 'OuterPosition', [0 0 1300 1000]);

print(gcf, 'figs/supplement/RSA.png', '-dpng', '-r400')



