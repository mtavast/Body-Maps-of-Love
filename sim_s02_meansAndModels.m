clear all
close all
addpath('/m/nbe/scratch/braindata/shared/toolboxes/cbrewer/');
addpath(genpath('external/'))

load_labels 

load output/sim/sim_data.mat % loads the data

%% MANTEL TEST ON MEDIAN AND MEAN
% We check that mean and median give same information:

mean_data_D=nanmean(simmat,3); 
median_data_D=nanmedian(simmat,3); 
vsim=nanvar(simmat,[],3);


%% Check histograms of all pairwise distances

ids=find(triu(ones(27),1));             % upper triangle
distancehistogram = zeros(351,length(good_subj_list));     % room for all pairwise distances for all participants
for s=1:351                             % pairwise distances
    for i=1:length(good_subj_list)      % participants
        temp = simmat(:,:,i);           % each participant's distance matrix
        temp2 = temp(ids);              % to a column vector
        distancehistogram(:,i) = temp2; % save
    end
end

figure(1242)
for s=1:117
    subplot(13,9,s);
    histogram(distancehistogram(s,:));      % pairwise distances are in rows, each column a participant
end

figure(1243)
for s=118:234
    subplot(13,9,s-117);
    histogram(distancehistogram(s,:));
end

figure(1244)
for s=235:351
    subplot(13,9,s-234);
    histogram(distancehistogram(s,:));
end


%% do median and mean look very similar? yes

[r_mM p_mM]=bramila_mantel(mean_data_D,median_data_D,5000,'spearman')

r_mM
p_mM
%% split sample reliability

cfg=[];
% reshape the data so that they are one vector per subject
tempids=find(triu(ones(27),1));
cfg.data=zeros(length(tempids),size(simmat,3));
for s=1:size(simmat,3);
    temp=simmat(:,:,s);
    cfg.data(:,s)=temp(tempids);
end

cfg.niter=5000;

[icc permu icc_pvals]=bramila_splitsample(cfg);

% plot 
figure(123)
histogram(icc)
xlabel('Spearman Correlation')
ylabel('Distribution of split-sample correlations across each sensations pair')
% as the pvalues cannot be trust in similarity matrices, we need multiple
% calls to bramila_mantel
Nsubjhalf=floor(size(simmat,3)/2);
for i=1:cfg.niter
    disp(['Mantel test for split sample ' num2str(i)])
    gmat1=nanmean(simmat(:,:,permu(i,1:Nsubjhalf)),3);
    gmat2=nanmean(simmat(:,:,permu(i,(Nsubjhalf+1:end))),3);
    [rtemp(i,1) icc_pvals_mantel(i,1)]=bramila_mantel(gmat1,gmat2,5000,'spearman');
end

save output/sim/split_sample icc permu icc_pvals icc_pvals_mantel rtemp

error('stop')

%% multidimensional scaling

disp('CMDscale')
sim_mds.cmd=cmdscale(mean_data_D,2); 

% MDS 3D
disp('cmd3D 3d')
sim_mds.cmd3D=cmdscale(mean_data_D,3);

save(['output/sim/sim_mds.mat'],'sim_mds','mean_data_D','labels');

%% PCA

[coeff, score, LATENT, TSQUARED, EXPLAINED ] = pca(mean_data_D)

for s=1:27
    figure(203)
    hold on
    plot(score(s,1),score(s,2),'o','MarkerSize',5)
    al='left';
    th=text(score(s,1),score(s,2),labels_en{s},'Rotation',0,'FontSize',12,'HorizontalAlignment',al);
end
box off
axis off
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color',[1 1 1]);
set(gca, 'YDir', 'reverse')
set(gca, 'XDir', 'reverse')

EXPLAINED / sum(EXPLAINED)

%% clustering  using dbscan
addpath('./external/chemometria.us.edu.pl/')
% using spatial distribution from original mean distances
EPS=[];

for minpts_loop=1:10; % optimising numbers of items
    % loop to identify best EPS that maximises number of clusters and
    % minimises outliers
    temp=[];
    for EPS=0:0.01:0.35    % optimising
        [class_from_data,type_from_data]=dbscan_D(mean_data_D,minpts_loop,EPS);
        NCFD=max(class_from_data); % num of clusters from data
        temp=[temp;NCFD EPS];
    end
    % use optimal EPS
    EPS=median(temp(find(temp(:,1)==max(temp(:,1))),2));
    [class_from_mean_data,type_from_mean_data]=dbscan_D(mean_data_D,minpts_loop,EPS);

    numoutliers(minpts_loop)=length(find(class_from_mean_data==-1))
end
minptsOPT=find(min(numoutliers)==numoutliers);

disp('Optimal solution is:')
minpts=minptsOPT

% given minpts, find fine-grained EPS optimal value
for EPS=0:0.001:0.35    % optimising
    [class_from_data,type_from_data]=dbscan_D(mean_data_D,minpts,EPS);
    NCFD=max(class_from_data); % num of clusters from data
    temp=[temp;NCFD EPS];
end
% use optimal EPS
EPS=median(temp(find(temp(:,1)==max(temp(:,1))),2));
[class_from_mean_data,type_from_mean_data]=dbscan_D(mean_data_D,minpts,EPS);
sim_cluster.DBSCAN.class_from_mean_data=class_from_mean_data;
sim_cluster.DBSCAN.type_from_mean_data=type_from_mean_data;
sim_cluster.DBSCAN.minpts=minpts;
sim_cluster.DBSCAN.EPS=EPS;

%% clustering with hierarchical clustering
% We decided beforehand to replicate analysis in Nummenmaa et al (2018).
% Note, in HC and k-means clustering the optimal cluster size is still from that paper, these are just
% for exploration purposes details not changed for this study
z=linkage(squareform(mean_data_D),'complete'); % furthest distance method
for c=2:27
    h_cluster_solution=cluster(z,'MaxClust',c); % constructs a maximum of N clusters using distance as a criterion.  Each node's height in the tree represents the distance between the two subnodes merged at that node.  cluster finds
                                                % the smallest height at which a "horizontal cut" through the tree will leave N or fewer clusters. 
    mhc(c)=min(histc(h_cluster_solution,unique(h_cluster_solution)));
end


opt_c=min(find(diff(mhc>5)==-1)); 

sim_cluster.HC.class_from_mean_data=cluster(z,'MaxClust',opt_c);
sim_cluster.HC.NC=opt_c;

%% k-means clustering
% using spatial distribution from original mean distances
mds_out=sim_mds.cmd;
for NC=2:27
    k_cluster_solution=kmeans(mds_out,NC);
    mkc(NC)=min(histc(k_cluster_solution,unique(k_cluster_solution)));
end

opt_c=min(find(diff(mkc>5)==-1)); % optimal numbers of clusters so that the smallest cluster is of size comparable to DBSCAN
sim_cluster.KMEANS.class_from_mean_data=kmeans(mds_out,opt_c);
sim_cluster.KMEANS.NC=opt_c;

save(['output/sim/sim_cluster.mat'],'sim_cluster');

%% Evaluating k-means
% These seem to suggest 3 cluster solution, however, it uses just the 2D
% MDS solution, and not all the information available from the distances
% Silhouette  
sil = evalclusters(mds_out,'kmeans','silhouette','klist',[2:6])
figure(11)
plot(sil)

% Calinski Harabasz
CH = evalclusters(mds_out,'kmeans','CalinskiHarabasz','klist',[2:6])
figure(12)
plot(CH)

% Davies-Bouldin
DB = evalclusters(mds_out,'kmeans','DaviesBouldin','klist',[2:6])
figure(13)
plot(DB)

%% plotting MDS results

%You can also plot MDS on median values if you want, uncomment from plot
    sim_mds.cmdmedian2 =cmdscale(median_data_D,2);
    mds_out2=sim_mds.cmdmedian2;

mds_out=sim_mds.cmd;

sensations=labels_en;

% do a 2DIM subplot
figure(202)
hold on

for s=1:27
    figure(202)
    hold on
    plot(mds_out(s,1),mds_out(s,2),'o','MarkerSize',7, 'MarkerFaceColor', 'r')
    %plot(mds_out2(s,1),mds_out2(s,2),'o','MarkerSize',7,'MarkerFaceColor', 'b') % uncomment to plot median
    al='right';
    al2='left';
    if s==9 | s==25 | s==16           % overlapping concepts
    text(mds_out(s,1)-0.003,mds_out(s,2),sensations{s},'Rotation',0,'FontSize',12,'HorizontalAlignment',al2); % added a little space between (0.003)
    else
    text(mds_out(s,1)+0.003,mds_out(s,2),sensations{s},'Rotation',0,'FontSize',12,'HorizontalAlignment',al);
    end
    %if s==9 | s==25 | s==16           % overlapping concepts
    %text(mds_out2(s,1)-0.003,mds_out2(s,2),sensations{s},'Rotation',0,'FontSize',12,'HorizontalAlignment',al2); % added a little space between (0.003)
    %else
    %text(mds_out2(s,1)+0.003,mds_out2(s,2),sensations{s},'Rotation',0,'FontSize',12,'HorizontalAlignment',al);
    %end
end
%d = daspect;
box off
axis off
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color',[1 1 1]);
set(gca, 'XDir','reverse')

%% Eigenvalues of MDS results

% for the first 6 dimensions
[Y, E] = cmdscale(mean_data_D,27)

[E E/max(abs(E))]

