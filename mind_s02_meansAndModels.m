clear all
close all
addpath('./external/')
load('output/mind/M10_data')
load_labels;
sensations_classes;
sensations=labels_en;
%% data = words x dimensions x good subjects
% 6 dimensions listed here:
Ndim=6;
dim_labels={ %'Bodily sensation strength','Mind sensation strength','Emotion intensity','Agency','Last time'};
    'Bodily saliency' % How much feels in body question
    'Mental saliency' % How much involves the mind
    'Valence' % emotional valence question
    'Controllability' % how much you can control question (used to be agency)
    'Last experienced' % last time question
    'Touch' % how strongly associate with touch
    };

data(find(data==0))=1; % converts 0 to 1 to avoid problems in binning
data(find(data==1000))=999; % converts 1000 to 999 to avoid inflating the values after the atanh transformation

mean_data=zeros(27,6);
median_data=zeros(27,6);
sem_data=zeros(27,6);
mean_dataZ=zeros(27,6);
median_dataZ=zeros(27,6);
sem_dataZ=zeros(27,6);

% for plotting distros
Nres=100; % resolution
distros=zeros(27,Nres+1,6);
distrosZ=distros;
maph=flipud(cbrewer('div','Spectral',11)); 


% compute means, medians and distributions
for s=1:length(sensations)
    thisdata=squeeze(data(s,:,:)); 
    thisdata=(thisdata-500)/500;   % transform the values for range -1 and 1
    ids=find(~isnan(max(thisdata)));
    thisdata=thisdata(:,ids); % let's keep only those with data (here everyone, as we exluced the ones who don't have data)
    xi=((0:Nres)-Nres/2)/(Nres/2);  
    xiZ=5*(((0:Nres)-(Nres/2))/(Nres/2)); 
    for subdim=1:6
        [fi]=ksdensity(thisdata(subdim,:)',xi); 
        distros(s,:,subdim)=fi/sum(fi);        
    end
    thisdataZ=atanh(thisdata-eps); 
    for subdim=1:6 
        [fi]=ksdensity(thisdataZ(subdim,:)',xiZ);
        %figure(subdim)   % uncomment if you want to compare kdensity plots to histograms 
        %subplot(3,9,s)
        %histogram(thisdataZ(subdim,:), 15);
        %title([labels_en(s) ' D ' num2str(subdim)])
        distrosZ(s,:,subdim)=fi/sum(fi);
    end
    mean_data(s,:)=mean(thisdata,2); 
    mean_dataZ(s,:)=mean(thisdataZ,2); 
    median_data(s,:)=median(thisdata,2); 
    median_dataZ(s,:)=median(thisdataZ,2);
    Nitems(s,1)=size(thisdata,2);
    std_data(s,:)=std(thisdata,0,2); 
    sem_data(s,:)=std(thisdata,0,2)./sqrt(size(thisdata,2));
    sem_dataZ(s,:)=std(thisdataZ,0,2)./sqrt(size(thisdataZ,2));
    sens_count(s,1)=size(thisdata,2);
end

%% visualize case of using normal mean
figure(2)
for subdim=1:6
    subplot(1,6,subdim)
    imagesc(xi,1:27,distros(:,:,subdim),[0 prctile(distros(:),99)])
    colormap(maph)
    set(gca,'YTick',[])
    xlabel(dim_labels{subdim})
    hold on
    for s=1:27
        if(subdim==1)
            text(-1.1,s,sensations{s},'HorizontalAlignment','right')
        end
        if(subdim==6)
            text(1.1,s,sensations{s},'HorizontalAlignment','left')
        end
        plot(mean_data(s,subdim),s,'o','Color',0*[1 1 1])
        plot(median_data(s,subdim),s,'x','Color',0*[1 1 1])
    end
    if(subdim==3)
        title('Probability distribution functions for each item and each dimension. o = mean. x = median')
    end
end
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color',[1 1 1]);
export_fig -m3 'figs/supplement/mind_PDFs_and_mean_data.png'

% case using z transformed means 
figure(3)
for subdim=1:6
    subplot(1,6,subdim)
    imagesc(xiZ,1:27,distrosZ(:,:,subdim))
    colormap(maph)
    set(gca,'YTick',[])
    xlabel([dim_labels{subdim} ' (Z-score)'])
    hold on
    for s=1:27
        if(subdim==1)
            text(-5.1,s,sensations{s},'HorizontalAlignment','right')
        end
        if(subdim==6)
            text(5.1,s,sensations{s},'HorizontalAlignment','left')
        end
        plot(mean_dataZ(s,subdim),s,'o','Color',0*[1 1 1])
        plot(median_dataZ(s,subdim),s,'x','Color',0*[1 1 1])
        
    end
    if(subdim==3)
        title('Probability distribution functions for each item and each dimension. o = mean. x = median')
    end
end
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color',[1 1 1]);
%% label sorting with dim 3 (emotions) for mean_data and mean_dataZ
[temp emomedian]=sort(median_data(:,3)); % dimension to be sorted by

emomean_mean_data=mean_data(emomedian,:);
emomean_median_data=median_data(emomedian,:);
emomean_distros=distros(emomedian,:,:);
emomean_sensations=sensations(emomedian);
%visualize case of using normal mean
figure(22)
for subdim=1:6
    subplot(1,7,subdim)
    imagesc(xi,1:27,emomean_distros(:,:,subdim),[0 prctile(emomean_distros(:),99)])
    colormap(maph)
    set(gca,'YTick',[])
    xlabel(dim_labels{subdim})
    hold on
    for s=1:27
        if(subdim==1)
            text(-1.1,s,emomean_sensations{s},'HorizontalAlignment','right')
        end
        if(subdim==6)
            text(1.1,s,emomean_sensations{s},'HorizontalAlignment','left')
        end
        plot(emomean_mean_data(s,subdim),s,'o','Color',0*[1 1 1])
        plot(emomean_median_data(s,subdim),s,'x','Color',0*[1 1 1])
    end
    if(subdim==3)
        title('Probability distribution functions for each item and each dimension. o = mean. x = median')
    end
end
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color',[1 1 1]);
%% visualize as violin plots using z  
% article version
load output/sim/sim_cluster.mat

classColors(2, :) = [0, 205, 108] / 256
classColors(1,:) = [175, 88, 186] / 256
classColors(3,:) = 0.3*[1 1 1];

sensations_shortened = {'self',                  
    'parental',              
    'romantic',             
    'friends',     
    'neighbours',  
    'of beauty',             
    'country',     
    'for God',               
    'of nature',             
    'universal',             
    'sexual',                
    'moral',                 
    'siblings',    
    'practical',             
    'of wisdom',                 
    'mother’s',
    'father’s',
    'true',                  
    'reciprocal',            
    'unconditional',         
    'companionate',          
    'passionate',            
    'altruistic',    
    'benevolent',            
    'for life',              
    'strangers',         
    'animals'} 

Fsize = 14
Fsize2 = 18

figure(22001)
for subdim=1:6
    if subdim==1
        subplot_1 = subplot(1,6,subdim)
    elseif subdim==2
        subplot_2 = subplot(1,6,subdim)
    elseif subdim==3
        subplot_3 = subplot(1,6,subdim)
    elseif subdim==4
        subplot_4 = subplot(1,6,subdim)
    elseif subdim==5
        subplot_5 = subplot(1,6,subdim)
    elseif subdim==6
        subplot_6 = subplot(1,6,subdim)
    end
    % to make the subplots and texts fit
    
    temp_distros = distrosZ(:,:,subdim); % temporary distributions for the given dimension
    if subdim==5
        [temp_median temp_sensations_order]=sort(median_dataZ(:,subdim), 'descend'); % last experienced has reversed order
    else
        [temp_median temp_sensations_order]=sort(median_dataZ(:,subdim)); % temp_sensations = order of the stimuli
    end
    temp_distros_sorted = temp_distros(temp_sensations_order, :);
    if subdim==1
        temp_sensations = sensations(temp_sensations_order);   % reorder the labels
    else
        temp_sensations = sensations_shortened(temp_sensations_order);   % reorder the labels, shortened versions
    end
    
    for s=27:-1:1
        p=patch([xiZ'; flipud(xiZ')],[s+(squeeze(30*temp_distros_sorted(s,:)))'; s*ones(length(xiZ),1)] ,1);
        colID=sim_cluster.DBSCAN.class_from_mean_data(temp_sensations_order(s));
        if(colID<=0)
            col=classColors(end,:);
        else
            col=classColors(colID,:);
        end
        set(p,'EdgeColor',[1 1 1]*.9)
        set(p,'FaceColor',col);
        hold on
    
        % labels
    th=text(-1*(4.1),s+.25,temp_sensations{s},'HorizontalAlignment','right','FontSize',Fsize);

    end
    if(subdim==3)
        %title('Probability distribution functions for each item and each dimension. o = mean. x = median')
    end
    set(gca,'YTick',[])
    xlabel(dim_labels{subdim}, 'FontSize', 18)
    axis([-4 4 0 30])
    axis off
    p  = patch([ 0 0],[0 101],1,'linestyle',':','edgecolor',[0 0 0]+.15,'linewidth',.3,'edgealpha',0.15); % vertical line
    text(0,-.4,'0','HorizontalAlignment','center','FontSize',Fsize2);
    text(4,-.4,'4','HorizontalAlignment','center','FontSize',Fsize2);
    text(-4,-.4,'-4','HorizontalAlignment','center','FontSize',Fsize2);
    text(0,-2.6,dim_labels{subdim},'HorizontalAlignment','center','FontSize',Fsize2);
end

subplot_1.Position = subplot_1.Position - [0. 0 0.028 0] % decrease width
subplot_1.Position = subplot_1.Position + [0.022 0 0 0] % push little bit on the right   

subplot_2.Position = subplot_2.Position - [0 0 0.028 0] % dexrease width
subplot_2.Position = subplot_2.Position + [0.035 0 0 0] % push little bit on the right           

subplot_3.Position = subplot_3.Position - [0 0 0.028 0] % dexrease width
subplot_3.Position = subplot_3.Position + [0.05 0 0 0] % push little bit on the right

subplot_4.Position = subplot_4.Position - [0 0 0.028 0] % dexrease width
subplot_4.Position = subplot_4.Position + [0.072 0 0 0] % push little bit on the right 

subplot_5.Position = subplot_5.Position - [0 0 0.028 0] % dexrease width
subplot_5.Position = subplot_5.Position + [0.085 0 0 0] % push little bit on the right

subplot_6.Position = subplot_6.Position - [0 0 0.028 0] % dexrease width
subplot_6.Position = subplot_6.Position + [0.1 0 0 0] % push little bit on the right 


h = gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

print(gcf, 'figs/main/distributions_sorted.png', '-dpng', '-r500')
exportgraphics(gcf,'figs/main/distributions_sorted.pdf','ContentType','vector')

%% visualize as violin plots 
load output/sim/sim_cluster.mat
figure(2200)
close
figure(2200)
for subdim=1:6
    subplot(1,6,subdim) %
    
    
    for s=27:-1:1 % käänteisessä järjestyksessä
        p=patch([xi'; flipud(xi')],[s+(squeeze(30*emomean_distros(s,:,subdim)))'; s*ones(length(xi),1)] ,1);
        colID=sim_cluster.DBSCAN.class_from_mean_data(emomedian(s));
        if(colID<=0)
            col=classColors(end,:); % musta, rivillä 6
        else
            col=classColors(colID,:); % muut värit
        end
        set(p,'EdgeColor',[1 1 1]*.9)
        set(p,'FaceColor',col);
        hold on
        if(subdim==1)
            text(-1.1,s+.25,emomean_sensations{s},'HorizontalAlignment','right','FontSize',6)
        end
        if(subdim==6)
            text(1.1,s+.25,emomean_sensations{s},'HorizontalAlignment','left','FontSize',6)
        end

    end
    if(subdim==3)
        title('Probability distribution functions for each item and each dimension. o = mean. x = median')
    end
    set(gca,'YTick',[])
    xlabel(dim_labels{subdim})
    axis([-1 1 0 30])
    axis off
    text(0,0,'500','HorizontalAlignment','center')
    text(1,0,'1000','HorizontalAlignment','center')
    text(-1,0,'0','HorizontalAlignment','center')
    text(0,-2,dim_labels{subdim},'HorizontalAlignment','center')
end
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color',[1 1 1]);
%% case using z transformed means
[temp emomedianZ]=sort(median_dataZ(:,3));

emomean_mean_dataZ=mean_dataZ(emomedianZ,:);
emomean_median_dataZ=median_dataZ(emomedianZ,:);
emomean_distrosZ=distrosZ(emomedianZ,:,:);
emomean_sensationsZ=sensations(emomedianZ);

figure(33)
for subdim=1:6
    subplot(1,6,subdim)
    imagesc(xiZ,1:27,emomean_distrosZ(:,:,subdim))
    colormap(maph)
    set(gca,'YTick',[])
    xlabel([dim_labels{subdim} ' (Z-score)'])
    hold on
    for s=1:27
        if(subdim==1)
            text(-5.1,s,emomean_sensationsZ{s},'HorizontalAlignment','right','FontSize',6)
        end
        if(subdim==6)
            text(5.1,s,emomean_sensationsZ{s},'HorizontalAlignment','left','FontSize',6)
        end
        plot(emomean_mean_dataZ(s,subdim),s,'o','Color',0*[1 1 1])
        plot(emomean_median_dataZ(s,subdim),s,'x','Color',0*[1 1 1])
        
    end
    if(subdim==3)
        title('Probability distribution functions for each item and each dimension. o = mean. x = median')
    end
end
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color',[1 1 1]);
%[sensations(emomedian) sensations(emomedianZ)]

%% visualize as violin plots using z  

load output/sim/sim_cluster.mat

[temp emomedianZ]=sort(median_dataZ(:,3));

emomean_mean_dataZ=mean_dataZ(emomedianZ,:);
emomean_median_dataZ=median_dataZ(emomedianZ,:);
emomean_distrosZ=distrosZ(emomedianZ,:,:);
emomean_sensationsZ=sensations(emomedianZ);

Fsize = 13
Fsize2 = 14

figure(22001)
for subdim=1:6
    subplot(1,6,subdim) 
    
    for s=27:-1:1
        p=patch([xiZ'; flipud(xiZ')],[s+(squeeze(30*emomean_distrosZ(s,:,subdim)))'; s*ones(length(xiZ),1)] ,1);
        colID=sim_cluster.DBSCAN.class_from_mean_data(emomedianZ(s));
        if(colID<=0)
            col=classColors(end,:);
        else
            col=classColors(colID,:);
        end
        set(p,'EdgeColor',[1 1 1]*.9)
        set(p,'FaceColor',col);
        hold on
        switch(4)
            case 4 % labels everywhere without alternatinv
                if(subdim==1)
                    th=text(-1*(4.1),s+.25,emomean_sensationsZ{s},'HorizontalAlignment','right','FontSize',Fsize)
                end              
        end
    end
    if(subdim==3)
        %title('Probability distribution functions for each item and each dimension. o = mean. x = median')
    end
    set(gca,'YTick',[])
    xlabel(dim_labels{subdim})
    axis([-4 4 0 30])
    axis off
    p  = patch([ 0 0],[0 101],1,'linestyle',':','edgecolor',[0 0 0]+.15,'linewidth',.3,'edgealpha',0.15); % vertical line
    text(0,-.4,'0','HorizontalAlignment','center','FontSize',Fsize2)
    text(4,-.4,'4','HorizontalAlignment','center','FontSize',Fsize2)
    text(-4,-.4,'-4','HorizontalAlignment','center','FontSize',Fsize2)
    text(0,-2.6,dim_labels{subdim},'HorizontalAlignment','center','FontSize',Fsize2)
end

set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color',[1 1 1]);

%print(gcf, 'figs/main/distributions_clustered.png', '-dpng', '-r300')


%% Multidimensional scaling
% PCA 2D
pcadata=[(mean_dataZ)];
[coeff, score, latent, tsquared, explained] = pca(pcadata);

[L1, T]=rotatefactors(coeff(:,1:2));
score_rot=pcadata*L1;
mind_mds.pca2D=score_rot;
% cmdscale function
mind_mds.cmd=cmdscale(squareform(pdist(mean_data)),2); % from the non-transformed mean data
% PCA 3D
pcadata=[mean_dataZ];
[coeff, score, latent, tsquared, explained] = pca(pcadata);
[L1, T]=rotatefactors(coeff(:,1:3));
mind_mds.pca3D=pcadata*L1;

save(['output/mind/mind_mds.mat'],'mind_mds','mean_data','mean_dataZ','median_data','median_dataZ','sensations','classID','sens_count');

error('stop here before plots')


%% assessing reliability of the tool
rng(0)
cfg=[];
for dim=1:6
    cfg.data=squeeze(atanh((data(:,dim,:)-500)/500));
    cfg.niter=5000;
    if(dim>1)
        cfg.permu=permu; % use same permutations
    end
    [icc(:,dim) permu pvals(:,dim)] = bramila_splitsample(cfg);
end
%% plotting reliabilities

for dim=1:6                         
    figure(100)                     
    subplot(2,6,dim)                
    histogram(icc(:,dim))            
    set(gca,'XLim',[0 1])           
    xlabel('Spearman correlation')  
    ylabel('Number of iterations')  
    title(dim_labels{dim})          
    axis square                     
    subplot(2,6,dim+6)              
    plot(icc(:,dim),-log10(pvals(:,dim)),'.')  
    axis([0 1 0 40])                           
    yticks=get(gca,'YTick');                  
    yticks=yticks(1:end)';                    
    set(gca,'YTick',yticks)                   
    set(gca,'YTickLabel',num2str(10.^-yticks))  
    xlabel('Spearman correlation')
    ylabel('P-value')
    title(['Largest p-value: ' num2str(max(pvals(:,dim)),3)]) %
    axis square                                               
    grid on
end

set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 .5]);
set(gcf,'color',[1 1 1]);
print(gcf, 'figs/supplement/reliability_of_rating_tool.png', '-dpng', '-r500')

disp(['Median split-half correlation value', num2str(round(median(icc), 2))])
