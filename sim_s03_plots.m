clear all
close all
load output/sim/sim_data.mat
load('output/sim/sim_mds.mat');
load('output/sim/sim_cluster.mat')
load_labels
sensations = labels_en;
addpath('/m/nbe/scratch/braindata/shared/toolboxes/cbrewer/');
addpath(genpath('external/'));
sensations_classes 

%% behavioural data correlations
behav_labels={
    'gender' %(F=1)
    'age'
    'mothertongue' 
    'education'
    'psychologist'
    'children'
    'religion' 
    'Items moved'
    'Total time spent on task with pauses'
    'Actual time spent on task'
    }

behav_data(:,3) = [] % remove the native languange column, all 0, can't be correlated
behav_data(:,8) = [] % remove the number of answers column, all 27, can't be correlated


% better labels
better_behav_labels={
'GEN' 
'Age'
'EDU'
'PSY'
'CHL'
'REL' 
'ORI'
'Time1'
'Time2'
}

corrplot(behav_data,'varnames',better_behav_labels, 'type', 'Spearman')

%% do a dendrogram base on sim data
sensations_classes
load_labels
close all
pd=squareform(mean_data_D);
z=linkage(pd,'complete');
[H,T,OUTPERM] = dendrogram(z,0)

axis([0 100 0 0.7])
set(gca,'XTick',[])
for s=1:27
    th=text(s,0,sensations{OUTPERM(s)},'Rotation',90,'BackgroundColor',classColors(classID(OUTPERM(s)),:),'FontSize',6,'HorizontalAlign','right','Color',[1 1 1]);
end

set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color',[1 1 1]);

%% plotting  task data

close all
figure(1)
M = max(max(mean_data_D))
h=imagesc(mean_data_D,[0 0.7]); % maximum distance in mean data is 0.6383
map=cbrewer('seq','Reds',9);
map=[flipud(map);1 1 1];
colormap(map)
hcb = colorbar;
ylabel(hcb, 'Mean distance (a.u.)')
set(h, 'AlphaData', tril(ones(27),0));

axis square

axis off
for n=1:27
    text(n+.5,n-.25,labels{n},'Rotation',45)
end
axis([-10 30 -5 30]) 

%% plot task data resorted by DBSCAN
close all
figure(1)
[aaa bbb]=sort(sim_cluster.DBSCAN.class_from_mean_data,'Descend'); 
FS = 18
h=imagesc(mean_data_D(bbb,bbb),[0 0.7]); 
map=cbrewer('seq','Oranges',9);         
map=[flipud(map);1 1 1];                
colormap(map)                           
hcb = colorbar('Location', 'WestOutside');     

ylabel(hcb, 'Mean distance','FontSize', FS+4)
set(h, 'AlphaData', tril(ones(27),0)); % remove upper triangle

axis square
axis off
for n=1:27
    text(n+.5,n-.25,labels_en{bbb(n)},'Rotation',45, 'FontSize', FS)
end
gtid=find(diff([0;aaa';inf]));
axis([-10 30 -5 30]) 
% prepare top triangle

dbscan_cluster_labels={
  'Love for ideas and non-human animals'
  'Interpersonal love'
  'Between clusters'
};
for g=1:length(gtid)-1
    hold on
    plot([gtid(g) gtid(g)]-.5,[gtid(g) gtid(g+1)]-.5,'k','LineWidth',2)
    templabel=strrep(dbscan_cluster_labels{g},' ',newline);
    text(-1,mean([gtid(g) gtid(g+1)]-.5),templabel,'Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom', 'FontSize', FS)
    plot([gtid(g) gtid(g+1)]-.5,[gtid(g+1) gtid(g+1)]-.5,'k','LineWidth',2)
    plot([-4 gtid(g+1)]-.5,[gtid(g+1) gtid(g+1)]-.5,'Color',[.5 .5 .5]/2,'LineWidth',1)
    plot([gtid(g) gtid(g)]-.5,[gtid(g) 105]-.5,'Color',[.5 .5 .5]/2,'LineWidth',1)
    text(mean([gtid(g) gtid(g+1)]-.5),102,templabel,'HorizontalAlignment','center','VerticalAlignment','top')
end
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color',[1 1 1]);

%print(gcf, ['figs/main/sim_mean_data_resorted_dbscan.png'], '-dpng', '-r300')
print(gcf, ['figs/main/sim_mean_data_resorted_dbscan.jpg'], '-djpeg', '-r400') % journal requirements

%% plot MDS as a network

tsne_out=sim_mds.cmd; % use MDS instead of t-sne

figure(2000)

classColors(2, :) = [0, 205, 108] / 256;
classColors(1,:) = [175, 88, 186] / 256;
classColors(3,:) = 0.3*[1 1 1];       
 
fonttsize = 14;
% network plot
NNeigh=3; % amount of neighbours to be plotted for each node
% first identify range of numbers we have
all_data_temp=[];
for cc=1:27
    temp=mean_data_D(cc,:);    % all distances for the current stimulus
    [aa bb]=sort(temp);        
    bb(1)=[]; % remove itself (index)
    aa(1)=[]; % remove itself (distances)
    for n=1:NNeigh
                    all_data_temp=[all_data_temp; mean_data_D(cc,bb(n))];  % Stores NNeigh closest distances
    end
end

% bin them in 4
th=prctile(all_data_temp,[0 33 66]);
binned_colors=[
       0 0 0

    0.33 0.33 0.33
          0.66 0.66 0.66
    
 
 
];

binned_thicknesses=[
2
1
1
1
];

% plot the lines
for cc=1:27
    temp=mean_data_D(cc,:);
    [aa bb]=sort(temp);
    bb(1)=[]; % remove itself (index)
    aa(1)=[]; % remove itself (distances)
    for n=1:NNeigh
        if(aa(n)<4) % do not plot links that are too distant, it can be removed
            t = linspace(0,1,101); 
            pt1 = [ tsne_out(cc,1);tsne_out(cc,2)];
            pt3 = [tsne_out(bb(n),1); tsne_out(bb(n),2)];
            pt2=(pt1+pt3)/2;
            
            
            %if(sqrt(sum((pt1-pt3).^2))>50)
            %    pt2(1)=1.05*pt2(1);
            %    pt2(2)=1.05*pt2(2);
            %else
            %    pt2=(pt1+pt3)/2;
            %end
           
            pts = kron((1-t).^2,pt1) + kron(2*(1-t).*t,pt2) + kron(t.^2,pt3);

            hold on
            %plot(pts(1,:),pts(2,:),'Color',[1 1 1]*(n-1)/NNeigh) % the problem with this line is that overlapping links have different color
            
            % find bin
            binID=0;
            for bins=1:length(th)
                if(mean_data_D(cc,bb(n))>= th(bins)) 
                    binID=bins;
                end
            end
            
            plot(pts(1,:),pts(2,:),'Color',binned_colors(binID,:),'LineWidth',binned_thicknesses(binID))
            
            
        end
    end
end

% plot the location of the stimuli
for c=1:max(sim_cluster.DBSCAN.class_from_mean_data)  
    if(c<=0) continue; end
   ids=find(sim_cluster.DBSCAN.class_from_mean_data==c);  
   hhh=plot(tsne_out(ids,1),tsne_out(ids,2),'o','Color',classColors(c,:),'MarkerSize',15); 
   set(hhh,'MarkerFaceColor',classColors(c,:)) ;
   hold on
end
ids=find(sim_cluster.DBSCAN.type_from_mean_data==0); 
plot(tsne_out(ids,1),tsne_out(ids,2),'ko','MarkerSize',15) 
ids=find(sim_cluster.DBSCAN.type_from_mean_data==-1);
hhh=plot(tsne_out(ids,1),tsne_out(ids,2),'ko','MarkerSize',15); 
   set(hhh,'MarkerFaceColor',0.3*[1 1 1]); 

% plot labels
for cc=1:27 
    %text(tsne_out(cc,1),tsne_out(cc,2),sensations{cc},'Rotation',0,'Color','Black','FontSize',14,'VerticalAlignment','bottom');
    if cc==25 | cc==16           % overlapping concepts
    text(tsne_out(cc,1)-0.005,tsne_out(cc,2),sensations{cc},'Rotation',0,'FontSize',fonttsize,'HorizontalAlignment','left'); % added a little space between (0.003)
    elseif cc==9          % overlapping concepts
    text(tsne_out(cc,1)-0.006,tsne_out(cc,2)+0.012,sensations{cc},'Rotation',0,'FontSize',fonttsize,'HorizontalAlignment','left'); 
    elseif cc==17          % overlapping concepts
    text(tsne_out(cc,1)-0.001,tsne_out(cc,2)-0.012,sensations{cc},'Rotation',0,'FontSize',fonttsize,'HorizontalAlignment','right'); 
    elseif cc==27          % overlapping concepts
    text(tsne_out(cc,1)-0.001,tsne_out(cc,2)-0.012,sensations{cc},'Rotation',0,'FontSize',fonttsize,'HorizontalAlignment','right'); 
    else
    text(tsne_out(cc,1)+0.005,tsne_out(cc,2),sensations{cc},'Rotation',0,'FontSize',fonttsize,'HorizontalAlignment','right');
    end
end


box off
axis off
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color',[1 1 1]);
set(gca, 'XDir','reverse');                                  % more intuitive layout
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color',[1 1 1]);

% Add legend
a = plot(NaN,NaN,'or', 'MarkerFaceColor', classColors(1,:), 'MarkerEdgeColor', classColors(1,:), 'MarkerSize', 15);
b = plot(NaN,NaN,'ob', 'MarkerFaceColor', classColors(2,:), 'MarkerEdgeColor', classColors(2,:), 'MarkerSize', 15);
c = plot(NaN,NaN,'ok', 'MarkerFaceColor', classColors(3,:), 'MarkerEdgeColor', classColors(3,:), 'MarkerSize', 15);
[hleg1, hobj1] = legend([a b c],{'Interpersonal love', 'Love for ideas and non-human animals', 'Outliers'}, 'FontSize', 14);
set(hleg1,'position',[0.1 0.8 0.3 0.15]);

%export_fig 
%exportgraphics(gcf,'figs/main/sim_network_mds..pdf','ContentType','vector')
%print(gcf, ['figs/main/sim_network_mds.png'], '-dpng', '-r500')
print(gcf, ['figs/main/sim_network_mds.jpg'], '-djpeg', '-r400')

%% Plot HC and K-Means
% Now the chosen clustering is based on the same cluster size as in
% Nummenmaa et al (though K-means suggested 3 cluster solution). 

figure(100)
titlefsize = 18

classColors = cbrewer('qual', 'Set1', 3)

subplot(2,1,1)
% plot the location of the stimuli
for c=1:max(sim_cluster.KMEANS.class_from_mean_data) % DBSCAN /KMEANS / HC
    if(c<=0) continue; end
   ids=find(sim_cluster.KMEANS.class_from_mean_data==c);  % DBSCAN / KMEANS/ HC
   hhh=plot(tsne_out(ids,1),tsne_out(ids,2),'o','Color','k','MarkerSize',13) 
   set(hhh,'MarkerFaceColor',classColors(c,:)) 
   hold on
end

% plot labels
for cc=1:27 
    %text(tsne_out(cc,1),tsne_out(cc,2),sensations{cc},'Rotation',0,'Color','Black','FontSize',14,'VerticalAlignment','bottom');
    if cc==9 | cc==25 | cc==16           % overlapping concepts
    text(tsne_out(cc,1)-0.005,tsne_out(cc,2),sensations{cc},'Rotation',0,'FontSize',12,'HorizontalAlignment','left'); % added a little space between (0.003)
    else
    text(tsne_out(cc,1)+0.005,tsne_out(cc,2),sensations{cc},'Rotation',0,'FontSize',12,'HorizontalAlignment','right');
    end
end

title('k-means clustering','FontSize', titlefsize)
set(gca, 'XDir','reverse')  

% H-Clustering
subplot(2,1,2)
% plot the location of the stimuli
for c=1:max(sim_cluster.HC.class_from_mean_data) % DBSCAN /KMEANS / HC
    if(c<=0) continue; end
   ids=find(sim_cluster.HC.class_from_mean_data==c);  % DBSCAN / KMEANS/ HC
   hhh=plot(tsne_out(ids,1),tsne_out(ids,2),'o','Color','k','MarkerSize',13) 
   set(hhh,'MarkerFaceColor',classColors(c,:)) 
   hold on
end


% plot labels
for cc=1:27 
    %text(tsne_out(cc,1),tsne_out(cc,2),sensations{cc},'Rotation',0,'Color','Black','FontSize',14,'VerticalAlignment','bottom');
    if cc==9 | cc==25 | cc==16           % overlapping concepts
    text(tsne_out(cc,1)-0.005,tsne_out(cc,2),sensations{cc},'Rotation',0,'FontSize',12,'HorizontalAlignment','left'); % added a little space between (0.003)
    else
    text(tsne_out(cc,1)+0.005,tsne_out(cc,2),sensations{cc},'Rotation',0,'FontSize',12,'HorizontalAlignment','right');
    end
end

title('Hierarchichal clustering', 'FontSize', titlefsize)
set(gca, 'XDir','reverse')  
set(gcf,'Units', 'Normalized', 'OuterPosition', [1 1 1 1]);