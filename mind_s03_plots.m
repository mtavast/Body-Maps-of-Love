clear all
close all
addpath(genpath('external'));
load('output/mind/M10_data');
load('output/mind/mind_mds.mat');
sensations_classes
load_labels

%% corr plot for the mean_dataZ to check it match mycorr plot
for i=1
figure(123456)
corrplot(mean_dataZ,'varnames',dim_labels,'type','Spearman','testR','on', 'alpha', 0.05)
hcp=gcf;
hcpC=get(hcp,'Children');
redmap=cbrewer('seq','Reds',9);
qualmap=cbrewer('qual','Set1',9);
set(gcf, 'OuterPosition', [0 0 900 900]);
for han=1:length(hcpC)
    temp=get(hcpC(han),'XLabel');
    set(temp,'Color',[0 0 0])
    temp=get(hcpC(han),'Children');
    if(strcmp(class(temp),'matlab.graphics.chart.primitive.Histogram'))
        set(temp,'FaceColor',qualmap(1,:))
        set(temp,'FaceAlpha',1)
        set(temp,'EdgeColor',redmap(end,:))
    else
        for tt=1:length(temp)
            set(hcpC(han),'XLim',[-3.5 3.5])
            set(hcpC(han),'YLim',[-3.5 3.5])
            set(temp(tt),'Color',redmap(end,:))
            set(temp(tt),'MarkerEdgeColor',qualmap(1,:))
        end
    end
end
set(gcf,'color',[1 1 1]);
end
%% myCorr plot 
%close all

FS=12.5; % font size
DF=8; % distance between labels
sim=load('output/sim/sim_cluster.mat')
load output/mind/mind_mds.mat

dbscan_cluster_labels={
    'Interpersonal love'
    'Love for ideas and non-human animals'
    'Outliers'
    };

figure(1)
data=mean_dataZ;
N=size(data,2); % dimensions
IF=6;           % for setting subplots
idxMat=zeros(N*IF); 
idxMat(:)=(1:N*N*IF*IF); 
idxMat=idxMat';

    dimlabels={ 
        'Bodily saliency' % How much feels in body question
        'Mental saliency' % How much involves the mind
        'Valence' % emotional valence question
        'Controllability' % how much you can control question (used to be agency)
        'Last experienced' % last time question
        'Touch' % how strongly associate with touch
    };
diagvals=1:(N+1):N^2

map=[0 0 0;0 0 0;0.8510    0.8510    0.8510];
nbins = 7
for n=1:N^2
    % subplot interpolation
    temp=idxMat(1:IF,1:IF);
    temp=temp(:);
    dx=mod(IF*(n-1),IF*N);
    dy=floor((n-1)/N)*IF*N*IF;
    subplot(N*IF,N*IF,temp+dx+dy) 
    [dim1 dim2]=ind2sub([N N],n);
    % if diagonal then plot histogram
    if(ismember(n,diagvals));
            
        [hhh hb]=hist(data(:,dim1),nbins); 
        % we can't use bar
        %H=bar(hb,hhh,1,'w'
        % let's use patch
        DB=min(diff(hb))/2;
        for bID=1:nbins    
            H=patch([hb(bID)-DB hb(bID)-DB hb(bID)+DB hb(bID)+DB],[0 hhh(bID) hhh(bID) 0],1);
            set(H,'FaceColor',[1 1 1])
            set(H,'EdgeColor',[ 0 0 0])
            hold on
        end
        % Scaling factor for density estimates
        difference_ = diff(hb);   % how wide are the bins
        scaling_factor = sum(hhh * difference_(1)) / 3;  % height of bins*width, divided by three 
        
        axis([-3 3 0 15])   % histogrammin axis
        hold on
        box on
        
        % plot also ksdensities of each sub cluster
        for sim_clu=1:max(sim.sim_cluster.DBSCAN.class_from_mean_data)
            xi=-3:.1:3;
            [fi]=ksdensity(data(find(sim.sim_cluster.DBSCAN.class_from_mean_data==sim_clu),dim1),xi); % rivit ja kolumnit oikein           
            plot(xi,scaling_factor*fi,'Color',classColors(sim_clu,:),'LineWIdth',1) % 
            %disp(['SCALED INTEGRAL ' num2str(trapz(xi,scaling_factor*fi))])
            if(n==1) % cluster texts and line
                text(-5+20.5,29-DF*(sim_clu),'—','FontSize',20,'Color',classColors(sim_clu,:),'FontWeight','bold')
                text(-3.6+20.5,28.5-DF*(sim_clu),dbscan_cluster_labels{sim_clu},'FontSize',FS)         
            end          
        end
        [fi]=ksdensity(data(find(sim.sim_cluster.DBSCAN.class_from_mean_data<=0),dim1),xi); % between clusters
        plot(xi,scaling_factor*fi,'Color',[0 0 0],'LineWIdth',1) 
        disp(['SCALED INTEGRAL ' num2str(trapz(xi,scaling_factor*fi))])
        sim_clu=3; 
        if(n==1) % between clusters text and line
            text(-5+20.5,29-DF*(sim_clu),'—','FontSize',20,'Color',[0 0 0],'FontWeight','bold')
            text(-3.6+20.5,28.5-DF*(sim_clu),dbscan_cluster_labels{sim_clu},'FontSize',FS)
        end
        
    else % we do a scatter plot   below the SCATTER PLOTS
      
        if(dim1>dim2) % removes upper triangle 
            axis off  %
            continue; % 
                      % 
        end           % 
        % let's try only the bottom triangle
        hold on
        box on
        shapes={
            'o'
            's'
            };
        for s=1:27 
            if(sim.sim_cluster.DBSCAN.class_from_mean_data(s)>0)
                col=classColors(sim.sim_cluster.DBSCAN.class_from_mean_data(s),:);
            else
                % out-of-clusters are in grey
                col=.3+[0 0 0];
            end
            markershape=shapes{1};
 
            pppp=plot(data(s,dim1),data(s,dim2),markershape,'MarkerEdgeColor',col,'Markersize',5);
        end
        % same axis for everybody
        axis([-1 1 -1 1]*3)        
        
        x=data(:,dim1);      % mean_dataZ
        y=data(:,dim2);      
        Fit = polyfit(x,y,1); % x = x data, y = y data, 1 = order of the polynomial.                  
        ttt=sort(x);                                                               
        yhat=polyval(Fit,ttt);                                                     
        plot([ttt(1) ttt(end)],[yhat(1) yhat(end)],'LineWidth',2,'Color',[0 0 0]); 
        [rsp ppp]=corr(x,y,'type','spearman'); % Spearmanin correlaatio                         
        star='';
        if(ppp<0.05)     
            star='*';    
        end              
        if(ppp<0.005)    
            star='**';   
        end              
        
        if(dim1<dim2)
        th=text(-2.7,2.5,[strrep(num2str(rsp,2),'-','–') star],'Color',[0 0 0],'FontSize',FS); 
        if(length(star)>0)              % bold significant
            set(th,'FontWeight','bold') 
        end                              
        end
        
    end
    
    % fix axes
    set(gca,'FontSize',FS)
    if(dx ~= 0)
        set(gca,'YTick',[0])
        set(gca,'YTickLabel',{' '})
    end
    if(dx == 0 && dy > 0)
        set(gca,'YTick',[-2:2:2])
        set(gca,'YTickLabel',{'–2' '0' '2'})
    end
    
    if(dy ~=  N*IF*IF*(N-1))
        set(gca,'XTick',[-3])
        set(gca,'XTickLabel',{' '})
    else
        set(gca,'XTick',[-2:2:2])
        set(gca,'XTickLabel',{'–2' '0' '2'})
    end
    % fix axis labels
    if(dx == 0)
        ylabel(strrep(dimlabels{dim2},' ',newline));
    end
    
    if(dy==N*IF*IF*(N-1))
        xlabel(strrep(dimlabels{dim1},' ',newline))
    end
    axis square
end

% prepare for printing
set(gcf,'Units', 'Pixels', 'OuterPosition', [0 0 850 1000]);

set(gcf,'color',[1 1 1]);

exportgraphics(gcf,'figs/main/mind_corrplot_subset_matlab.pdf','ContentType','vector')
export_fig 'figs/main/mind_corrplot_subset_matlab.png';

%% run a PCA and a dendrogram based on the means across subjects 3D plot
% and varimax rotation
pcadata=[mean_data];
[coeff, score, latent, tsquared, explained] = pca(pcadata);

[L1, T]=rotatefactors(coeff(:,1:3));
score_rot=pcadata*L1;
figure(1001)
subplot(1,2,1)
imagesc(L1)
colorbar
title('PCA loadings')
xlabel('PCs');
set(gca,'YTick',[1:6])
set(gca,'YTickLabel',dim_labels)

for s=1:size(score_rot,1)
    figure(1001)
    subplot(1,2,2)
    hold on
    plot3(score_rot(s,1),score_rot(s,2),score_rot(s,3),'o','MarkerSize',10,'Color',map(classID(s),:),'MarkerFaceColor',map(classID(s),:))
end
hold on
for cc=1:27
    text(score_rot(cc,1),score_rot(cc,2),score_rot(cc,3), labels_en{cc})
    hold on
end

xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
grid on
view(3)
camlight
lighting PHONG
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color',[1 1 1]);


%% MDS plots 
load('output/sim/sim_mds.mat');

close all

mds_out = sim_mds.cmd

mds_out=mds_out./max(abs(mds_out(:)))*380; % scaling
    
IF=2;

% Sim MDS plot
for s=1:27
    figure(2000)
    %subplot(6*IF,6*IF,[1:20*(IF^2)])
    hold on
    plot(mds_out(s,1),mds_out(s,2),'o','MarkerSize',5,'Color',classColors(classID(s),:),'MarkerFaceColor',classColors(classID(s),:))
    al='left';
    th=text(mds_out(s,1)+3,mds_out(s,2)+1,sensations{s},'Rotation',0,'Color',classColors(classID(s),:),'FontSize',11,'HorizontalAlignment',al);
end
set(gca, 'XDir','reverse')      
d = daspect;
box off
axis off
close Figure 2000
    
% plot the dimensions below
temp=round(mds_out-min(mds_out))+40;       
mapmap=flipud(cbrewer('div','RdBu',11));   
mapmap(6,:)=[1 1 1]                        % white
mean_data_centered=zscore(mean_dataZ);     
for d=1:6                                  
    figure(2020+d)                         
    tempMat=zeros(max(temp)+40);           
    for s=1:27                             
        tempMat(temp(s,1), temp(s,2))=mean_data_centered(s,d); 
    end
    tempMatB=imgaussfilt(tempMat,17,'FilterSize',81);          
        
    for s=1:27
        MulFac(s)=tempMat(temp(s,1), temp(s,2))/tempMatB(temp(s,1), temp(s,2)); 
    end
 
    imagesc(median(MulFac)*tempMatB',[-3 3]) 
    axis xy                                  
    
    colormap(mapmap)                         
    %colorbar                                
    xlabel(dim_labels{d},'FontSize',140)
    set(gca,'XTickLabel',[])
    set(gca,'XTick',[]) % xtikit pois
    set(gca,'YTick',[])
    set(gca,'YTickLabel',[])
    set(gca, 'XDir','reverse')             
    set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    set(gcf,'color',[1 1 1]);
    exportgraphics(gcf,['figs/main/gaussian_filtered_mean_values_mind_MDS' dim_labels{d} '.pdf'],'ContentType','vector')
    print(gcf, ['figs/main/gaussian_filtered_mean_values_mind_MDS' dim_labels{d} '.png'], '-dpng', '-r300')
    % Uncomment below if you want to see the range of values for good colormap range
    %disp(['MAX VALUE ' num2str(max(max(testi))) ' MIN VALUE ' num2str(min(min(testi)))])
end

%% Double checking that the stimuli are in right positions (they are)
    for s=1:27                             
        tempMat(temp(s,1), temp(s,2))= 1000;
    end

    figure(2028)
    imagesc(tempMat')
    colormap(map)
    colorbar
    axis xy
    set(gca, 'XDir','reverse')    
