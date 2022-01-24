clear all
close all
addpath('/m/nbe/scratch/braindata/shared/toolboxes/cbrewer/');
addpath(genpath('external'));
load('output/mind/M10_data')
load_labels
sensations=labels_en;

%% Loop to get only the asnwers where lapse > 950

never = data(:,5,:) > 950; 

% add as a column
data(:,7,:) = never; 

% replace with nans, if never experienced
for s=1:length(data)
    for k=1:27
        if data(k,7,s) == 1;
           data(k,:,s) = NaN;
        end
    end
end

data(:,7,:) = [];  

%% Percentage of participants who reported over 950 in the lapse dimension
nsubjects = length(good_subj_list)

never2 = zeros(27,1,nsubjects)

for i=1:27
    for s=1:nsubjects
        if never(i,1,s) == 1
            never2(i,1,s) = 0;
        else
            never2(i,1,s) = 1;
        end
    end
end

for i=1:27
    experienced(i,1) = sum(never2(i,:,:)); 
end

not_experienced_percentage = ((nsubjects - experienced) / nsubjects) * 100;

%% data = words x dimensions x good subjects
% 6 dimensions listed here:
Ndim=6;
dimlabels={ 
    'Bodily saliency' % How much feels in body question
    'Mental saliency' % How much involves the mind
    'Valence' % emotional valence question
    'Controllability' % how much you can control question (used to be agency)
    'Last experienced' % last time question
    'Touch' % how strongly associate with touch
};

sensations_classes 
map=cbrewer('qual','Set1',max(classID));
map(6,:)=[0 0 0];
data(find(data==0))=1; % just to avoid problems in binning
data(find(data==1000))=999; % converts 1000 to 999 to avoid inflating the values after the atanh function 

mean_data=zeros(27,6);
median_data=zeros(27,6);
sem_data=zeros(27,6);
mean_dataZ=zeros(27,6);
median_dataZ=zeros(27,6);
sem_dataZ=zeros(27,6);

% for plotting distros
Nres=100;
distros=zeros(27,Nres+1,6);
distrosZ=distros;
maph=flipud(cbrewer('div','Spectral',31));


% compute means, medians and distributions
for s=1:length(sensations)
    thisdata=squeeze(data(s,:,:));
    thisdata=(thisdata-500)/500;
    ids=find(~isnan(max(thisdata)));
    thisdata=thisdata(:,ids); % let's keep only those with data
    disp(['LENGTH OF STIMULUS ' num2str(s) ' IS ' num2str(length(thisdata))])
    xi=((0:Nres)-Nres/2)/(Nres/2);
    xiZ=5*(((0:Nres)-(Nres/2))/(Nres/2));
    for subdim=1:6
        [fi]=ksdensity(thisdata(subdim,:)',xi);
        distros(s,:,subdim)=fi/sum(fi);
    end
    thisdataZ=atanh(thisdata-eps);
    for subdim=1:6
        [fi]=ksdensity(thisdataZ(subdim,:)',xiZ);
        distrosZ(s,:,subdim)=fi/sum(fi);
    end
    mean_data(s,:)=nanmean(thisdata,2);
    mean_dataZ(s,:)=nanmean(thisdataZ,2);
    median_data(s,:)=nanmedian(thisdata,2);
    median_dataZ(s,:)=nanmedian(thisdataZ,2);
    Nitems(s,1)=size(thisdata,2);
    std_data(s,:)=std(thisdata,0,2);
    sem_data(s,:)=std(thisdata,0,2)./sqrt(size(thisdata,2));
    sem_dataZ(s,:)=std(thisdataZ,0,2)./sqrt(size(thisdataZ,2));
    sens_count(s,1)=size(thisdata,2);
end

%% Correlation matrix
figure(123412)
corrplot(mean_dataZ,'varnames',dim_labels,'type','Spearman','testR','on')

%% label sorting with dim 3 (emotions) for mean_data and mean_dataZ
% new labels with n
sensations_n = {}
for i=1:length(sensations)
    temp = sensations{i};  % label
    temp2 = append(temp, ' (n=', num2str(experienced(i))) % how many over 950 in last experienced
    temp3 = append(temp2, ')')
    sensations_n{i} = temp3 
end

[temp emomedian]=sort(median_data(:,3));
maph=flipud(cbrewer('div','Spectral',11)); 

emomean_mean_data=mean_data(emomedian,:);
emomean_median_data=median_data(emomedian,:);
emomean_distros=distros(emomedian,:,:);
emomean_sensations=sensations_n(emomedian);
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
[temp emomedianZ]=sort(median_dataZ(:,3));

emomean_mean_dataZ=mean_dataZ(emomedianZ,:);
emomean_median_dataZ=median_dataZ(emomedianZ,:);
emomean_distrosZ=distrosZ(emomedianZ,:,:);
emomean_sensationsZ=sensations_n(emomedianZ);

[sensations(emomedian) sensations(emomedianZ)]

load output/sim/sim_cluster.mat
figure(22001)
for subdim=1:6
    subplot_= subplot(1,6,subdim)
    subplot_.Position = subplot_.Position + [0.032 0 0 0] % push little bit on the rightif   
    for s=27:-1:1
        odd_or_even=mod(s,2);
        p=patch([xiZ'; flipud(xiZ')],[s+(squeeze(30*emomean_distrosZ(s,:,subdim)))'; s*ones(length(xiZ),1)] ,1);
        colID=sim_cluster.DBSCAN.class_from_mean_data(emomedianZ(s));
        if(colID<=0)
            col=classColors(end,:);
        else
            col=classColors(colID,:);
        end
        %set(p,'EdgeAlpha',0)
        set(p,'EdgeColor',[1 1 1]*.9)
        set(p,'FaceColor',col);
        hold on
        switch(4)
            case 1 % labels on both side alternating
                if(subdim==1)
                    th=text(-1*(4.1+odd_or_even*7),s+.25,emomean_sensationsZ{s},'HorizontalAlignment','right','FontSize',12)
                    if(odd_or_even) text(-1*(4+odd_or_even*7),s+.9,'___________','Color',[1 1 1]*.88,'Interpreter','none','HorizontalAlignment','left');end
                end
                if(subdim==5)
                    text(4.1+odd_or_even*7,s+.25,emomean_sensationsZ{s},'HorizontalAlignment','left','FontSize',12)
                    if(odd_or_even) text((4),s+.9,'__________','Color',[1 1 1]*.88,'Interpreter','none');end
                end
            case 2 % only on the left alternatinv
                if(subdim==1)
                    th=text(-1*(4.1+odd_or_even*7),s+.25,emomean_sensationsZ{s},'HorizontalAlignment','right','FontSize',12)
                    if(odd_or_even) text(-1*(4+odd_or_even*7),s+.9,'___________','Color',[1 1 1]*.88,'Interpreter','none','HorizontalAlignment','left');end
                end
            case 3 % alternating between dimensions
                if(subdim==1 && odd_or_even==0)
                    th=text(-1*(4.1),s+.25,emomean_sensationsZ{s},'HorizontalAlignment','right','FontSize',12)
                end
                if(subdim==5 && odd_or_even==1)
                    th=text(4.1,s+.25,emomean_sensationsZ{s},'HorizontalAlignment','left','FontSize',12)
                end
            case 4 % labels everywhere without alternatinv
                if(subdim==1)
                    th=text(-1*(4.1),s+.25,emomean_sensationsZ{s},'HorizontalAlignment','right','FontSize',12)
                end
                %if(subdim==6)
                %    th=text(4.1,s+.25,emomean_sensationsZ{s},'HorizontalAlignment','left','FontSize',12)
                %end
                
                
        end

    end
    set(gca,'YTick',[])
    xlabel(dim_labels{subdim})
    axis([-4 4 0 30])
    axis off

    text(0,-.4,'0','HorizontalAlignment','center','FontSize',12)
    text(4,-.4,'4','HorizontalAlignment','center','FontSize',12)
    text(-4,-.4,'-4','HorizontalAlignment','center','FontSize',12)
    text(0,-2.6,dim_labels{subdim},'HorizontalAlignment','center', 'FontSize', 12)
end

set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color',[1 1 1]);
sgtitle('After excluding data with Last experienced > 950', 'FontWeight','bold', 'FontSize',20)

print(gcf, 'figs/supplement/distributions_never_experienced.png', '-dpng', '-r500')



%% myCorr plot 
%close all
sim=load('output/sim/sim_cluster.mat')
FS=13.5; % font size
DF=7; % distance between labels

dbscan_cluster_labels={
    'Interpersonal love'
    'Love for ideas and non-human animals'
    'Between clusters'
    };
%sensations_classes

figure(1)
data=mean_dataZ;
N=size(data,2); 
IF=6;           
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
            
        [hhh hb]=hist(data(:,dim1),nbins); % dim1 takes correct dimension checked
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
            if(n==1) % cluster texts and line
                text(-5+20.5,22-DF*(sim_clu),'—','FontSize',20,'Color',classColors(sim_clu,:),'FontWeight','bold')
                text(-3.6+20.5,21.5-DF*(sim_clu),dbscan_cluster_labels{sim_clu},'FontSize',FS)      
            end          
        end
        [fi]=ksdensity(data(find(sim.sim_cluster.DBSCAN.class_from_mean_data<=0),dim1),xi); % between clusters
        plot(xi,scaling_factor*fi,'Color',[0 0 0],'LineWIdth',1) 
        sim_clu=3; 
        if(n==1) % between clusters text and line
            text(-5+20.5,22-DF*(sim_clu),'—','FontSize',20,'Color',[0 0 0],'FontWeight','bold')
            text(-3.6+20.5,21.5-DF*(sim_clu),dbscan_cluster_labels{sim_clu},'FontSize',FS)
        end
        
    else % we do a scatter plot   
      
        if(dim1>dim2) % remove upper triangle
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
        y=data(:,dim2);      % 
        Fit = polyfit(x,y,1); % x = x data, y = y data, 1 = order of the polynomial.                   
        ttt=sort(x);                                                               
        yhat=polyval(Fit,ttt);                                                     
        plot([ttt(1) ttt(end)],[yhat(1) yhat(end)],'LineWidth',2,'Color',[0 0 0]);  
        [rsp ppp]=corr(x,y,'type','spearman');                        
        star='';
        if(ppp<0.05)     
            star='*';    
        end              
        if(ppp<0.005)    
            star='**';   
        end              
        
        if(dim1<dim2)
        th=text(-2.7,2.5,[strrep(num2str(rsp,2),'-','–') star],'Color',[0 0 0],'FontSize',FS); 
        if(length(star)>0)              
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


sgtitle('After excluding data with Last experienced > 950', 'FontWeight','bold', 'FontSize',20)
% prepare for printing
set(gcf,'Units', 'Pixels', 'OuterPosition', [0 0 1000 1200]);

set(gcf,'color',[1 1 1]);

print(gcf, 'figs/supplement/corrplot_never_experienced.png', '-dpng', '-r500')

