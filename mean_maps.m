%% Mean, standard deviaton, confidence intervals bodymaps

clear all
close all

addpath(genpath('external/'))

%% t-test
cfg=[];
cfg.outdata = './output/body/';
cfg.datapath = './output/body/';
cfg.list = [cfg.outdata 'whitelist.txt'];
cfg.doaverage = 0;
cfg.niter=0; % do not run cluster correction
  
bspm = bodySPM_ttest_varmaps(cfg);


%% Sort to the same order as the main plot

bspmsort=bspm; % will be used for sorting the stimuli
load_labels
tvalues= []
tvalues.labels = labels_en

% loop to find sum of significant t-values in each stimulus 
for i=1:27
       temp_id = find(bspm.ttest.tval(:,:,i) > bspm.ttest.tTH(2));
       length(temp_id)
       temp_body = bspm.ttest.tval(:,:,i);
       tvalues.sums_secondmethod(i,:)=sum(temp_body(temp_id));
end

tvalues = struct2table(tvalues);
[~, idx] = sortrows(tvalues, 2, 'descend');

% Sort the labels
labelssort = [];
bspmsort.ttest.ci_l = bspmsort.ttest.ci_l(:,:,idx); % sort the ci_l
bspmsort.ttest.ci_h = bspmsort.ttest.ci_h(:,:,idx);     % sort the ci_h
bspmsort.ttest.mean = bspmsort.ttest.mean(:,:,idx);     % sort the means
bspmsort.ttest.std = bspmsort.ttest.std(:,:,idx);     % sort the stds

labelssort = tvalues.labels(idx)

%% Plot 

plotname = {'95% CI for the population mean (lower limit)',
    '95% CI for the population mean (upper limit)',
    'Mean',
    'Standard deviation'}

for plot = 1:4
    
    % get the right data
    if plot == 1;
        tempdata = bspmsort.ttest.ci_l; % '95% CI for the population mean (lower limit)'
    elseif plot == 2;
        tempdata = bspmsort.ttest.ci_h; % '95% CI for the population mean (upper limit)'
    elseif plot == 3
        tempdata = bspmsort.ttest.mean; % 'Mean'
    elseif plot == 4
        tempdata = bspmsort.ttest.std;  % 'Standard deviation'
    end
            
    cfg.bspm = bspmsort % use the sorted version
    mask=uint8(imread('bodySPM_base3.png'));
    in_mask=find(mask>128); % list of pixels inside the mask
    base=uint8(imread('bodySPM_base2.png'));
    base2=base(10:531,33:203,:); % single image base
    labels = labelssort


        NC=size(tempdata,3); % number of conditions
        M=max(abs(tempdata(:)))
        NumCol=64
        th=0; % threshold for effect sizes

        non_sig=0; % proportion of non significant colors
        hotmap=hot(NumCol-non_sig);
        coldmap=flipud([hotmap(:,3) hotmap(:,2) hotmap(:,1) ]);
        hotcoldmap=[
            coldmap
            zeros(2*non_sig,3);
            hotmap
            ];


        % plotting
        plotcols = 9; %set as desired
        plotrows = 3; % number of rows is equal to number of conditions+1 (for the colorbar)
        H=size(tempdata,1); % height of picture
        W=size(tempdata,2); % width of picture
        allmaps=zeros(plotrows*H,plotcols*W); 
        allbases=repmat(base2,plotrows,plotcols);
        allmask=repmat(mask,plotrows,plotcols);
        figure
        imagesc(allbases);
        axis equal
        allover=zeros(size(allmaps));
        for n=1:NC
            rr=ceil(n/plotcols);
            cc=mod(n-1,plotcols)+1;
            axis('off');
            set(gcf,'Color',[1 1 1]);
            hold on;
            over2=tempdata(:,:,n);
            allover((1:H)+(rr-1)*H,(1:W)+(cc-1)*W)=over2;
        end
        hold on
        fh=imagesc(allover,[-M M]);
        set(fh,'AlphaData',allmask)
        axis('off');
            axis equal
            colormap(hotcoldmap);
            for n=1:NC
                if n == 6 
                    rr=ceil(n/plotcols);
                    cc=mod(n-1,plotcols)+1;
                    templabel=sprintf(strrep(labels{n},' ','\n'));

                    for xx =-1:1
                        for yy=-1:1
                    th=text(cc*W-W/2+yy*5,5*xx+rr*H-H/4,templabel);
                    if(xx==0 && yy==0) set(th,'Color','white'); else
                        set(th,'Color',[0 0 0]); end
                    set(th,'HorizontalAlignment','center');
                    set(th,'FontSize',13);
                    set(th,'FontWeight','bold');

                        end
                    end
                    xx=0;
                    yy=0;
                    th=text(cc*W-W/2+yy*5,5*xx+rr*H-H/4,templabel);
                    if(xx==0 && yy==0) set(th,'Color','white'); else
                        set(th,'Color',[0 0 0]); end
                    set(th,'HorizontalAlignment','center');
                    set(th,'FontSize',13);
                    set(th,'FontWeight','bold');
                elseif n == 10 
                    rr=ceil(n/plotcols);
                    cc=mod(n-1,plotcols)+1;
                    templabel=sprintf(strrep(labels{n},' ','\n'));

                    for xx =-1:1
                        for yy=-1:1
                    th=text(cc*W-W/2+yy*5,5*xx+rr*H-H/4-50,templabel);
                    if(xx==0 && yy==0) set(th,'Color','white'); else
                        set(th,'Color',[0 0 0]); end
                    set(th,'HorizontalAlignment','center');
                    set(th,'FontSize',13);
                    set(th,'FontWeight','bold');

                        end
                    end
                    xx=0;
                    yy=0;
                    th=text(cc*W-W/2+yy*5,5*xx+rr*H-H/4-50,templabel);
                    if(xx==0 && yy==0) set(th,'Color','white'); else
                        set(th,'Color',[0 0 0]); end
                    set(th,'HorizontalAlignment','center');
                    set(th,'FontSize',13);
                    set(th,'FontWeight','bold');
                elseif mod(n,2) == 1    
                    rr=ceil(n/plotcols);
                    cc=mod(n-1,plotcols)+1;
                    templabel=sprintf(strrep(labels{n},' ','\n'));

                    for xx =-1:1
                        for yy=-1:1
                    th=text(cc*W-W/2+yy*5,5*xx+rr*H-H/4-55,templabel);
                    if(xx==0 && yy==0) set(th,'Color','white'); else
                        set(th,'Color',[0 0 0]); end
                    set(th,'HorizontalAlignment','center');
                    set(th,'FontSize',13);
                    set(th,'FontWeight','bold');

                        end
                    end
                    xx=0;
                    yy=0;
                    th=text(cc*W-W/2+yy*5,5*xx+rr*H-H/4-55,templabel);
                    if(xx==0 && yy==0) set(th,'Color','white'); else
                        set(th,'Color',[0 0 0]); end
                    set(th,'HorizontalAlignment','center');
                    set(th,'FontSize',13);
                    set(th,'FontWeight','bold');

                 else
                    rr=ceil(n/plotcols);
                    cc=mod(n-1,plotcols)+1;
                    templabel=sprintf(strrep(labels{n},' ','\n'));

                    for xx =-1:1
                        for yy=-1:1
                    th=text(cc*W-W/2+yy*5,5*xx+rr*H-H/4+50,templabel);
                    if(xx==0 && yy==0) set(th,'Color','white'); else
                        set(th,'Color',[0 0 0]); end
                    set(th,'HorizontalAlignment','center');
                    set(th,'FontSize',13);
                    set(th,'FontWeight','bold');

                        end
                    end
                    xx=0;
                    yy=0;
                    th=text(cc*W-W/2+yy*5,5*xx+rr*H-H/4+50,templabel);
                    if(xx==0 && yy==0) set(th,'Color','white'); else
                        set(th,'Color',[0 0 0]); end
                    set(th,'HorizontalAlignment','center');
                    set(th,'FontSize',13);
                    set(th,'FontWeight','bold');

            end
            end       

           clrbr=colorbar;
           set(clrbr, 'Limits', [0 M]);
           ttx=xlabel(clrbr, [char(10) plotname{plot}])
           set(ttx, 'FontWeight', 'bold', 'FontSize', 16);  

           h=gcf;

    set(h,'PaperOrientation','landscape');
    set(h,'PaperUnits','normalized');
    set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    print(gcf,['figs/supplement/' plotname{plot} '.png'], '-dpng', '-r500')
end