%% Effect size plots
% First, arrange in amount of significant t-values
clear all; close all;
addpath(genpath('external/'))
load('output/body/bspm_ttest.mat');

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
bspmsort.ttest.tval = bspmsort.ttest.tval(:,:,idx); % sort the tvalues
bspmsort.ttest.es = bspmsort.ttest.es(:,:,idx);     % and effect sizes
labelssort = tvalues.labels(idx)


%% Plot the effect sizes

cfg.bspm = bspmsort % use the sorted version
mask=uint8(imread('bodySPM_base3.png'));
in_mask=find(mask>128); % list of pixels inside the mask
base=uint8(imread('bodySPM_base2.png'));
base2=base(10:531,33:203,:); % single image base
labels = labelssort


	NC=size(cfg.bspm.ttest.es,3); % number of conditions
	M=max(abs(cfg.bspm.ttest.es(:)))
	M=round(prctile(abs(cfg.bspm.ttest.es(:)),99.9)/.1)*.1 
	NumCol=64
    th=cfg.bspm.ttest.tTH(2)/sqrt(min(cfg.bspm.ttest.df+1)); % threshold for effect sizes

	non_sig=round(th/M*NumCol); % proportion of non significant colors
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
	H=size(bspm.ttest.es,1); % height of picture
    W=size(bspm.ttest.es,2); % width of picture
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
		over2=cfg.bspm.ttest.es(:,:,n);
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
       ttx=xlabel(clrbr,'Effect Size')
       set(ttx, 'FontWeight', 'bold', 'FontSize', 16);  
       
       h=gcf;

set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

%exportgraphics(gcf,'figs/main/sensations_body.pdf','ContentType','vector')
%print(gcf, 'figs/main/sensations_body.png', '-dpng', '-r500')
print(gcf, 'figs/main/sensations_body.jpg', '-djpeg', '-r400') % journal format
%% Set up for descriptive plots
% Loads each subjects data directly

% The t-test results should be considered as a visualization method that
% is affected by the variance and intensity of the paint values. However,
% testing againts zero is somewhat arbitrary, as the paint values can't be
% negative. Also, the distributional assumptions of the t-test complicate
% the interpretation. For these reasons, the data is also plotted as below.
% This descriptive figure visualizes how many participants painted at least something
% to each body location, thus providing an alternative and easy-to-interpret view of the
% group level data. 

clear all
close all

cfg.datapath = '/m/nbe/scratch/love/Love final data/Rakkaus_keho/subjects';
lista = ['./output/body/whitelist.txt'];
whitelist = textread(lista,'%s');
f = {'name'};
NSu = length(whitelist)
whitelist = cell2struct(whitelist,f,NSu);


% only final subjects
subjects = whitelist;

% the base image used for painting (in our case only one sided since we
% subtract values)
base=uint8(imread('baseindividual.png'));
base2=base(10:531,33:203,:); % single image base
mask=imread('maskindividual.png');


for s=1:length(subjects)
    % skip dot and dotdot folders
    if(strcmp(subjects(s).name(1),'.')) continue; end 

    % Data loading
    % let's load the subject's answers into a variable a
    data=load_subj([cfg.datapath '/' subjects(s).name],2);
    data(1) = [];
    NC=length(data); % number of conditions
    
    % Painting reconstruction
    % 'data' now contains all mouse movements. What we need are the mouse
    % locations while the button was pressed (i.e. during painting)
    % Furthermore, the painting tool has a brush size. We recreate that
    % using image filter
    
    for n=1:NC;
        T=length(data(n).paint(:,2)); % number of mouse locations
        over=zeros(size(base,1),size(base,2)); % empty matrix to reconstruct painting
        for t=1:T
            y=ceil(data(n).paint(t,3)+1);
            x=ceil(data(n).paint(t,2)+1);
            if(x<=0) x=1; end
            if(y<=0) y=1; end
            if(x>=900) x=900; end % hardcoded for our experiment, you need to change it if you changed layout
            if(y>=600) y=600; end % hardcoded for our experiment, you need to change it if you changed layout
            over(y,x)=over(y,x)+1;
        end
        % Simulate brush size with a gaussian disk
        %h1=fspecial('gaussian',[15 15],5); 
        % Note: The commented line above is the standard way to create the brush size of the online tool.
        % It recreates the brush size, and smooths the values.
        % However, for this analysis it is slighly problematic, as the shape of this filter is actually square. 
        % In the t-test this should have little effect, as the values towards the edges get miniscule values.
        % However, here I will transform every pixel that has ANY VALUE AT
        % ALL, to have a value of 1, it is better to use an actual circular
        % filter. 
        h=fspecial('disk', 7)    % A square filter, same diameter as the length of the square in gaussian [15 15]
        over=imfilter(over,h);
        % we subtract left part minus right part of painted area
        % values are hard-coded to our web layout
        %over2=over(10:531,33:203,:)-over(10:531,696:866,:);
        % no substraction
        over2=over(10:531,33:203,:);
        resmat(:,:,n)=over2;
    end
  
resmats.data{s} = resmat;
resmats.subject(s)= s;

end


%% Checking that there are no NaN's in the data

for i=1:length(subjects)
    numberofnan(i,1) = sum(isnan(resmats.data{1,i}),'all');
end

%% Number of empty bodies
mask=imread('maskindividual.png');
in_mask=find(mask>128);


for i=1:length(subjects)
    for s=1:27
    temp = resmats.data{1,i}(:,:,s);
    temp2 = temp(in_mask);
    empties(i,s) = sum(temp2,'all');
    end
end

for s=1:27
    for i=1:length(subjects)
        if empties(i,s) == 0  % if the sum is 0, the body is empty
        empties2(i,s) = 0;
        else
        empties2(i,s) = 1;
        end
    end
end

percentage_nonempty = mean(empties2)' * 100 % agrees with the method in body_s03_ttest

%% Loop to calculate the percentage of participants who painted in each body area

NSti = 27
sumvalues = zeros(522,171,NSti);
for n=1:NSti;                           % every stimulus
    for s=1:length(subjects);                        % every subject
    temp = resmats.data(s);             % paint data for each subject and stimulus
    temp = temp{1,1};                   % same as above
        for a = 1:522;                  % each row
            for b = 1:171;              % each column
                if temp(a,b,n) ~= 0;    % if has a value, set to 1 (otherwise 0)
                   temp(a,b,n) = 1;
                end
            end
        end
    sumvalues(:,:,n) = sumvalues(:,:,n) + temp(:,:,n);     % add all values
    end
end

percentvalues = (sumvalues./NSu) .* 100;                             % convert to percentages

% The same order as in the effect size plot
addpath(genpath('external/'))
load('output/body/bspm_ttest.mat');

bspmsort=bspm; % will be used for sorting the stimuli
load_labels
tvalues= []
tvalues.labels = labels_en;
labels = labels_en;

% loop to find sum of significant t-values in each stimulus 
for i=1:27
       temp_id = find(bspm.ttest.tval(:,:,i) > bspm.ttest.tTH(2));
       length(temp_id)
       temp_body = bspm.ttest.tval(:,:,i);
       tvalues.sums_secondmethod(i,:)=sum(temp_body(temp_id));
end


tvalues = struct2table(tvalues);
[~, idx] = sortrows(tvalues, 2, 'descend');

percentvalues=percentvalues(:,:,idx); % reorder by total sum

labels= labels(idx)                   % reorder also the labels

labelsP = num2cell(percentage_nonempty)

labelsP = labelsP(idx)                % and the percentage labels

%% Plotting

M = 100  % max range for the colormap, set to 100 %
NumCol=64;
th=0*M; % No threshold
non_sig=round(th/M*NumCol); % proportion of non significant colors
hotmap=hot(NumCol-non_sig);
%hotmap=hot(NumCol);
coldmap=flipud([hotmap(:,3) hotmap(:,2) hotmap(:,1) ]);
hotcoldmap=[
 coldmap
 zeros(2*non_sig,3)
 hotmap
];

    plotcols = 9; 
	plotrows = 3; 
	H=size(percentvalues,1);
    W=size(percentvalues,2);
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
		over2=percentvalues(:,:,n);
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
                
                % Percentage labels
                txt = num2str(round(labelsP{n},1))
                txt = append(txt, '%') % round to one decimal
                th2=text(cc*W-W/3.9-1*5,5*-1+rr*H-H/50,txt) % percentages
                set(th2,'HorizontalAlignment','left');%
                set(th2,'FontSize',10);%
                set(th2,'FontWeight','bold', 'FontAngle', 'italic');%
            elseif n == 10 
                rr=ceil(n/plotcols);
                cc=mod(n-1,plotcols)+1;
                templabel=sprintf(strrep(labels{n},' ','\n'));

                for xx =-1:1
                    for yy=-1:1
                th=text(cc*W-W/2+yy*5,5*xx+rr*H-H/4-60,templabel);
                if(xx==0 && yy==0) set(th,'Color','white'); else
                    set(th,'Color',[0 0 0]); end
                set(th,'HorizontalAlignment','center');
                set(th,'FontSize',13);
                set(th,'FontWeight','bold');

                    end
                end
                xx=0;
                yy=0;
                th=text(cc*W-W/2+yy*5,5*xx+rr*H-H/4-60,templabel);
                if(xx==0 && yy==0) set(th,'Color','white'); else
                    set(th,'Color',[0 0 0]); end
                set(th,'HorizontalAlignment','center');
                set(th,'FontSize',13);
                set(th,'FontWeight','bold');
                
                % Percentage labels
                txt = num2str(round(labelsP{n},1))
                txt = append(txt, '%') % round to one decimal
                th2=text(cc*W-W/3.9-1*5,5*-1+rr*H-H/50,txt) % percentages
                set(th2,'HorizontalAlignment','left');%
                set(th2,'FontSize',10);%
                set(th2,'FontWeight','bold', 'FontAngle', 'italic');%
             elseif n == 14 
                rr=ceil(n/plotcols);
                cc=mod(n-1,plotcols)+1;
                templabel=sprintf(strrep(labels{n},' ','\n'));

                for xx =-1:1
                    for yy=-1:1
                th=text(cc*W-W/2+yy*5,5*xx+rr*H-H/4-10,templabel);
                if(xx==0 && yy==0) set(th,'Color','white'); else
                    set(th,'Color',[0 0 0]); end
                set(th,'HorizontalAlignment','center');
                set(th,'FontSize',13);
                set(th,'FontWeight','bold');

                    end
                end
                xx=0;
                yy=0;
                th=text(cc*W-W/2+yy*5,5*xx+rr*H-H/4-10,templabel);
                if(xx==0 && yy==0) set(th,'Color','white'); else
                    set(th,'Color',[0 0 0]); end
                set(th,'HorizontalAlignment','center');
                set(th,'FontSize',13);
                set(th,'FontWeight','bold');
                
                % Percentage labels
                txt = num2str(round(labelsP{n},1))
                txt = append(txt, '%') % round to one decimal
                th2=text(cc*W-W/3.9-1*5,5*-1+rr*H-H/50,txt) % percentages
                set(th2,'HorizontalAlignment','left');%
                set(th2,'FontSize',10);%
                set(th2,'FontWeight','bold', 'FontAngle', 'italic');%
           elseif n == 18 |n == 20 | n == 22 | n == 24 | n == 26
                rr=ceil(n/plotcols);
                cc=mod(n-1,plotcols)+1;
                templabel=sprintf(strrep(labels{n},' ','\n'));

                for xx =-1:1
                    for yy=-1:1
                th=text(cc*W-W/2+yy*5,5*xx+rr*H-H/4-10,templabel);
                if(xx==0 && yy==0) set(th,'Color','white'); else
                    set(th,'Color',[0 0 0]); end
                set(th,'HorizontalAlignment','center');
                set(th,'FontSize',13);
                set(th,'FontWeight','bold');

                    end
                end
                xx=0;
                yy=0;
                th=text(cc*W-W/2+yy*5,5*xx+rr*H-H/4-10,templabel);
                if(xx==0 && yy==0) set(th,'Color','white'); else
                    set(th,'Color',[0 0 0]); end
                set(th,'HorizontalAlignment','center');
                set(th,'FontSize',13);
                set(th,'FontWeight','bold');
                
                % Percentage labels
                txt = num2str(round(labelsP{n},1))
                txt = append(txt, '%') % round to one decimal
                th2=text(cc*W-W/3.9-1*5,5*-1+rr*H-H/50,txt) % percentages
                set(th2,'HorizontalAlignment','left');%
                set(th2,'FontSize',10);%
                set(th2,'FontWeight','bold', 'FontAngle', 'italic');%
            elseif mod(n,2) == 1    
                rr=ceil(n/plotcols);
                cc=mod(n-1,plotcols)+1;
                templabel=sprintf(strrep(labels{n},' ','\n'));

                for xx =-1:1
                    for yy=-1:1
                th=text(cc*W-W/2+yy*5,5*xx+rr*H-H/4-65,templabel);
                if(xx==0 && yy==0) set(th,'Color','white'); else
                    set(th,'Color',[0 0 0]); end
                set(th,'HorizontalAlignment','center');
                set(th,'FontSize',13);
                set(th,'FontWeight','bold');

                    end
                end
                xx=0;
                yy=0;
                th=text(cc*W-W/2+yy*5,5*xx+rr*H-H/4-65,templabel);
                if(xx==0 && yy==0) set(th,'Color','white'); else
                    set(th,'Color',[0 0 0]); end
                set(th,'HorizontalAlignment','center');
                set(th,'FontSize',13);
                set(th,'FontWeight','bold');
                
                % Percentage labels
                txt = num2str(round(labelsP{n},1))
                txt = append(txt, '%') % round to one decimal
                th2=text(cc*W-W/3.9-1*5,5*-1+rr*H-H/50,txt) % percentages
                set(th2,'HorizontalAlignment','left');%
                set(th2,'FontSize',10);%
                set(th2,'FontWeight','bold', 'FontAngle', 'italic');%
                
             else
                rr=ceil(n/plotcols);
                cc=mod(n-1,plotcols)+1;
                templabel=sprintf(strrep(labels{n},' ','\n'));

                for xx =-1:1
                    for yy=-1:1
                th=text(cc*W-W/2+yy*5,5*xx+rr*H-H/4+60,templabel);
                if(xx==0 && yy==0) set(th,'Color','white'); else
                    set(th,'Color',[0 0 0]); end
                set(th,'HorizontalAlignment','center');
                set(th,'FontSize',13);
                set(th,'FontWeight','bold');

                    end
                end
                xx=0;
                yy=0;
                th=text(cc*W-W/2+yy*5,5*xx+rr*H-H/4+60,templabel);
                if(xx==0 && yy==0) set(th,'Color','white'); else
                    set(th,'Color',[0 0 0]); end
                set(th,'HorizontalAlignment','center');
                set(th,'FontSize',13);
                set(th,'FontWeight','bold');
                % Percentage labels
                txt = num2str(round(labelsP{n},1))
                txt = append(txt, '%') % round to one decimal
                th2=text(cc*W-W/3.9-1*5,5*-1+rr*H-H/50,txt) % percentages
                set(th2,'HorizontalAlignment','left');%
                set(th2,'FontSize',10);%
                set(th2,'FontWeight','bold', 'FontAngle', 'italic');%
                
        end
        end    
       
       clrbr=colorbar;
       set(clrbr, 'Limits', [0 100]);
       ttx=xlabel(clrbr,'Percentage of participants reporting sensations')
       set(ttx, 'FontWeight', 'bold', 'FontSize', 16);
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

exportgraphics(gcf,'figs/main/sensations_body.pdf','ContentType','vector')
print(gcf, 'figs/supplement/sensations_body_descriptive.png', '-dpng', '-r500')
