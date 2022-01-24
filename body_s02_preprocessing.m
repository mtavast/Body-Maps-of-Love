clear all
close all

addpath('./external/bodyspm')
cfg.outdata = './output/body/';
cfg.datapath = '/m/nbe/scratch/love/Love final data/Rakkaus_keho/subjects';

cfg.Nstimuli = 27;
cfg.Nempty = 1;
cfg.list = [cfg.outdata 'body_list.txt'];


cfg.doaverage = 0;
cfg.Nbatches=1
cfg.posneg=1;                           % note
cfg.overwrite=1;
bspm = bodySPM_preprocess_batch(cfg);
save([cfg.outdata '/bspm.mat'])

%% quality control and filtering
if(0)
for ss=1:cfg.Nstimuli
	figure(ss)
	loglog(squeeze(bspm.allTimes(ss,1,:)),squeeze(bspm.allTimes(ss,2,:)),'.'); 
	axis square
	saveas(gcf,[cfg.outdata '/allTimes_' num2str(ss)  '.png'])
end
end

% find amount of stimuli painted by each subjects
npainted=sum(~isnan(squeeze(bspm.allTimes(:,2,:)))); % preprocessing stores the amount of pixels painted (allTimes(:,2)), if NaN, then no paint
figure(1)
subplot(2,2,1)
stem(npainted);
subplot(2,2,2)
stem(sort(npainted));


ids=find(bspm.tocheck<=14);
ids=find(npainted>=14); % All that have painted over half of the bodies
subjects=textread(cfg.list,'%s');
out=[];
for i = 1:length(ids)
    out{i,1}=subjects{ids(i)};
end


fileID=fopen([cfg.outdata '/whitelist.txt'],'w');
for i=1:length(out);
    if(~strcmp(out{i}(1),'F'))
        disp([out{i} ' has invalid ID for this study']);
        continue;
    end
    fprintf(fileID,'%s\n',[out{i}]);
end
fclose(fileID);

%% Checking outliers 
% for the total amount of paint, just to see if there are any exremely out
% of the sample participants

% Check the paint-intensities
allTimes=bspm.allTimes;                       % amount of painted stored here
histogram(allTimes(:,2,:))
amount_of_paint = allTimes(:,2,:);            
amount_of_paint = squeeze(amount_of_paint);   % remove the extra dimension

% zscore across all values because there are theoretical reasons why some stimuli could have "real" outliers
% e.g. for a religious person love for god could elicits a lot stronger sensations than for a majority of our participants
paint = reshape(amount_of_paint, [length(amount_of_paint)*height(amount_of_paint) 1]);        % reshape because zscore function does not work with NaN's
zpaint = []
zpaint = (amount_of_paint - mean(paint, 'omitnan')) ./ std(paint, 'omitnan')

% maybe percentiles would be better?

figure(2); histogram(zpaint)                                    % Check the histogram of Z-scores
[maxi indx] = max(zpaint, [], 'all', 'linear'); 
maxi                                                            % there is one clear outlier with a z-score above 16
amount_of_paint(indx)                                           % the outlier is subject 29 stimuli 18
subjects(29)                                                    % subject 29 is FI312418


%% remove the outlier participant from whitelist

out2=[];
ids2=[];
for i = 1:length(out)
    if ~isequal(out{i,1},'FI312418');    
        out2= [out2; cellstr(out{i})]
        ids2= [ids2; ids(i)]                      
    end
end

out = out2                               % for the behavioral data below
ids = ids2'

fileID=fopen([cfg.outdata '/whitelist.txt'],'w');
for i=1:length(out2);
    fprintf(fileID,'%s\n',[out2{i}]);
end
fclose(fileID);

%% behav report

bdata=load([cfg.outdata '/behav_data.mat']);
behav=bdata.behav_data;
behav_labels=bdata.behav_labels;

if(length(out)~=length(bdata.out))
    behav=behav(ids,:); % only final subjects 
end

disp( '---------------------------------')
disp( '>>> Summary statistics report <<<')
disp( '---------------------------------')
disp(['    Final subjects: ' num2str(size(behav,1))])
disp(['Males/females/other or dont want to define: ' num2str(length(find(behav(:,1)==0))) '/' num2str(length(find(behav(:,1)==1))) '/' num2str(length(find(behav(:,1)==2)))])
disp(['18-25: ' num2str(length(find(behav(:,2)==1)))])
disp(['26-35: ' num2str(length(find(behav(:,2)==2)))])
disp(['36-45: ' num2str(length(find(behav(:,2)==3)))])
disp(['46-55: ' num2str(length(find(behav(:,2)==4)))])
disp(['56-65: ' num2str(length(find(behav(:,2)==5)))])
disp(['66 or more: ' num2str(length(find(behav(:,2)==6)))])
disp(['# of parents: ' num2str(length(find(behav(:,6)==1)))])
disp(['perus/keski/korkea: ' num2str(length(find(behav(:,4)==0))) '/' num2str(length(find(behav(:,4)==1))) '/' num2str(length(find(behav(:,4)==2)))])
disp(['Atheist/agnostic/religious: ' num2str(length(find(behav(:,7)==0))) '/' num2str(length(find(behav(:,7)==1))) '/' num2str(length(find(behav(:,7)==2)))])
histogram(behav(:,2))
behav_data = behav;                     
save([cfg.outdata '/behav_data_whitelist.mat'],'behav_data','behav_labels','out','allTimes')

%% Individual subject body map quality check

lista = ['./output/body/body_list.txt'];
bodylist = textread(lista,'%s')
bodylist(131) = cellstr('FI990065') % adding the subject for code reproducibility, the subject is removed from the final data based on the data quality here (removed in script body_s01)
f = {'name'}
bodylist = cell2struct(bodylist,f,130)

subjects = bodylist


% the base image used for painting (in our case only one sided since we
% subtract values)
base=uint8(imread('baseindividual.png'));
base2=base(10:531,33:203,:); % single image base
labels={'Itserakkaus' % original Finnish stimuli
'Vanhempain rakkaus'
'Romanttinen rakkaus'
'Rakkaus ystäviin'
'Rakkaus lähimmäisiin'
'Kauneuden rakkaus'
'Rakkaus kotimaahan'
'Rakkaus jumalaan'
'Luonnonrakkaus'
'Universaali rakkaus'
'Seksuaalinen rakkaus'
'Moraalinen rakkaus'
'Rakkaus sisaruksiin'
'Käytännöllinen rakkaus'
'Viisauden rakkaus'
'Äidin rakkaus lapseen'
'Isän rakkaus lapseen'
'Tosirakkaus'
'Vastavuoroinen rakkaus'
'Ehdoton rakkaus'
'Kumppanuusrakkaus'
'Intohimoinen rakkaus'
'Altr. (epäitsekäs) rakkaus'
'Hyväntahtoinen rakkaus'
'Rakkaus elämään'
'Rakkaus muukalaisiin'
'Rakkaus ei-ihmiseläimiin'};

mask=imread('maskindividual.png');


for s=1:length(subjects)
    % skip dot and dotdot folders
    if(strcmp(subjects(s).name(1),'.')) continue; end 

    % Data loading
    % let's load the subject's answers into a variable a
    data=load_subj([cfg.datapath '/' subjects(s).name],2);
    data(1) = []
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
        h=fspecial('gaussian',[15 15],5);
        over=imfilter(over,h);
        % we subtract left part minus right part of painted area
        % values are hard-coded to our web layout
        over2=over(10:531,33:203,:)-over(10:531,696:866,:);
        resmat(:,:,n)=over2;
    end
  
    % visualize subject's data
        
    M=max(abs(resmat(:))); % max range for colorbar
    NumCol=64;
    hotmap=hot(NumCol);
    coldmap=flipud([hotmap(:,3) hotmap(:,2) hotmap(:,1) ]);
    hotcoldmap=[
        coldmap
        hotmap
    ];

	plotcols = 15; %set as desired
    plotrows = ceil((NC+1)/plotcols); % number of rows is equal to number of conditions+1 (for the colorbar)


    for n=1:NC
        figure(s)
        set(gcf, 'Position', get(0, 'Screensize'))
        subplot(plotrows,plotcols,n)
        imagesc(base2);
        axis('off');
        set(gcf,'Color',[1 1 1]);
        hold on;
        over2=resmat(:,:,n);
        fh=imagesc(over2,[-M,M]);
        axis('off');
        axis equal
        colormap(hotcoldmap);
        set(fh,'AlphaData',mask)
        if floor(n/2)==n/2 
           title(labels(n),'FontSize',8, 'Position',[85 -40 0])
        else 
           title(labels(n),'FontSize',8, 'Position', [85 -120 0])
        end
        if(n==NC) 
            subplot(plotrows,plotcols,n+1)
            fh=imagesc(ones(size(base2)),[-M,M]);
            axis('off');
            colorbar;
            fname = './output/individualbodies'
            saveas(gcf, fullfile(fname, [subjects(s).name '.png']), 'png');
            close(gcf)
        end
    end
end

