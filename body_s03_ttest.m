%%
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
  
bspm = bodySPM_ttest_love(cfg);

%% transform t-tests into effect sizes
% replace tvalues with effect sizes
for cc=1:27
   bspm.ttest.es(:,:,cc)= bspm.ttest.tval(:,:,cc)/sqrt(bspm.ttest.df(cc)+1);
end
save([cfg.outdata 'bspm_ttest'], 'bspm')

%% compute reliability of means 
% prepare data as vectors
base=uint8(imread('bodySPM_base2.png'));
mask=uint8(imread('bodySPM_base3.png'));
mask=mask*.85;
base2=base(10:531,33:203,:);
in_mask=find(mask>128);
subjects=textread(cfg.list,'%s');
Nsubj=size(subjects,1);


for ns=1:Nsubj;
    disp(['Processing subject ' subjects{ns} ' which is number ' num2str(ns) ' out of ' num2str(Nsubj)]);
    matname=[cfg.datapath '/' subjects{ns} '.mat'];
    load (matname) % 'resmat','resmat2','times'
    if(cfg.doaverage==0)
        tempdata=reshape(resmat,[],size(resmat,3));
    else 
        tempdata=reshape(resmat2,[],size(resmat2,3));
    end
	if(ns==1)
        [ length(in_mask) Nsubj size(tempdata,2) ]
		alldata=zeros(length(in_mask),Nsubj,size(tempdata,2));
	end
	alldata(:,ns,:)=tempdata(in_mask,:);
end

% for pixel 1 we run it to get the permu
cfg=[]
cfg.data=(squeeze(alldata(1,:,:)))'; % all answers to pixel 1 across stimuli and participants
cfg.niter=100;
[icc permu iccp]=bramila_splitsample(cfg);
cfg.permu=permu;

% for each pixel we compute a median value of reliability
for pixel=1:size(alldata,1); % for every pixel in mask, 50291 pixels
    cfgtemp=cfg;
    disp([num2str(pixel) ' out of ' num2str(size(alldata,1))])
    cfgtemp.data=(squeeze(alldata(pixel,:,:)))';
    [icc permu iccp]=bramila_splitsample(cfgtemp); % Spearman correlation (split samples) for amount of paint in the pixel in different stimuli
    pvm(pixel,1)=median(icc); % median split sample from 100 permutations
    pvmp(pixel,1)=median(iccp);% median p-values
end

save output/body/split_sample.mat pvm pvmp permu mask

%% First plot with original Finnish stimuli

load_labels
cfg=[]
cfg.bspm = bspm
cfg.type = 'ttest'
cfg.labels = labels_en
cfg.NumCol = 64
bodySPMLOVE(cfg) 

%% LOCAL FUNCTION FOR PLOTTING THE RESULTS
function bodySPMLOVE(cfg)
% Input parameter contains ttest OR ttest2 OR glm results
% cfg.bspm
% cfg.type = ttest 
% cfg.M = optional max tvalue
mask=uint8(imread('bodySPM_base3.png'));
in_mask=find(mask>128); % list of pixels inside the mask
base=uint8(imread('bodySPM_base2.png'));
base2=base(10:531,33:203,:); % single image base
labels={'Itserakkaus'
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
'Rakkaus ei-ihmiseläimiin';
}

if (strcmp(cfg.type,'ttest'))
	NC=size(cfg.bspm.ttest.tval,3); % number of conditions
	
    if(isfield(cfg,'M')) 
        M=cfg.M;
    else
        M=max(abs(cfg.bspm.ttest.tval(:)));
        M=round(prctile(abs(cfg.bspm.ttest.tval(:)),99.9)/5)*5;
    end
    
    
	NumCol=cfg.NumCol;
	th=cfg.bspm.ttest.tTH(2);
	if(isempty(th)) 
		% using uncorrected T-value threshold
		th=2.9719; 
	end

	non_sig=round(th/M*NumCol); % proportion of non significant colors
	hotmap=hot(NumCol-non_sig);
	coldmap=flipud([hotmap(:,3) hotmap(:,2) hotmap(:,1) ]);
	hotcoldmap=[
		coldmap
		zeros(2*non_sig,3);
		hotmap
		];
	if(0)
	% reshaping the tvalues into images
	tvals_for_plot=zeros(size(mask,1),size(mask,2),NC);
	for condit=1:NC
		temp=zeros(size(mask));
		temp(in_mask)=tdata(:,condit);
		temp(find(~isfinite(temp)))=0; % we set nans and infs to 0 for display
		max(temp(:))
		tvals_for_plot(:,:,condit)=temp;
	end
	end
	% plotting
	plotcols = ceil((NC+1)/2); %set as desired
	plotrows = ceil((NC+1)/plotcols); % number of rows is equal to number of conditions+1 (for the colorbar)
	for n=1:NC
		figure(1100)
		subplot(plotrows,plotcols,n)
		imagesc(base2);
		axis('off');
		set(gcf,'Color',[1 1 1]);
		hold on;
		%over2=tvals_for_ploat(:,:,n);
		over2=cfg.bspm.ttest.tval(:,:,n);
		fh=imagesc(over2,[-M,M]);
		axis('off');
		axis equal
		colormap(hotcoldmap);
		set(fh,'AlphaData',mask)
		if floor(n/2)==n/2 % (MIKKE) tämä koodi lisätty erottelemaan otsikoiden rivit, muuten päällekäin
           title(labels(n),'FontSize',8, 'Position',[85 -40 0])
        else 
           title(labels(n),'FontSize',8, 'Position', [85 -120 0])
        end
		if(n==NC)
			subplot(plotrows,plotcols,n+1)
			fh=imagesc([-M-eps,M+eps]);
			axis('off');
            clrbr=colorbar;
            set(clrbr, 'Limits', [0 10]);
            ttx=xlabel(clrbr,'Colorbar')
            set(ttx, 'FontWeight', 'bold');
            
		end
	end
else
	error('not yet implemented')
end

end


