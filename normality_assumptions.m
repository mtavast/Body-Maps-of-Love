%% Normality assumptions
% for each pixel
addpath('external/bodyspm/')

cfg=[];
cfg.outdata = './output/body/';
cfg.datapath = './output/body/';

cfg.Nstimuli = 27;
cfg.Nempty = 1;
cfg.list = [cfg.outdata 'whitelist.txt'];
cfg.doaverage = 0;
cfg.Nbatches=5;
cfg.posneg=0;
cfg.overwrite=1;
cfg.niter=0; % do not run cluster correction


base=uint8(imread('bodySPM_base2.png'));
mask=uint8(imread('bodySPM_base3.png'));
mask=mask*.85;
base2=base(10:531,33:203,:);
in_mask=find(mask>128); % palauttaa maskista niiden indeksit jotka yli 128

subjects=textread(cfg.list,'%s');
Nsubj=length(subjects);
tocheck=zeros(Nsubj,1);

for ns=1:Nsubj;
    disp(['Processing subject ' subjects{ns} ' which is number ' num2str(ns) ' out of ' num2str(Nsubj)]);
    matname=[cfg.datapath '/' subjects{ns} '.mat'];
    load (matname) % 'resmat','resmat2','times'
    if(cfg.doaverage==0)
        tempdata=reshape(resmat,[],size(resmat,3)); % kolumnikohtaisesti uudelleenjärjestää
    else 
        tempdata=reshape(resmat2,[],size(resmat2,3));
    end
	if(ns==1)
        [ length(in_mask) Nsubj size(tempdata,2) ]
		alldata=zeros(length(in_mask),Nsubj,size(tempdata,2));
	end
	alldata(:,ns,:)=tempdata(in_mask,:); % tallentaa koehenkilön maskin sisällä olevan datan 
end

%% Are the pixel distributions normal?
reject = 0;
dontreject = 0;
count = 0;
alldata(isnan(alldata)) = 0;     % NaN's are interpreted as 0 (no report of sensation)
for i=1:size(alldata,1)
    for j=1:size(alldata,3)
        temp = alldata(i,:,j);
        temp = temp(~isnan(temp));
        if kstest(temp) == 1
            count = count + 1;
        else
            dontreject = dontreject + 1;
        end
    end
end