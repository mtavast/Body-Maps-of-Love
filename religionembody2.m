clear all
close all

% choose the group
group = 'religious'
%group = 'atheists'

%%
cfg.outdata = './output/religion/';
cfg.datapath = '/m/nbe/scratch/love/Love final data/Rakkaus_keho/subjects';

cfg.Nstimuli = 27;
cfg.Nempty = 1;
cfg.list = [cfg.outdata 'body_list_' group '.txt']; 

cfg.doaverage = 0;
cfg.Nbatches=1
cfg.posneg=1;
cfg.overwrite=1;

bspm = bodySPM_preprocess_batch(cfg); 

save([cfg.outdata '/bspm_ ' group '.mat']) 

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
npainted=sum(~isnan(squeeze(bspm.allTimes(:,2,:))));
figure(1)
subplot(2,2,1)
stem(npainted);
subplot(2,2,2)
stem(sort(npainted));


ids=find(bspm.tocheck<=14);
ids=find(npainted>=15);
subjects=textread(cfg.list,'%s');
out=[];
for i = 1:length(ids)
    out{i,1}=subjects{ids(i)};
end

fileID=fopen([cfg.outdata '/whitelist_' group '.txt'],'w'); 
for i=1:length(out);
    if(~strcmp(out{i}(1),'F'))
        disp([out{i} ' has invalid ID for this study']);
        continue;
    end
    fprintf(fileID,'%s\n',[out{i}]);
end
fclose(fileID);

