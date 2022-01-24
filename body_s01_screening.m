clear all
close all

addpath('./external/bodyspm')

cfg.outdata = './output/body';
cfg.datapath = '/m/nbe/scratch/love/Love final data/Rakkaus_keho/subjects';
cfg.Nstimuli = 27; % In the bodySPM_parseSubjects function, it would flag participants with lower than NStiumuli number files in directory. Full subjects here need to have 33 files, but we don't have a single subject with filenumers 27-32.
cfg.Nempty = 1; 
cfg.phenodata = 0;

bspm=bodySPM_parseSubjects_batch(cfg);


%% filter subjects and output report

disp( '---------------------------------')
disp( '>>> Summary statistics report <<<')
disp( '---------------------------------')
disp(['    Total number of subjects who registered for the experiment: ' num2str(size(bspm.data_filter,1))])
ids=find(bspm.data_filter(:,1)>=33); % Counts the number of csv files in the folder. Needs to be 33 if the subject has completed the task.
disp(['    Total number of subjects who completed the whole task: ' num2str(size(ids,1))])
Nvariables=8;
% behav variables
%   1: sex : 0 (male) 1 (female) 2 (other or do not want to answer)
%   2: age : 0 (no answer), 1 (18-25), ... 6 (66 or more)
%   3: mothertongue : 0 (Finnish), 1 (other)
%   4: education: 0 (perusaste), 1 (keskiaste), 2 (korkea-aste)
%   5: psychologist: 0 (no), 1 (yes)
%   6: children: 0 (no), 1 (yes)
%   7: religion: 0 (atheist), 1 (agnostic), 2 (religious)
%   8: orientation: 0 (hetero), 1 (homo), 2 (bi), 3 (other or do not want to answer)
behav_labels={
    'gender' %(F=1)
    'age'
    'mothertongue'
    'education'
    'psychologist'
    'children'
    'religion' 
    'orientation'
};

behav_data=zeros(length(ids),Nvariables); % luo tyhjän listan
blacklist=[];
for i = 1:length(ids) % length(ids) lopullisten koehenkilöiden määrä
        out{i,1}=bspm.subjects(ids(i)).name;
    % load phenodata
   
    thisdata=zeros(1,Nvariables);
    
    temp=textread([cfg.datapath '/' out{i,1} '/data.txt' ],'%s','delimiter','\n');   

    % go through each line
    
    for t=1:Nvariables % 8 variables
            thisdata(1,t)=str2num(temp{t});
    end
    
    if(thisdata(2)<1) % changed to 1 
        disp(['Skipping (no age registered)']);
        blacklist=[blacklist; i];
        
    end
    
    if(thisdata(3)==1) % filter not native Finnish speakers
         disp(['Skippining (not native Finnish speaker)']);
         blacklist=[blacklist; i]; % lisätty i
            
    end 
    
% Only one subject was removed based on data quality, subject {'FI990065'}
% The subject in question drew question marks on couple of the bodies.
% Individual maps can be recreated with the code at the end of body_s02_preprocessing script.
   
    if strcmp(out{i, 1}, 'FI990065') == 1 
         disp(['Skipping due to paint data quality']);
         blacklist=[blacklist; i]; % lisätty i
    end 

    behav_data(i,:)=thisdata;
end


if(length(blacklist)==0) 
    disp(['All subjects registered their age'])
else
   out(blacklist)=[];
   behav_data(blacklist,:)=[];
end

save([cfg.outdata '/behav_data.mat'],'behav_data','behav_labels','blacklist','out')

fileID=fopen([cfg.outdata '/body_list.txt'],'w'); % writes the final subject IDs
for i=1:length(out); 
	if(~strcmp(out{i}(1),'F')) 
		disp([out{i} ' has invalid ID for this study']); 
		continue;
	end
	fprintf(fileID,'%s\n',[out{i}]);
end
fclose(fileID);
