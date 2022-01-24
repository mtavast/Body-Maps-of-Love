clear all
close all

addpath(genpath('/m/nbe/scratch/braindata/eglerean/code/bodyspm/'));

cfg.outdata = './output/religion';
cfg.datapath = '/m/nbe/scratch/love/Love final data/Rakkaus_keho/subjects';
cfg.Nstimuli = 27;
cfg.Nempty = 1;
cfg.phenodata = 0;

bspm=bodySPM_parseSubjects_batch(cfg);

groups = {'religious', 'atheists'}

%% filter subjects and output report
for member = groups
    member
    disp( '---------------------------------')
    disp( '>>> Summary statistics report <<<')
    disp( '---------------------------------')
    disp(['    Total number of subjects who registered for the experiment: ' num2str(size(bspm.data_filter,1))])
    ids=find(bspm.data_filter(:,1)>=33); % Counts all the files in the folder. Needs to be 33 if the subject has completed all stimuli.
    disp(['    Total number of subjects who completed the whole task: ' num2str(size(ids,1))])
    Nvariables=8;
    % behav variables

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
            %disp(['Skipping ' list(id).name ' (no age registered)']);
            blacklist=[blacklist; i];

        end

        if(thisdata(3)==1) % filter not native Finnish speakers
             %disp(['Skipping ' list(id).name ' (not native Finnish speaker)']);
             blacklist=[blacklist; i]; % lisätty i

        end
        
        % One subject was removed based on data quality, subject {'FI990065'}
        %  The subject in question drew question marks on couple of the bodies.
        % Individual maps can be recreated with the code at the end of body_s02_preprocessing script.
        if strcmp(out{i, 1}, 'FI990065') == 1 
             disp(['Skipping due to paint data quality']);
             blacklist=[blacklist; i]; % lisätty i
        end 

        
        if isequal(member{1}, 'religious')
            if(thisdata(7)<2) % religious are 2, others are blacklistet
                 %disp(['Skipping ' list(id).name ' (not native Finnish speaker)']);
            blacklist=[blacklist; i]; % lisätty i
            end
        end
        
        if isequal(member{1}, 'atheists')
            if(thisdata(7)>0) % atheists are 0, others are blacklisted
                %disp(['Skipping ' list(id).name ' (not native Finnish speaker)']);
            blacklist=[blacklist; i]; % uncomment for atheists
            end
        end 
    
        behav_data(i,:)=thisdata;
    end

    if(length(blacklist)==0) 
        disp(['    All subjects registered their age (minimum age = ' num2str(min(behav_data(:,2))) ')'])
    else
       out(blacklist)=[];
       behav_data(blacklist,:)=[];
    end

    save([cfg.outdata '/religion_data_' member{1} '.mat'],'behav_data','behav_labels','blacklist','out')
    disp(length(out));
    fileID=fopen([cfg.outdata '/body_list_' member{1} '.txt'],'w');
    for i=1:length(out); 
        if(~strcmp(out{i}(1),'F')) 
            disp([out{i} ' has invalid ID for this study']); 
            continue;
        end
        fprintf(fileID,'%s\n',[out{i}]);
    end
    fclose(fileID);
end