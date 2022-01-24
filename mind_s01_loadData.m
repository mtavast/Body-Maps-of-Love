%% load labels and categories, screening
clear all
close all

datapath='/m/nbe/scratch/love/Love final data/Rakkaus_arvio/subjects';
outdatapath='./data/mind/M10/';
prefix='FI';
    
sensations=textread('sensations.txt','%s','delimiter','\n');

    list=dir([datapath '/' prefix '*']);
    
    % blacklist subjects who did not answer to all questions, with age == 0, and non-native Finnish
    % speakers
    
    blacklist=[];
    
    Nvariables=8;
%   array contains 8 variables
%   sex : 0 (male) 1 (female) 2 (other or do not want to define)
%   age : 0 (no answer), 1 (18-25), ... 6 (66 or more)
%   mothertongue : 0 (Finnish), 1 (other)
%   education: 0 (perusaste), 1 (keskiaste), 2 (korkea-aste)
%   psych: 0 (no), 1 (yes)
%   children: 0 (no), 1 (yes)
%   religion: 0 (atheist), 1 (agnostic), 2 (believer)
%   orientation: 0 (heterosexual), 1 (homosexual), 2 (bisexual), 3 (other or do not want to
%   define)
    
    behav_data=zeros(length(list),Nvariables+1);
    for id=1:length(list)
        thisdata=zeros(1,Nvariables+1);
        
        temp=textread([datapath '/' list(id).name '/data.txt' ],'%s','delimiter','\n'); % behavioral data
        resplist=dir([datapath '/' list(id).name '/*.csv']);   % all csv files in the subjects folder
        resp={};
        for rl=1:length(resplist)
            if(regexp(resplist(rl).name,'^[0-9]')) % if the file name starts with number this is true
                resp{end+1}=resplist(rl).name;     % then add to the list
            end
        end
        % number of responses (need to have answered to all)
        Nresp=length(resp);
        thisdata(end)=Nresp;
        if(Nresp < 27)
            disp(['Skipping ' list(id).name ' (not enough responses)']);
            blacklist=[blacklist;id];  
            
        end
        
        % go through each line
        
        for t=1:Nvariables % 8 variables 
                thisdata(1,t)=str2num(temp{t});
        end
        
        if(thisdata(2)<1) % filter out the ones that did not respond to age
            disp(['Skipping ' list(id).name ' (age is zero)']);
            blacklist=[blacklist; id];
            
        end
        
        % to remove non-native Finnish speakers 
        if(thisdata(3)==1) 
            disp(['Skipping ' list(id).name ' (not native Finnish speaker)']);
            blacklist=[blacklist; id];
            
        end    
        behav_data(id,:)=thisdata;
    end
    
    blacklist=unique(blacklist);
    whitelist=ones(length(list),1);
    whitelist(blacklist)=0;
    
    good_subj_list=list(find(whitelist));
    Nsubj = length(good_subj_list);

count = 0
for i=1:length(behav_data)
    if behav_data(i,9) == 27
        count = count +1;
    end
end
disp(['Number of participants who completed the dimension experiment ', num2str(count)])

    %% load the actual data
    % data = words x dimensions x good subjects
    % 6 dimensions listed here:
    Ndim=6;
    dim_labels={ %'Bodily sensation strength','Mind sensation strength','Emotion intensity','Agency','Last time', 'Touch'};
		'Bodily saliency' % How much feels in body question
		'Mental saliency' % How much feels in the mind
		'Valence' % emotional valence question
		'Controllability' % how much you can control question 
		'Last experienced' % last time question
        'Touch' % how strongly associate with touch
	};	
    data=NaN(length(sensations),Ndim,Nsubj); 
    for s=1:Nsubj
        resplist=dir([datapath '/' good_subj_list(s).name '/*.csv']); % subject's data
        resp={};
        for rl=1:length(resplist)
            if(regexp(resplist(rl).name,'^[0-9]'))                    % if starts with number
                resp{end+1}=resplist(rl).name;                        % what stimulus
                respIDtemp=strsplit(resplist(rl).name,'_');           % 
                respID=str2num(respIDtemp{1});                        %
                thisdata=load([datapath '/' good_subj_list(s).name '/' resplist(rl).name]);
                data(respID,:,s)=thisdata(end,:); % it needs the 'end' because it's the last answer recorded for when subj hit back button
            end
        end
    end
    
    save(['output/mind/M10_data'], 'data', 'sensations', 'dim_labels','behav_data','good_subj_list','blacklist','whitelist');


%% some subject stats
clear all
close all
load output/mind/M10_data.mat
disp(['Initial number of participants: ' num2str(size(behav_data,1))])

behav=behav_data(find(whitelist),:);
disp(['Number of participants after screening for age, native language and number of responses: ' num2str(size(behav,1))])
disp(['Males/females/other or dont want to define: ' num2str(length(find(behav(:,1)==0))) '/' num2str(length(find(behav(:,1)==1))) '/' num2str(length(find(behav(:,1)==2)))])
disp(['18-25: ' num2str(length(find(behav(:,2)==1)))])
disp(['26-35: ' num2str(length(find(behav(:,2)==2)))])
disp(['36-45: ' num2str(length(find(behav(:,2)==3)))])
disp(['46-55: ' num2str(length(find(behav(:,2)==4)))])
disp(['56-65: ' num2str(length(find(behav(:,2)==5)))])
disp(['66 or more: ' num2str(length(find(behav(:,2)==6)))])
disp(['perus/keski/korkea: ' num2str(length(find(behav(:,4)==0))) '/' num2str(length(find(behav(:,4)==1))) '/' num2str(length(find(behav(:,4)==2)))])
disp(['# of parents: ' num2str(length(find(behav(:,6)==1)))])
disp(['Atheist/agnostic/religious: ' num2str(length(find(behav(:,7)==0))) '/' num2str(length(find(behav(:,7)==1))) '/' num2str(length(find(behav(:,7)==2)))])
age_Q=prctile(behav(:,2),[0 25 50 75 100]);
disp(['Age (min, q1, median, q2, max): ' num2str(age_Q)])
