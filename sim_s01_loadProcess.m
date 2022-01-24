clear all
close all

addpath('/m/nbe/scratch/braindata/shared/toolboxes/cbrewer/');

datapath='/m/nbe/scratch/love/Love final data/Rakkaus_sanat/subjects';
countryfile='/m/nbe/scratch/braindata/tavastm1/Backups/LoveDatabackups/countries.txt';
prefix='FI';
load_labels
sensations=labels;

list=dir([datapath '/' prefix '*']);

% IMPORTANT NOTE: there is a subject that causes error to preprocessing below. 
% The subject has no data at all, i.e. no outdata.csv or presentation file. 
% Don't know the reason for this, maybe the subject close the browser before the first stimuli?
% I interpret that the subject did not even start the task.
% Here, the participant is removed manually, so that the preprocessing script does not cause errors. 
list(70) = [];   % the participant with no data is #70 in our list.

% blacklist subjects with age = 0, non finnish speakers or with no responses
blacklist=[];

Nvariables=8;
% array contains 8 variables
%   sex : 0 (male) 1 (female) 2 (other or do not want to answer)
%   age : 0 (no answer), 1 (18-25), ... 6 (66 or more)
%   mothertongue : 0 (Finnish), 1 (other)
%   education: 0 (perusaste), 1 (keskiaste), 2 (korkea-aste)
%   psychologist: 0 (no), 1 (yes)
%   children: 0 (no), 1 (yes)
%   religion: 0 (atheist), 1 (agnostic), 2 (religious)
%   orientation: 0 (heterosexual), 1 (homosexual), 2 (bisexual), 3 (other or do not want to
%   answer)

behav_data=zeros(length(list),Nvariables+2);

for id=1:length(list) % jokaiselle koehenkilolle
    thisdata=zeros(1,Nvariables+2); %luo tyhjan listan
    
    temp=textread([datapath '/' list(id).name '/data.txt' ],'%s','delimiter','\n'); % hakee datan, tsekattu, etta oikein yhella khlla
    resplist=dir([datapath '/' list(id).name '/*.csv']); % luo struct, joka vie itse dataan
    resp={};
    for rl=1:length(resplist)
        if(regexp(resplist(rl).name,'^[0-9]'))
            resp{end+1}=resplist(rl).name;
        end
    end
    
    
    % go through each line
    
    for t=1:Nvariables % 11 variables
             thisdata(1,t)=str2num(temp{t}); % tuo tempin thisdataan       
    end
    
    if(thisdata(2)<1) %  filter out the ones that did not respond to age
        disp(['Skipping ' list(id).name ' (age missing)']);
        blacklist=[blacklist; id]; % lisaa idn uudelle riville
        
    end
            % to remove non-native Finnish speakers 
    if(thisdata(3)==1) % filter out the ones that responded 1 (other) to native languange
         disp(['Skipping ' list(id).name ' (not native Finnish speaker)']);
         blacklist=[blacklist; id];
            
    end
     
    behav_data(id,:)=thisdata; 
end

blacklist=unique(blacklist);
whitelist=ones(length(list),1);
whitelist(blacklist)=0; 
behav_data_unfiltered=behav_data;
behav_data(blacklist,:)=[]; % remove thos in blacklist
good_subj_list=list(find(whitelist));
Nsubj = length(good_subj_list);

%% load the actual data 
% 
DOPLOT=0;
moved=zeros(27,Nsubj);
for s=1:Nsubj
    fid = fopen([datapath '/' good_subj_list(s).name '/outdata.csv'],'rt');
    indata = textscan(fid, '%s', 'HeaderLines',1);
    fclose(fid);
    yourdata = indata{1};
    if(length(yourdata)==0) 
        disp('found a bad subject')
        
        continue; 
    end
    
    datastart=str2num(yourdata{1}); % used only for the timestamp
    
    data=str2num(yourdata{end}); % ottaa viimeisen rivin, joka sisältää lopullisen datan
    x=data(1:2:54); % changed from 200
    y=data(2:2:54); % changed from 200
    
    
    % QC percentage of X that are under pixel coordinate 350
    behav_data(s,9)=sum(x>350); % how many items were moved inside the square
    behav_data(s,10)=data(end)-datastart(end); % seconds from beginning to end
    behav_data(s,11)=15*length(yourdata); % time spent on task in seconds. This is an approximation with an error of 15 seconds
    
    
    if(DOPLOT)
        for t=1:length(yourdata)
            data=str2num(yourdata{t});
            x=data(1:2:54);
            y=data(2:2:54);
            figure(s);
            plot(x,y,'.')
            axis([0 1000 0 1000])
            pause(1/25)

        end
        
          
    end
    moved(find(x>=350),s)=1;        % jos ruudussa = 1
    y(find(x<350))=NaN;              %IGNORATAAN JOS EI RUUDUSSA, NaN Ignorataan myöhemmin laskuissa
    x(find(x<350))=NaN;
    temp=squareform(pdist([x' y'])); % Euclidinen etäisyys
    temp=temp/max(abs(temp(:)));     % Skaalaataan välille 0-1, niin että suurin etäisyys =1
    simmat(:,:,s)=temp;
    
end

% more subject cleaning
new_blacklist=find(behav_data(:,9)<27); % Pitää olla liikuttanut kaikkia sanoja

good_subj_list(new_blacklist)=[];
behav_data(new_blacklist,:)=[];
simmat(:,:,new_blacklist)=[]; % Poistaa ne, jotka eivät olleet

%% some subject stats

behav=behav_data;
disp(['Number of participants after screening for age and number of words moved: ' num2str(size(behav,1))])
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
disp(['Hetero/gay/bi/other or dont want to define: ' num2str(length(find(behav(:,8)==0))) '/' num2str(length(find(behav(:,8)==1))) '/' num2str(length(find(behav(:,8)==2))) '/' num2str(length(find(behav(:,8)==3)))])
age_Q=prctile(behav(:,2),[0 25 50 75 100]);
disp(['Age (min, q1, median, q2, max): ' num2str(age_Q)])


        
%% Behavioral correlations

% check if behaviour explains any distance between items
R=ones(27,27,11); 
P=ones(27,27,11);
for s1=1:27;
    for s2=(s1+1):27;
        temp=squeeze(simmat(s1,s2,:)); % Kaikkien koehenkilöiden ärsyke s ja äryke s+1 välinen etäisyys
        [rr pp]=corr(temp(~isnan(temp)),behav_data(~isnan(temp),:),'type','spearman'); % järjestyskorrelaatio etäisyyden ja beh. välillä
        R(s1,s2,:)=rr; % korrelaatiokertoimet
        R(s2,s1,:)=rr; 
        P(s1,s2,:)=pp; % p-arvot 
        P(s2,s1,:)=pp; 
	end
end

map=cbrewer('div','RdBu',9)
for bID=1:11 % jokaisesta behavioraalisesta muuttujasta
    figure(100)
    subplot(2,7,bID)
    rmat=R(:,:,bID); % korrelaatiokerroinmatriisi ko. behavioraalisen muuttujan ja kaikkien ärsykkeiden välillä
    pmat=P(:,:,bID); % sama kuin yllä, mutta p-arvoille
    q=mafdr(pmat(find(triu(ones(27),1))),'BHFDR','true'); % FDR-korjaus
    mask=zeros(27);
    mask(find(triu(ones(27),1)))=double(q<=0.05); % näyttää alle 0.05
    Nposvals=sum(mask(:));
    
    mask=mask+mask';
    
    imagesc(rmat.*mask,[-.6 .6]) 
    colormap(map)
    axis square
    title(num2str(Nposvals)) % otsikoksi tulee positiivisten määrä FDR-korjauksen jälkeen
    if(Nposvals>0 )          % tehdään erikseen kuvat niistä, joissa positiivisia korrelaatioita
        figure(bID)
        imagesc(rmat.*mask,[-.6 .6])
        colormap(map)
        axis square
        title(num2str(Nposvals))
        tempIDs=find(sum(mask)>0);
        set(gca,'XTick',tempIDs)
        set(gca,'XTickLabel',sensations(  tempIDs))
        set(gca,'XTickLabelRotation',90);
        set(gca,'YTick',tempIDs)
        set(gca,'YTickLabel',sensations(  tempIDs))
        colorbar
    end
    
end





msim=nanmean(simmat,3);
medsim=nanmedian(simmat,3);
vsim=nanvar(simmat,[],3);


p=squareform(msim);



save output/sim/sim_data.mat simmat good_subj_list behav_data behav_data_unfiltered

error('stop')  

%% Hierarcical clustering and multidimensional scaling

z=linkage(p,'complete');
dendrogram(z,0,'orientation','left','labels',sensations)
sensations = labels_en
figure(123)
[Y,e] = cmdscale(msim);
for n=1:length(sensations)
    
    plot(Y(n,1), Y(n,2),'o','Color','Red','MarkerSize',10, 'MarkerFaceColor', 'Red');
    hold on
    text(Y(n,1), Y(n,2)+sign(randn)*0.005,sensations{n},'Color','Red');
    
end
grid on

figure(1234)
for n=1:length(sensations)
    
    plot3(Y(n,1), Y(n,2),Y(n,3),'o','Color','Red','MarkerSize',10,'MarkerFaceColor','Red');
    hold on
    text(Y(n,1), Y(n,2),Y(n,3)+sign(randn)*0.005,sensations{n},'Color','Red');
    
    
end
grid on

%% The script above removes the participants based on demographic variables before based on if they have completed the task
% This can be used to check how many participants would
% be removed based ONLY on demographic criteria (4 participants).

% Number of words moved by participants who were excluded based on languange/age:
%Skipping FI (not native Finnish speaker): 27 
%Skipping FI(age missing): 6 
%Skipping FI (not native Finnish speaker): 27
%Skipping FI (not native Finnish speaker): 27
%Skipping FI (not native Finnish speaker): 10
%Skipping FI (not native Finnish speaker): 27
%Skipping FI (not native Finnish speaker): 1 
%Skipping FI (not native Finnish speaker): 0 
%Skipping FI (not native Finnish speaker): 0
%Skipping FI (not native Finnish speaker): 10

    fid = fopen([datapath '/' 'FI' '/outdata.csv'],'rt'); % put the participant numbers here
    indata = textscan(fid, '%s', 'HeaderLines',1);
    fclose(fid);
    yourdata = indata{1};

    datastart=str2num(yourdata{1}); % used only for the timestamp
    
    data=str2num(yourdata{end}); 
    x=data(1:2:54); 
    y=data(2:2:54); 
    sum(x>350) % number of moved words

