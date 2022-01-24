%% Set up for calculating the sum of paint
% Loads each subjects data directly
clc
clear all
close all

%group = 'religious'
group = 'atheists'
%%

cfg.datapath = '/m/nbe/scratch/love/Love final data/Rakkaus_keho/subjects';
lista = ['./output/religion/whitelist_' group '.txt'];

whitelist = textread(lista,'%s');
NSu = length(whitelist);
f = {'name'};
whitelist = cell2struct(whitelist,f,NSu);

% only final subjects
subjects = whitelist;

% the base image used for painting 
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
        h=fspecial('gaussian',[15 15],5);
        over=imfilter(over,h);
        % no substraction, as there is only one body
        over2=over(10:531,33:203,:);
        resmat(:,:,n)=over2;
    end
  
resmats.data{s} = resmat;
resmats.subject(s)= s;

end

save(['./output/religion/resmats_' group '.mat'], 'resmats')

%% Loop to calculate means

NSti = 27
sumvalues = zeros(NSti, NSu);
for n=1:NSti;                           % every stimulus
    for s=1:NSu;                        % every subject
    temp = resmats.data(s);             % paint data for each subject and stimulus
    temp = temp{1,1};                   % same as above
    temp_stimulus = temp(:,:,n);
    sumvalues(n,s) = sum(temp_stimulus, 'All');    % add all values
    end
end

means = mean(sumvalues, 2) % rowmeans

save(['./output/religion/means_' group '.mat'], 'means')

