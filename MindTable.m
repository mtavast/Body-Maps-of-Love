clear all
close all
addpath('./external/table2latex/')
load('output/mind/M10_data')
load_labels;
sensations_classes;
sensations=labels_en;
Ndim=6;
dim_labels={ %'Bodily sensation strength','Mind sensation strength','Emotion intensity','Agency','Last time'};
    'Bodily saliency' % How much feels in body question
    'Mental saliency' % How much involves the mind
    'Valence' % emotional valence question
    'Controllability' % how much you can control question (used to be agency)
    'Last experienced' % last time question
    'Touch' % how strongly associate with touch
    };

mean_data=zeros(27,6);
sem_data=zeros(27,6);
mean_dataZ=zeros(27,6);
sem_dataZ=zeros(27,6);

% compute means, medians and distributions
for s=1:length(sensations)
    thisdata=squeeze(data(s,:,:)); 
    ids=find(~isnan(max(thisdata)));
    thisdata=thisdata(:,ids); 
    mean_data(s,:)=mean(thisdata,2); 
    median_data(s,:)=median(thisdata,2); 
    Nitems(s,1)=size(thisdata,2);
    std_data(s,:)=std(thisdata,0,2); 
    sem_data(s,:)=std(thisdata,0,2)./sqrt(size(thisdata,2));
    sens_count(s,1)=size(thisdata,2);
end

mean_data = round(mean_data);
median_data = round(median_data);
sem_data = round(sem_data);
std_data = round(std_data);

datatable = table();
datatable(:,1) = labels_en;
for i=1:6
    datatable(:,end+1) = array2table(mean_data(:,i));
    datatable(:,end+1) = array2table(std_data(:,i));
    datatable(:,end+1) = array2table(sem_data(:,i));
end

datatable.Properties.VariableNames = {'Labels',
    'M Bodily saliency', 
    'SD Bodily saliency', 
    'SE Bodily saliency', 
    'M Mental saliency', 
    'SD Mental saliency',
    'SE Mental saliency',
    'M Valence', 
    'SD Valence', 
    'SE Valence', 
    'M Controllability',
    'SD Controllability', 
    'SE Controllability', 
    'M Last experienced',
    'SD Last experienced', 
    'SE Last experienced', 
    'M Touch',
    'SD Touch',
    'SE Touch'};

table2latex(datatable)
%% Mean table

mean_data_with_labels = array2table(mean_data);
mean_data_with_labels(:,7) = labels_en;

mean_data_with_labels.Properties.VariableNames = { %'Bodily sensation strength','Mind sensation strength','Emotion intensity','Agency','Last time'};
    'Bodily saliency' % How much feels in body question
    'Mental saliency' % How much involves the mind
    'Valence' % emotional valence question
    'Controllability' % how much you can control question (used to be agency)
    'Last experienced' % last time question
    'Touch' % how strongly associate with touch
    'Labels'};

% median table

median_data_with_labels = array2table(median_data);
median_data_with_labels(:,7) = labels_en;

median_data_with_labels.Properties.VariableNames = { %'Bodily sensation strength','Mind sensation strength','Emotion intensity','Agency','Last time'};
    'Bodily saliency' % How much feels in body question
    'Mental saliency' % How much involves the mind
    'Valence' % emotional valence question
    'Controllability' % how much you can control question (used to be agency)
    'Last experienced' % last time question
    'Touch' % how strongly associate with touch
    'Labels'};

