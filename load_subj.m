function data=load_subj(folder,option)

% load_subj(folder,option)
%   Slightly modified from:
%   Function from original embody study: https://version.aalto.fi/gitlab/eglerean/embody/-/blob/master/matlab/load_subj.m
%   
%   If you use it, please cite: Nummenmaa L., Glerean E., Hari R., Hietanen, J.K. (2014) Bodily maps of emotions, Proceedings of the National Academy of Sciences of United States of America doi:10.1073/pnas.1321664111
%   http://www.pnas.org/content/111/2/646.abstract


    list=csvread([folder '/' 'presentation_1.csv']); % we have only one batch here, so presentation_1 file is sufficient
    N=length(list);
    if(option==0)
        disp('here')
        for n=0:N-1;
            file=[folder '/' num2str(n) '.csv'];
            fid = fopen(file);
            line=textscan(fid,'%s','CollectOutput',1,'Delimiter',';');
            data(:,n+1)=line{1};
        end
        data=data';
    end
    if(option==1)
        disp('option 1 not implemented')
    end
    if(option>=2)
        for n=1:N;
            file=[folder '/' num2str(n) '_1.csv'];
            
            line=dlmread(file,',');
            delim=find(-1==line(:,1));
            data(n+1).mouse=line(1:delim(1)-1,:);
            data(n+1).paint=line((delim(1)+1):(delim(2)-1),:);
            data(n+1).mousedown=line((delim(2)+1):(delim(3)-1),:);
            data(n+1).mouseup=line((delim(3)+1):end,:);
        end
    end
end
