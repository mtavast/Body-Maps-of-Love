%% plot split group consistency of tool
clear all
close all
load output/sim/split_sample.mat
disp(['Min correlation ' num2str(min(icc))])
disp(['Max correlation ' num2str(max(icc))])
disp(['Min p-value ' num2str(min(icc_pvals_mantel))])
disp(['Max p-value ' num2str(max(icc_pvals_mantel))])
hp=histogram(icc_pvals_mantel)
pbins=get(hp,'BinEdges');
subplot(2,1,1)
h=histogram(icc,'BinLimits',[0.87,0.97]) ;
axis([0.87,0.97 0 500])  ;
xlabel('Split group consistency [Spearman correlation]');
ylabel('Number of iterations');
title('Consistency of the similarity tool 5000 split sample correlations of the means')

temptick=0.87:0.01:.97; 
set(gca,'XTick',temptick);
set(gca,'Xticklabel',num2str(temptick'));

temptick=0:200:500;
set(gca,'YTick',temptick);
set(gca,'Yticklabel',num2str(temptick'));

set(gca,'FontSize',20);

subplot(2,1,2);
[hc xe ye]=histcounts2(icc,icc_pvals_mantel,get(h,'BinEdges'),pbins);
imagesc(xe,ye,hc');

title('Heatmap of the p-values for each of the 5000 iterations')
xlabel('Split group consistency [Spearman correlation]')
ylabel('P-value (mantel test)');
axis xy
hold on

temptick=0.87:0.01:.97; 
set(gca,'XTick',temptick);
set(gca,'Xticklabel',num2str(temptick'));

temptick=get(gca,'YTick');
temptick=(1:10).*10^-8; % 1:3
set(gca,'YTick',temptick);
set(gca,'Yticklabel',num2str(temptick'))
axis([0.87,0.97 [1e-08 9.2254e-08]]) 
set(gca,'Ylim',[1e-08 9.2254e-08])
set(gca,'FontSize',20)
colormap(hot)
hc=colorbar('west');
ylabel(hc,'Number of iterations')
set(hc,'Color',[1 1 1])
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 .75 1]);
set(gcf,'color',[1 1 1]);

disp(['Median Spearman correlation ', num2str(median(icc))])
disp(['Median Spearman correlation p-values Mantel test', num2str(median(icc_pvals_mantel))])

print(gcf, 'figs/supplement/sim_tool_consistency.png', '-dpng', '-r500')

%% Individual level correlations

clear all
close all

load output/sim/sim_data.mat
Nsubj=size(simmat,3);                   
subj_sim=eye(Nsubj);                    
idsTop=find(triu(ones(27),1));          
subj_items=zeros(Nsubj);                
subj_p=zeros(Nsubj);
for s1=1:Nsubj                         
    disp(num2str(s1/Nsubj))
    for s2=(s1+1):Nsubj                 
        temp1=simmat(:,:,s1);           
        temp2=simmat(:,:,s2);
        ids=find(triu(~isnan(temp1.*temp2),1));
        
        subj_items(s1,s2)=length(ids);
        [r p]=corr(temp1(ids),temp2(ids),'type','Spearman');   
        subj_sim(s1,s2)=r;                 % save the correlation coefficient
        subj_sim(s2,s1)=r;
        subj_p(s1,s2)=p;                   % save the pvalues
        subj_p(s2,s1)=p;
    end
end
subj_D=1-subj_sim;

z=linkage(squareform(subj_D),'complete');
[H,T,OUTPERM] = dendrogram(z,0); % H = 
figure(123)
imagesc(subj_sim(OUTPERM,OUTPERM),[0 0.7]);colorbar

q=mafdr(subj_p(find(triu(ones(Nsubj),1))),'BHFDR','true');    % FDR correctiom

qmat=zeros(Nsubj);
qmat(find(triu(ones(Nsubj),1)))=q;
qmat=qmat+qmat';
qmask=qmat<0.05;                                              

figure(4)
imagesc(qmask(OUTPERM,OUTPERM));colorbar                      

%% correlations

m  = find(triu(ones(Nsubj),1))
v  = subj_sim(m).';

histogram(v)
disp(['Median correlation ' num2str(median(v))])
%% p-values

m  = find(triu(ones(Nsubj),1))
v  = subj_p(m).';

histogram(v)
disp(['Mean p ' num2str(mean(v))])
disp(['Median p ' num2str(median(v))])
disp(['Percent of p values under 0.05 ' num2str(sum(v < 0.05) / length(v) * 100)])
%%
addpath('./external/')

clear all
close all
if(0)
    load output/sim/sim_data.mat
    Nsubj=size(simmat,3);                   
    subj_sim=eye(Nsubj);                    
    idsTop=find(triu(ones(27),1));          
    subj_items=zeros(Nsubj);                
    subj_p=zeros(Nsubj);
    correlations=zeros(1,1)
    pvalues=zeros(1,1)
    for s1=1:Nsubj                          
        disp(num2str(s1/Nsubj))
        for s2=(s1+1):Nsubj                 
            temp1=simmat(:,:,s1);           
            temp2=simmat(:,:,s2);
            [rtemp icc_pvals_mantel]=bramila_mantel(temp1,temp2,5000,'spearman'); 
            correlations = [correlations; rtemp];
            pvalues = [pvalues; icc_pvals_mantel];
        end
    end

    correlations(1) = []; % delete the first zero
    pvalues(1) = []; 

    disp(['Median Spearman correlation ', num2str(median(correlations))])
    disp(['Median Spearman correlation p-values Mantel test', num2str(median(pvalues))])

    save(['./output/sim/individual_correlations_mantel_test.mat'], 'correlations', 'pvalues')
end

load('./output/sim/individual_correlations_mantel_test.mat')
disp(['Median Spearman correlation ', num2str(median(correlations))])
disp(['Median Spearman correlation p-values Mantel test', num2str(median(pvalues))])
