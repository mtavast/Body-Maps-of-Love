%% plot split group consistency (for pixels) results 
clc;clear all; close all
load output/body/split_sample.mat
FS = 22
Ticksize = 16
subplot(1,2,1)
body_consistency=zeros(size(mask));   
body_consistency(find(mask>128))=pvm; % pvm = Spearman correlation values
imagesc(mask)
hold on
imagesc((body_consistency),([.2 1])) 
hotmap=hot(64);
hotmap=[[1 1 1]/2;hotmap]
colormap(hotmap)
h=colorbar
set(h, 'FontSize', Ticksize)
ylabel(h,['Median split group consistency' newline '[Spearmann corr.]'], 'FontSize', FS)
axis equal
axis off

subplot(1,2,2)
body_consistency_median_p=zeros(size(mask));
body_consistency_median_p(find(mask>128))=-log10(pvmp+eps*eps); 

pvmQ=mafdr(pvmp,'BHFDR','true'); % BHFDR correction
imagesc((body_consistency_median_p))
colormap(hotmap)
h=colorbar

pticks=get(h,'Ticks');
set(h,'TickLabels',num2str(10.^-pticks'), 'FontSize', Ticksize)
axis equal
axis off

ylabel(h,['Median split group consistency p-values'], 'FontSize', FS)
    set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color',[1 1 1]);
print(gcf, 'figs/supplement/body_consistency_of_tool.png', '-dpng', '-r300')

%% Compute split-half reliability for order of bodies based on the strenght of embodiment
% This is code for computing Spearman correlations across 5000 split-halfs.
% The correlation is for the order of the strenght of the painting, that is, the overall amount of paint in each stimuli.
clc; close all; clear all;
addpath('./external');
base=uint8(imread('baseindividual.png'));
base2=base(10:531,33:203,:); % single image base
mask=imread('maskindividual.png'); 
in_mask=find(mask>128);

cfg.list = ['./output/body/whitelist.txt'];
subjects=textread(cfg.list,'%s');
Nsubj = length(subjects);
for ns=1:Nsubj;
    disp(['Processing subject ' subjects{ns} ' which is number ' num2str(ns) ' out of ' num2str(Nsubj)]);
    matname=['./output/body/' subjects{ns} '.mat'];                 
    load (matname) % 'resmat','times'                  % loads the data from output  
    tempdata=reshape(resmat,[],size(resmat,3));        
	alldata(:,ns,:)=tempdata(in_mask,:);               % only the data in mask
end

sums=zeros(27,Nsubj);
for s=1:Nsubj
    temp = sum(alldata(:,s,:), 'omitnan');
    sums(:,s) = reshape(temp,27,1);
end
% In the sums variable, the rows are the stimuli and the columns the
% participant. Each cell is the overall amount of paint (inside the mask) 
% for the given stimuli for the given participant. 

cfg.niter = 5000;
cfg.data = sums;
[icc2 iccpermu2 iccpvals2] = bramila_splitsample(cfg)
% The splitsample calculates rank correlation between average values of
% half of the data vs. average values of other half of the data.

%% Plot split-half group for amount of paint

hp=histogram(iccpvals2)
pbins=get(hp,'BinEdges');
subplot(2,1,1)
disp(['Min correlation ' num2str(min(icc2))])
disp(['Max correlation ' num2str(max(icc2))])
disp(['Min p-value ' num2str(min(iccpvals2))])
disp(['Max p-value ' num2str(max(iccpvals2))])
h=histogram(icc2,'BinLimits',[0.85,0.99]);
axis([0.85,0.99 0 500]);
xlabel('Split group consistency [Spearman correlation]');
ylabel('Number of iterations');
title('Consistency of the body painting 5000 split sample correlations of the total sum of paint');

temptick=0.85:0.01:.99; 
set(gca,'XTick',temptick);
set(gca,'Xticklabel',num2str(temptick'))

temptick=0:200:500;
set(gca,'YTick',temptick);
set(gca,'Yticklabel',num2str(temptick'))

set(gca,'FontSize',20)

subplot(2,1,2)
[hc xe ye]=histcounts2(icc2,iccpvals2,get(h,'BinEdges'),pbins);
imagesc(xe,ye,hc');

title('Heatmap of the p-values for each of the 5000 iterations')
xlabel('Split group consistency [Spearman correlation]')
ylabel('P-value');
axis xy
hold on

temptick=0.85:0.01:.99; 
set(gca,'XTick',temptick);
set(gca,'Xticklabel',num2str(temptick'));

temptick=get(gca,'YTick');
temptick= min(iccpvals2):0.0000005:max(iccpvals2); % 1:3
set(gca,'YTick',temptick);
set(gca,'Yticklabel',num2str(temptick'));
axis([min(icc2),max(icc2) [min(iccpvals2) max(iccpvals2)]]); 
set(gca,'Ylim',[min(iccpvals2) max(iccpvals2)]); 
set(gca,'FontSize',20);
colormap(hot)
hc=colorbar('west');
ylabel(hc,'Number of iterations')
set(hc,'Color',[1 1 1])
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 .75 1]);
set(gcf,'color',[1 1 1]);
print(gcf, 'figs/supplement/split_half_amountofpaint.png', '-dpng', '-r500')
disp(['Median Spearman correlation ', num2str(median(icc2))])
disp(['Median Spearman correlation p-values', num2str(median(iccpvals2))])