%% Plot the strength of the embodiment
clc; clear all; close all;
% load the sums
religious_means = open('./output/religion/means_religious.mat');
atheists_means = open('./output/religion/means_atheists.mat');

means_data(:,1) = religious_means.means;
means_data(:,2) = atheists_means.means;

load_labels
labels = labels_en 
addpath('/m/nbe/scratch/braindata/shared/toolboxes/cbrewer/');

figure(2)
hold on

for s=1:27
    plot(means_data(s,1), means_data(s,2))
end

hold on

for i=1:27
    plot(means_data(i,1),means_data(i,2),'o', 'MarkerSize', 4, 'MarkerFaceColor', 'k','MarkerEdgeColor', 'k' )
    if i == 22 | i==9 | i==23 | i==11 | i==6 | i==27 | i==15
        text(means_data(i,1)-5,means_data(i,2), labels_en{i},'Rotation',0,'FontSize',12, 'HorizontalAlignment','right')
    elseif i == 1
        text(means_data(i,1)+5,means_data(i,2)-5, labels_en{i},'Rotation',0,'FontSize',12, 'HorizontalAlignment','right')
    elseif i==24 
        text(means_data(i,1)-3,means_data(i,2)-5, labels_en{i},'Rotation',0,'FontSize',12, 'HorizontalAlignment','right')
    elseif i == 4
        text(means_data(i,1)+7,means_data(i,2)+8, labels_en{i},'Rotation',0,'FontSize',12, 'HorizontalAlignment','right')
    elseif i == 26
        text(means_data(i,1)+5,means_data(i,2)-5, labels_en{i},'Rotation',0,'FontSize',12)
    elseif i == 13
        text(means_data(i,1)+5,means_data(i,2), labels_en{i},'Rotation',0,'FontSize',12)
    elseif i == 14
        text(means_data(i,1)-5,means_data(i,2), labels_en{i},'Rotation',0,'FontSize',12, 'HorizontalAlignment','right')
    elseif i == 18
        text(means_data(i,1)+5,means_data(i,2)-6, labels_en{i},'Rotation',0,'FontSize',12)
    elseif i == 19
        text(means_data(i,1)+5,means_data(i,2)-7, labels_en{i},'Rotation',0,'FontSize',12)
    elseif i == 5
        text(means_data(i,1)+5,means_data(i,2)-5, labels_en{i},'Rotation',0,'FontSize',12)
    elseif i == 25
        text(means_data(i,1)-5,means_data(i,2), labels_en{i},'Rotation',0,'FontSize',12, 'HorizontalAlignment','right')
    elseif i == 12
        text(means_data(i,1)-5,means_data(i,2)+4, labels_en{i},'Rotation',0,'FontSize',12, 'HorizontalAlignment','right')
    else
        text(means_data(i,1)+5,means_data(i,2), labels_en{i},'Rotation',0,'FontSize',12)
    end
end
hold on

Fit = polyfit(means_data(:,1),means_data(:,2),1); % x = x data, y = y data, 1 = order of the polynomial.            
ttt=sort(means_data(:,1));  
yhat2=polyval(Fit,means_data(:,1)); 
yhat=polyval(Fit,ttt);                                                     
plot([ttt(1) ttt(end)],[yhat(1) yhat(end)],'LineWidth',2,'Color',[0 0 0]); 
hold on
map = cbrewer('seq','Reds',27)
map = flipud(map)
error = (yhat2 - means_data(:,2)).^2
[errors idx] = sort(error, 'descend')

% lines to see easily the distance to the least squares line
%for ii=1:27
%    plot([sums(ii,1) sums(ii,1)], [yhat2(ii) sums(ii,2)], 'LineWidth',2,'Color', map(idx==ii,:))
%end
hold on
axis square

label_h = ylabel('Mean paint value, atheists (n = 39)', 'FontSize', 20)
xlabel('Mean paint value, religious (n = 34)', 'FontSize', 20)
xlim([0 800])

set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 .75 1]);
print(gcf, 'figs/supplement/religious_embody.png', '-dpng', '-r500')
