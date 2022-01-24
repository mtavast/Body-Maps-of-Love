%% DBSCAN sanity check
% For searching the DBSCAN solutions, we used the same parameters as
% Nummenmaa et al, as we used the same task. 
% Here just a double check, that the chosen epsilon range
% makes sense in this study as well. 

clc; clear all; close all;
addpath('external/knee_pt/')

load('output/sim/sim_mds.mat');

m  = find(triu(ones(27),1));
v  = mean_data_D(m).';   % all the pairwise distances

v_sorted = sort(v);      
plot(v_sorted);

kneepoint = knee_pt(v_sorted); % function that finds the location of a knee in the curve, the value close to our range (0.35)

disp(['Value of the kneepoint ' num2str(v_sorted(kneepoint))])