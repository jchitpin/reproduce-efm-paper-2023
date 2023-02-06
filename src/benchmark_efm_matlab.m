%% DESCRIPTION
% This script does the following:
% (1) Simulate various unimolecular, closed-loop, fully-connected networks.
% (2) Benchmark EFM enumeration using FluxModeCalculator.
% (3) Export simulated stoichiometry matrices and a random set of steady
%     state fluxes for benchmarking with the Julia code.
% (4) Export the data for plotting.

%% USER PARAMETERS
% Set working directory of this script. Example:
cd('/home/jchitpin/Documents/PhD/Projects/reproduce-efm-paper-2023/src/')
addpath('functions') % do not change

%% FluxModeCalculator - Ensure this MATLAB package is installed and test it with the following example
clc
clear

% Generate random stoichiometric matrix
n=[10 20];
S=[diag(floor(rand(n(1),1)*3)+1),round(0.7*randn(n(1),n(2)-n(1)))];
A=1;
while rank(A)<n(1)
    A=round(randn(n(1)*1.5,n(1))*0.3);
end
S=A*S;
rev=rand(1,n(2))>0.6;

% Show help of both functions
help calculate_flux_modes
help bin2num_flux_modes

% Calculate binary EFMs
[efm_bin,S_unc,id_unc,T,stats]=calculate_flux_modes(S,rev);

% Calculate coefficients of binary EFMs and check for consistency
[efm,err]=bin2num_flux_modes(efm_bin,S);

%% Import BMC computational sphingolipid model exported from text file (using MATLAB import data toolstrip)
filename = '../data/stoich-corrected.csv';
delimiter = ',';
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);
stoich = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID dataArray ans;

%% Compute EFMs using FluxModeCalculator
% NOTE: sometimes FluxModeCalculator fails to compute EFMs for no reason; running the example script above solves this bug
samples = 100;
srev = false(1, size(stoich,2));
elapsed_sphingo = zeros(samples, 1);
for i=1:length(elapsed_sphingo)
    tic
    [efm_bin,S_unc,id_unc,T,stats]=calculate_flux_modes(stoich, srev);
    [efm,err]=bin2num_flux_modes(efm_bin, stoich);
    elapsed_sphingo(i) = toc;
end

%% Export EFMs to CSV
mean(elapsed_sphingo) % mean run time in seconds
writetable(...
    table(...
    mean(elapsed_sphingo),...
    'VariableNames',...
    {join(['Mean_time_in_seconds_after_', num2str(samples), '_samples'])})...
    , '../data/time-sphingo-fmc.csv',...
    'Delimiter', ','...
)

%% Randomly simulate fully-connected, unimolecular networks, enumerate binary EFMs, and save output
% Simulation code is far from perfect but the error checking is stringent.
% With the following rng(2), n = [20 35] for some reason fails to produce a
% fully-connected network with an EFM spanning each metabolite. For that
% reason, that network cannot be benchmarked with the Julia code. See line
% 122 which is commented out but verifies this problem.
rng(2);

n = [...
    [20 30]
    [20 35]
    [20 40]
    [20 45]
    [20 50]
    [20 55]
    [20 60]
    [20 65]
    [20 70]
    [20 75]
    [20 80]
];

% Benchmark FluxModeCalculator on simulated networks
S = cell(1, size(n,1));
E = zeros(1, size(n,1));
efm = {};
elapsed_time = zeros(samples, size(n,1));
for j = 1:size(elapsed_time,2)
    S{j} = num2cell(simulate_stoich(n(j,1:2)));
    for i =1:size(elapsed_time,1)
        elapsed_time(i,j) = benchmark_efms(cell2mat(S{j}));
    end
    rev = zeros(1,size(cell2mat(S{j}),2));
    [efm_bin,~,~,~,~] = calculate_flux_modes(cell2mat(S{j}),rev);
    E(j) = size(efm_bin,2);
    efm{j} = num2cell(efm_bin);
end

% Export stoichiometry matrices to file
for i = 1:length(S)
    writetable(...
        cell2table(S{i}),...
        join(['../data/stoich-', mat2str(n(i,1)), '-', mat2str(n(i,2)), '.csv']),...
        'Delimiter', ',',...
        'WriteVariableNames', 0....
    );
end

% Export EFM matrices to file
for i = 1:length(efms)
    writetable(...
        cell2table(efm{i}),...
        join(['../data/stoich-', mat2str(n(i,1)), '-', mat2str(n(i,2)), '-efm-matrix.csv']),...
        'Delimiter', ',',...
        'WriteVariableNames', 0....
    );
end

% Export a random set of steady state fluxes for each stoichiometry matrix
% Do this by selecting one of each enumerated EFM and aggregating fluxes
for i = 1:length(S)
    if sum(sum(cell2mat(efm{i}), 2) == 0) ~= 0
        disp(i)
        error('Error. Change rng seed because simulated network is not fully connected.')
    else
        writetable(...
            array2table(sum(cell2mat(efm{i}), 2)),...
            join(['../data/ss-flux-', mat2str(n(i,1)), '-', mat2str(n(i,2)), '.csv']),...
            'Delimiter', ',',...
            'WriteVariableNames', 0....
        );
    end
end

% Export run time
x = elapsed_time;
SEM = std(x) / sqrt(size(x,1)); % standard error
ts = tinv([0.025 0.975], size(x,1)-1); % T-score
CI = [mean(x); mean(x)]' + SEM'*ts; % 95-confidence interval
err = SEM'*ts; % 95-confidence interval for plotting
ttable = table(...
    (1:size(n,1))',...
    mean(x)',...
    err(:,2),...
    E',...
    'VariableNames',...
    {'x', 'y', 'y_ci', 'num_efms'}...
);
writetable(...
    ttable,...
    '../data/benchmark-stoich-fluxmodecalculator.csv',...
    'Delimiter', '\t'...
);

exit;
