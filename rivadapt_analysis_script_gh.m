%% adapt rivalry full analysis script

if ~ispc
    addpath ~/Desktop/Local/Analysis_Such/EEG/
    addpath ~/Desktop/Local/NeuroMatlabToolboxes/Paths/
    AddEEGPaths
    addpath ~/google_drive/Projects/Experiments/Rivalry_Adaptation_SSVEP/analysis/
    eeglab
    
else
    addpath C:\Analysis\Analysis_Such\EEG
    addpath C:\Analysis\NeuroMatlabToolboxes\Paths
    AddEEGPaths
    addpath 'C:\OneDrive - UC Davis\Projects\Experiments\Rivalry_Adaptation_SSVEP\analysis'
    eeglab
    
end

%% create experiment

clear all

expt_rivadapt = rivadapt_experiment_params;

%% extract the trials for analysis and save in results folder
 
ds_to_analyze = 1; % 0-without rereference, 1-with rerereference
rivadapt_evalERPs(expt_rivadapt, ds_to_analyze)

%% actual full analysis of data

freq_conds = [7.2 12 4.8];
ds_to_analyze = 1; % 0-without rereference, 1-with rerereference
use_snr = true; % true- SNR, false- amp
do_plots = [2];
mean_period = [0 .8];
if all(do_plots<2)
    rivadapt_plotERPs(expt_rivadapt, ds_to_analyze, freq_conds, ...
        use_snr, do_plots, mean_period);
else
    ampfreq = rivadapt_plotERPs(expt_rivadapt, ds_to_analyze, freq_conds, ...
        use_snr, do_plots, mean_period);
end

%% behavioral analysis
pract = false; % true: for practice runs, false: for exp runs
rivadapt_behavioral(expt_rivadapt, pract)

%% get the initial 60 s adaptation epochs to compare fused and rivaling gratings

ds_to_analyze = 1; % 0-without rereference, 1-with rerereference, 2-with artifacts w/o ref
epoch_type = 'initadapt'; % initadapt or topup
rivadapt_extractadaptationepochs(expt_rivadapt,ds_to_analyze, epoch_type)
