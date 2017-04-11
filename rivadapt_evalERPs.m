function rivadapt_evalERPs(expt_rivadapt, ds_to_analyze)


ntotsubj = numel(expt_rivadapt);

if ~exist('ds_to_analyze', 'var')
    ds_to_analyze = 1; % 0-without rereference, 1-with rerereference
end

init_dur = .02148; % value of the initial offset from the light sensor

to_analyze = [expt_rivadapt.to_analyze];
session_dir = expt_rivadapt(1).session_dir;
trig_conds2 = expt_rivadapt(1).trig_conds;
conds_names2 = expt_rivadapt(1).conds_names;
conditions2 = expt_rivadapt(1).conditions;
trial_duration = expt_rivadapt(1).trial_duration;

subj_ii = 0;

for ns = 1:ntotsubj
    if to_analyze(ns)
        % enter the directory
        expt_sub = expt_rivadapt(ns);
        sess_dir = fullfile(session_dir, expt_sub.session_name);
        cd(sess_dir)
        
        subj_ii = subj_ii + 1;

        % load the eeglab dataset file
        if ds_to_analyze
            ds_file = expt_sub.dataset_reref;
        else
            ds_file = expt_sub.dataset;
        end
        data_file = fullfile(sess_dir, 'Data', 'eeglab', ...
            [ds_file '.set']);
        EEG = pop_loadset(data_file);
        fs = EEG.srate;
        chanlabels = {EEG.chanlocs.labels};
        nchans = numel(chanlabels);
        
        ntestconds = numel(trig_conds2);
        if ds_to_analyze
            nt = expt_sub.ntrials_percond_reref;
        else
            nt = expt_sub.ntrials_percond;
        end
        
        epochinds = round([init_dur trial_duration+init_dur]*fs);
                
        epochs_allconds = zeros(ntestconds, 2, numel(chanlabels), diff(epochinds)+1);
        sepoch = size(epochs_allconds);
        
        for ntc = 1:ntestconds
            % for each type of test condition: baseline, fusion adapt, rivalry
            % adapt
            trig_conds = trig_conds2{ntc};
            cond_names = conds_names2{ntc};
            condition = conditions2{ntc};
            
            nconds = numel(trig_conds);
            
            % get number of valid trials
            ntrials = nt(ntc);
            trials = 1:ntrials;
            
            % concatenate all frequencies to analyze
            
            epoch_testcond = zeros(sepoch(2:end));
            
            for cc = 1:nconds
                epochs = Get_Epochs(EEG, trig_conds{cc}, [], epochinds);
                
                epochs = epochs(trials);
                
                epochs3 = permute(reshape(cell2num(epochs), ...
                    [size(epochs{1}, 1), numel(epochs), size(epochs{1}, 2)]), [2 1 3]);
                
                epoch_testcond(cc, :, :) = squeeze(mean(epochs3));
            end
            
            epochs_allconds(ntc, :, :, :) = epoch_testcond;
        end
        
        
        if ~isdir(fullfile(pwd, 'Results')), mkdir('Results'); end
        
        save(fullfile(pwd, 'Results', ['ERP_' ds_file]), 'epochs_allconds', ...
            'condition', 'cond_names', 'ntrials', 'fs', 'chanlabels')
        
    end
end

end