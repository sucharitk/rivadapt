function rivadapt_extractadaptationepochs(expt_rivadapt, ds_to_analyze, epoch_type)

if ~exist('ds_to_analyze', 'var')
    ds_to_analyze = 1; % 0-without rereference, 1-with rerereference
end

to_analyze = [expt_rivadapt.to_analyze];
session_dir = expt_rivadapt(1).session_dir;

switch epoch_type
    case 'initadapt'
        adapt_trigs = {'10', '11'};
    case 'topup'
        adapt_trigs = {'20', '30'; '21', '31'};
end

ntotsubj = numel(expt_rivadapt);
nv = 0;
for ns = 1:ntotsubj
    if to_analyze(ns)
        nv = nv + 1;
        % enter the directory
        expt_sub = expt_rivadapt(ns);
        sess_name = expt_sub.session_name;
        sess_dir = fullfile(session_dir, sess_name);
        cd(sess_dir)
        
        switch ds_to_analyze
            case 0
                ds_file = expt_sub.dataset;
            case 1
                ds_file = expt_sub.dataset_reref;
            case 2
                ds_file = expt_sub.dataset_withart;
        end
        
        fname = fullfile('Data/eeglab', [ds_file '.set']);
        EEG = pop_loadset(fname);
        
        switch epoch_type
            case 'initadapt'
                for nat = 1:numel(adapt_trigs)
                    epochs{nat} = Get_Epochs(EEG, adapt_trigs{nat}, [], [], true);
                end
            case 'topup'
                for nat = 1:numel(adapt_trigs)
                    epochs{nat} = Get_Epochs(EEG, adapt_trigs{nat}, [], [], true);
                end
        end
               
        fname_out = fullfile('Results\', ['adaptepochs_' ds_file]);
        chanlocs = EEG.chanlocs;
        fs = EEG.srate;
        
        save(fname_out, 'epochs', 'chanlocs', 'fs');
        
    end
end

end