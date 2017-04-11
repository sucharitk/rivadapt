function rivadapt_behavioral(expt_rivadapt, pract)

ntotsubj = numel(expt_rivadapt);
to_analyze = [expt_rivadapt.to_analyze];
session_dir = expt_rivadapt(1).session_dir;
cd(fullfile(session_dir, 'BR32_IOCDAdapt_EEG', 'Results'))

subjnum = 0;
condids = [2 3 0 1]; % baseline, fusion, rivalry

for ns = 1:ntotsubj
    
    if to_analyze(ns)
        subjnum = subjnum + 1;
        % enter the directory
        
        if pract
            load(expt_rivadapt(ns).behav_files_practice);
        else
            load(expt_rivadapt(ns).behav_files);
        end
        psycho = [results.psycho];
        
        ntrialcond = zeros(1, 3);
        nruns = numel(psycho);
        
        if nruns ~= 15
            sprintf('%s %d', expt_rivadapt(ns).behav_files, nruns)
        end
        cond = cell(1, 3);
        for nr = 1:nruns
            trials = psycho(nr).trials;
            cid = condids(trials.adapt_cond+1);
            nts = trials.curTrial;
            
            cond{cid} = ...
                [cond{cid} trials.shimmer_lvl == trials.test_texture];
            ntrialcond(cid) = ntrialcond(cid)+nts;
        end
        perf = zeros(1, 3);
        for cc = 1:3
            perf(cc) = sum(cond{cc})/numel(cond{cc});
        end
        subj_res(subjnum).cond = cond;
        subj_res(subjnum).perf = perf;
    end
end

totperf = [subj_res.perf];
totperf = reshape(totperf', numel(totperf)/subjnum, subjnum)';

if pract
    nanmean(totperf, 2)
else
    if subjnum==1
        plot(totperf)
    else
        errorbar(nanmean(totperf, 1), nanstd(totperf, 1)/sqrt(subjnum))
    end
    anova_rm(totperf(:, [2 1 3]))
    [h, p] = ttest(diff(totperf(:, [2 3]), [], 2))
end
end