function ampfreq = rivadapt_plotERPs_gh(expt_rivadapt, ds_to_analyze, freq_conds, ...
    use_snr, do_plots, starttime)

if ~exist('ds_to_analyze', 'var')
    ds_to_analyze = 1; % 0-without rereference, 1-with rerereference
end
selected_chans = true; % evaluate combined power at selected channels or all channels

load('chanlocs64')

to_analyze = [expt_rivadapt.to_analyze];
session_dir = expt_rivadapt(1).session_dir;
if ~exist('freq_conds', 'var')
    freq_conds = [expt_rivadapt(1).fund_freqs expt_rivadapt(1).inter_freqs];
end
conditions = expt_rivadapt(1).conditions;
post_electrodes = expt_rivadapt(1).post_electrodes;

numf = numel(freq_conds);

fm = 1.4;
sig_freqs = [-.07 +.07];

nse_freqs = [-1.25 +1.25]+sig_freqs;

if ~exist('starttime', 'var'), starttime = [0 1]; end

clip_flag = true;

ntotsubj = numel(expt_rivadapt);

clc

remove_mastoids = true;
snr_cutoff_selectrodes = 2;%.5; % minimum snr to for analyzing the electrodes

nv = 0;
for ns = 1:ntotsubj
    if to_analyze(ns)
        nv = nv + 1;
        % enter the directory
        expt_sub = expt_rivadapt(ns);
        sess_name = expt_sub.session_name;
        sess_dir = fullfile(session_dir, sess_name);
        cd(sess_dir)
        
        if ds_to_analyze
            ds_file = expt_sub.dataset_reref;
        else
            ds_file = expt_sub.dataset;
        end
        
        dat_file = fullfile(pwd, 'Results', ['ERP_' ds_file]);
        
        load(dat_file)
        
        chanlabels64 = {chanlocs.labels};
        chanloc_inds = get_channels_from_labels(chanlabels64, chanlabels);
        
        if remove_mastoids
            mm = strcmp(chanlabels, 'M1') | strcmp(chanlabels, 'M2');
            if any(mm)
                chanlabels = chanlabels(~mm);
                chanrem_inds = get_channels_from_labels(chanlabels64, {'M1', 'M2'});
                chanloc_inds = chanloc_inds & ~chanrem_inds;
                epochs_allconds = epochs_allconds(:, :, ~mm, :);
            end
        end
        chanloc_cond = chanlocs(chanloc_inds);
        
        NFFT = 2^(nextpow2(2*fs+1)+4);
        
        for iocd = 1:2
            for cc = 1:3
                basecond = squeeze(epochs_allconds(cc, iocd, :, :));
                
                [freq_amps, freqs] = FFT_at_freq(basecond, fs, ...
                    NFFT, freq_conds, fm, clip_flag);
                
                for ff = 1:numf
                    frinds = freqs(ff, :)>=freq_conds(ff)+sig_freqs(1) & ...
                        freqs(ff, :)<=freq_conds(ff)+sig_freqs(2);
                    sigval = squeeze(freq_amps(ff, :, frinds));
                    sigff(nv, ff, iocd, cc, :) = max(sigval, [], 2);
                    
                    frinds = freqs(ff, :)<freq_conds(ff)+nse_freqs(1) | ...
                        freqs(ff, :)>freq_conds(ff)+nse_freqs(2);
                    nseff = squeeze(freq_amps(ff, :, frinds));
                    
                    snrff(nv, ff, iocd, cc, :) = squeeze(sigff(nv, ff, iocd, cc, :))./mean(nseff, 2);
                    
                    if iocd==1 && cc==1
                        chanloc_cond_labels = {chanloc_cond.labels};
                        if selected_chans
                            chanloc_postinds = get_channels_from_labels(...
                                chanloc_cond_labels, post_electrodes);
                        else
                            chanloc_postinds = get_channels_from_labels(...
                                chanloc_cond_labels, {chanlocs.labels});
                        end
                        chanloc_postinds = find(chanloc_postinds);
                        
                        sel_snr = squeeze(snrff(nv, ff, iocd, cc, chanloc_postinds));
                        %                     selectrodes = chanloc_postinds(sel_snr>snr_cutoff_selectrodes);
                        selectrodes = chanloc_postinds;
                        if isempty(selectrodes)
                            selectrodes = strcmp({chanloc_cond.labels}, 'Oz');
                        end
                        
                        subj_erps(nv, ff, :, :, :) = ...
                            squeeze(mean(epochs_allconds(:, :, selectrodes, :), 3));
                        
                        %             save('Results\selectrodes', 'selectrodes')
                    end
                end
            end
        end
        
    end
end

subj_erps = subj_erps(:, :, [1 3 2], :, :);
conditions = conditions([1 3 2]);
startind = round(starttime*fs+1);
% NFFT = 2^10;
spectplotlim = [1.25 15.75];

if any(do_plots == 1)
    % plot the overall spectral response    
    merp = squeeze(mean(mean(mean(subj_erps(:, :, :, :, :), 2), 3), 4));
    merp = merp(:, startind(1):startind(2));
    [ft1, ff1] = AbsFFT(merp, fs);
    
    valf = ff1>spectplotlim(1) & ff1<spectplotlim(2);
    figure, plot(ff1(valf), ft1(:, valf), 'LineWidth', 3)
    title('individual subject frequecny spectrum'), xlabel('frequency (hz)'), ylabel('FFT amplitude')
    figure, shadedErrorBar(ff1(valf), mean(ft1(:, valf), 1), std(ft1(:, valf)/sqrt(nv)))
    title('subject averaged frequecny spectrum'), xlabel('frequency (hz)'), ylabel('FFT amplitude')
    ax = axis;
    ax(1:2) = spectplotlim;
    axis(ax)
    [freq_amps, freqs] = FFT_at_freq(merp, fs, NFFT, freq_conds, fm, clip_flag);
    
    for ff = 1:numf
        frinds = freqs(ff, :)>=freq_conds(ff)+sig_freqs(1) & ...
            freqs(ff, :)<=freq_conds(ff)+sig_freqs(2);
        sigval = squeeze(freq_amps(ff, :, frinds));
        sigavg = max(sigval, [], 2);
        
        frinds = freqs(ff, :)<freq_conds(ff)+nse_freqs(1) | ...
            freqs(ff, :)>freq_conds(ff)+nse_freqs(2);
        nseavg = squeeze(freq_amps(ff, :, frinds));
        snr_subj = sigavg./mean(nseavg, 2);
        sprintf('Avg SNR for %g Hz = %g, sem = %g', freq_conds(ff), mean(snr_subj), std(snr_subj)/sqrt(nv))
    end
    
    figure
    use_snr2 = true;
    iocd=1;
    for ff = 1:numf
        if use_snr2
            toplot = squeeze(mean(mean(mean(snrff(:, ff, iocd, :, :), 1), 2), 4));
            minmax = [1.5 3.5];
        else
            toplot = squeeze(mean(mean(mean(sigff(:, ff, iocd, :, :), 1), 2), 4));
            minmax = [0 1];
        end
        subplot(1, 3, ff)
        topoplot(toplot, chanloc_cond, ...
            'maplimits', minmax);
        title(num2str(freq_conds(ff)))
        colorbar;
        
    end
    
    figure % plot the difference of iocd topographies for reviewer's satisfaction
    for cc=1:3
        subplot(1, 3, cc)
        use_snr2 = true;
        for ff = numf:numf
            if ff<3, f2=[1 2]; else f2 =ff; end
            if use_snr2
                toplot = squeeze(diff(mean(mean(mean(snrff(:, f2, :, cc, :), 1), ...
                    2), 4), [], 3));
                minmax = [-.4 .7];
            else
                toplot = squeeze(diff(mean(mean(mean(sigff(:, f2, :, cc, :), 1), ...
                    2), 4), [], 3));
                minmax = [0 1];
            end
            
            topoplot(toplot, chanloc_cond, ...
                'maplimits', minmax);
            title(num2str(freq_conds(ff)))
            colorbar;
        end
    end
end


%% For figures 3 and 4 showing SSVEP adaptation
NFFT = 2^14;

if any(do_plots==2)
    clear sigff nseff
    for nc1 = 1:3
        for nc2 = 1:2
            for ff = 1:numf
                epoch = squeeze(subj_erps(:, 3, nc1, nc2, :));
                if nv==1, epoch = epoch'; end
                [ft1, ff1] = FFT_at_freq(epoch(:, startind(1):startind(2)), ...
                    
            %%% method 1 for SNR
            frinds = ff1>=freq_conds(ff)+sig_freqs(1) & ...
                ff1<=freq_conds(ff)+sig_freqs(2);
            
            sigval = squeeze(ft1(1, :, frinds));
            sigff = max(sigval, [], 2);
            
            frinds2 = ff1<freq_conds(ff)+nse_freqs(1) | ...
                ff1>freq_conds(ff)+nse_freqs(2);
            nseff = squeeze(ft1(:, :, frinds2));
            nseff = mean(nseff, 2);
            
            
            if use_snr
                ampfreq(ff, nc1, nc2, :) = sigff./nseff;
            else
                ampfreq(ff, nc1, nc2, :) = sigff;
            end
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % average in two panels for intermods and funds
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fr_cols = 'brgmcykbrgmcykbrgmcyk';
    figure('Name', 'Combined conditions, fund and intermods');
    switch use_snr
        case true
            if starttime(2)==2
                ax = [-.1 .3 1 4];
            else
                ax = [-.1 .3 1 2.5];
            end
        case false
            ax = [-.1 .3 .4 1.6];
    end
    
    clear mampfreq2 mampfreq1
    for nc1 = 1:3
        subplot(121)
        hold on
        mampfreq2(nc1, :, :) = squeeze(squeeze(ampfreq(3, nc1, :, :)));
        mm = squeeze(mean(mampfreq2(nc1, :, :), 3));
        ss = std(squeeze(diff(mampfreq2(nc1, :, :), [], 2))); ss = [ss; ss];
        errorbar([0 .2], mm, ss/sqrt(nv), fr_cols(nc1), 'LineWidth', 5)
        axis(ax)
        xlabel('IOCD'), ylabel('SNR')
        title('intermodulation')
        
        subplot(122)
        hold on
        mampfreq1(nc1, :, :) = squeeze(mean(squeeze(ampfreq([1 2], nc1, :, :)), 1));
        mm = squeeze(mean(mampfreq1(nc1, :, :), 3));
        ss = std(squeeze(diff(mampfreq1(nc1, :, :), [], 2))); ss = [ss; ss];
        errorbar([0 .2], mm, ss/sqrt(nv), fr_cols(nc1), 'LineWidth', 5)
        axis(ax)
        xlabel('IOCD'), ylabel('SNR')
        title('fundamental')
        
        % do r-manova for baseline condition
        iocdfact = repmat([0 1], 1, 2*nv);
        freqfact = [zeros(2*nv, 1); ones(2*nv, 1)];
        subjfact = repmat(1:nv, 2, 2);
        rm_anova2([mampfreq1(:); mampfreq2(:)], subjfact(:), iocdfact', freqfact, ...
            {'iocd', 'freqtype'});
        
        [~, p1, ~, tstat1] = ttest(squeeze(diff(mampfreq1(nc1, :, :), [], 2)));
        %     disp(tstat1)
        df_im(nc1, :) = squeeze(diff(mampfreq2(nc1, :, :), [], 2));
        df_fund(nc1, :) = squeeze(diff(mampfreq1(nc1, :, :), [], 2));
        [~, p2, ~, tstat2] = ttest(df_im(nc1, :));
        %     disp(tstat2)
        [~, p3] = ttest(squeeze(diff(mampfreq2(nc1, :, :)./mampfreq1(nc1, :, :), [], 2)));
        sprintf('%s: ttest pvals for funds=%g, intermods=%g, ratio=%g, intermod_tstat=%g', ...
            conditions{nc1}, p1, p2, p3, tstat2.tstat)
    end
    legend({'baseline', 'conflict adaptation', 'matched adaptation'})
    anova_rm(df_im');
    anova_rm(df_fund');
    
    % intermodulation frequency - comparing
    df1 = squeeze(diff(diff(ampfreq(3, [1 2], :, :), [], 3), [], 2));
    df2 = squeeze(diff(diff(ampfreq(3, [3 2], :, :), [], 3), [], 2));
    df3 = squeeze(diff(diff(ampfreq(3, [1 3], :, :), [], 3), [], 2));
    [~, p1, ~, tstat1] = ttest(df1);
    [~, p2, ~, tstat2] = ttest(df2);
    [~, p3, ~, tstat3] = ttest(df3);
    sprintf('intermod comparison of diff scores: \nbas>riv t=%g p=%g, fus>riv t=%g p=%g, bas>fus t=%g p=%g', ...
        tstat1.tstat, p1, tstat2.tstat, p2, tstat3.tstat, p3)
    
    
    [~, p1, ~, tstat1] = ttest(squeeze(sum(diff(ampfreq(3, [1 2], :, :), [], 2), 3)));
    [~, p2, ~, tstat2] = ttest(squeeze(sum(diff(ampfreq(3, [3 2], :, :), [], 2), 3)));
    [~, p3, ~, tstat3] = ttest(squeeze(sum(diff(ampfreq(3, [1 3], :, :), [], 2), 3)));
    sprintf('intermod comparison of sum scores: \nbas>riv t=%g p=%g, fus>riv t=%g p=%g, bas>fus t=%g p=%g', ...
        tstat1.tstat, p1, tstat2.tstat, p2, tstat3.tstat, p3)
    
    % fundamental frequency - comparing the differences between three
    % adaptation conditions of the difference between 0 and 0.2 - to see if
    % there was a difference in how much they dropped
    [~, p1, ~, tstat1] = ttest(squeeze(diff(diff(sum(ampfreq([1 2], [1 2], :, :), 1), [], 3), [], 2)));
    [~, p2, ~, tstat2] = ttest(squeeze(diff(diff(sum(ampfreq([1 2], [3 2], :, :), 1), [], 3), [], 2)));
    [~, p3, ~, tstat3] = ttest(squeeze(diff(diff(sum(ampfreq([1 2], [1 3], :, :), 1), [], 3), [], 2)));
    sprintf('fund comparison of diff scores: \nbas>riv t=%g p=%g, fus>riv t=%g p=%g, bas>fus t=%g p=%g', ...
        tstat1.tstat, p1, tstat2.tstat, p2, tstat3.tstat, p3)
    
    %% plot the difference and mean scores
    clear mc
    mc(1, :, :) = squeeze(diff(squeeze(mean(ampfreq([1 2], :, :, :), 1)), 1, 2));
    mc(2, :, :) = squeeze(diff(squeeze(mean(ampfreq(3, :, :, :), 1)), 1, 2));
    % mc now contains 2x3xn - freq_type x cond x subj
    xx = [[1; 2; 3]-.15, [1; 2; 3]+.15];
    mm1 = mean(mc, 3)'; ss1 = std(mc, [], 3)'/sqrt(nv);
    
    nel = numel(mc);
    freq_type = repmat([1; 2], nel/2, 1);
    adapt_type = repmat([1; 1; 2; 2; 3; 3], nel/6, 1);
    subj_num = repmat(1:nv, nel/nv, 1); subj_num = subj_num(:);
    sprintf('2x3 rmanova between frequency type and diff scores for 3 conds')
    rm_anova2(mc(:), subj_num, freq_type, adapt_type, {'freq_type', 'adapt_type'})
    
    mc2 = mc(:, 1:2, :);
    nel = numel(mc2);
    freq_type = repmat([1; 2], nel/2, 1);
    adapt_type = repmat([1; 1; 2; 2], nel/4, 1);
    subj_num = repmat(1:nv, nel/nv, 1); subj_num = subj_num(:);
    sprintf('2x2 rmanova between frequency type and adaptation type (only fused and rivalry) of diff scores')
    rm_anova2(mc2(:), subj_num, freq_type, adapt_type, {'freq_type', 'adapt_type'})
        
    figure('Name', 'Difference scores for the two types of frequencies')
    if use_snr
        ax = [0 4 -0.6 .4];
    else
        ax = [0 4 -0.18 0.6];
    end
    subplot(121)
    errorbar(mm1(:, 2), ss1(:, 2), 'LineWidth', 5)
    % anova_rm(squeeze(mc(2, : , :)));
    hold on
    bar(mm1(:, 2));
    axis(ax)
    subplot(122)
    errorbar(mm1(:, 1), ss1(:, 1), 'LineWidth', 5)
    hold on
    bar(mm1(:, 1))
    axis(ax)
end

%%%%%%%%%%%
%% sliding window analysis
%%%%%%%%
if any(do_plots==4)
    fsamp = 200; tsamp = 1/fsamp;
    starttime = [0:tsamp:.8; .8:tsamp:1.6]';
    ntp = size(starttime, 1);
    
    clear mm ss
    for st = 1:ntp
        startind = round(starttime(st, :)*fs+1);
        
        clear sigff nseff
        for nc1 = 1:3
            for nc2 = 1:2
                for ff = 1:numf
                    epoch = squeeze(subj_erps(:, 3, nc1, nc2, :));
                    if nv==1, epoch = epoch'; end
                    [ft1, ff1] = FFT_at_freq(epoch(:, startind(1):startind(2)), ...
                        fs, NFFT, freq_conds(ff), fm, true);
                    
                    frinds = ff1>=freq_conds(ff)+sig_freqs(1) & ...
                        ff1<=freq_conds(ff)+sig_freqs(2);
                    
                    sigval = squeeze(ft1(1, :, frinds));
                    sigff = max(sigval, [], 2);
                    
                    frinds2 = ff1<freq_conds(ff)+nse_freqs(1) | ...
                        ff1>freq_conds(ff)+nse_freqs(2);
                    nseff = squeeze(ft1(:, :, frinds2));
                    nseff = mean(nseff, 2);
                    
                    
                    if use_snr
                        ampfreq(ff, nc1, nc2, :) = sigff./nseff;
                    else
                        ampfreq(ff, nc1, nc2, :) = sigff;
                    end
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % average in two panels for intermods and funds
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for nc1 = 1:3
            % funds
            mampfreq1 = squeeze(mean(squeeze(ampfreq([1 2], nc1, :, :)), 1));
            mm(1, st, nc1, :) = mean(mampfreq1, 2);
            ss(1, st, nc1) = std(diff(mampfreq1));
            
            % intermods
            mampfreq2 = squeeze(squeeze(ampfreq(3, nc1, :, :)));
            mm(2, st, nc1, :) = mean(mampfreq2, 2);
            ss(2, st, nc1) = std(diff(mampfreq2));
            
            % ratio of intermod to fund
            mm(3, st, nc1, :) = mean(mampfreq2./mampfreq1, 2);
            ss(3, st, nc1) = std(diff(mampfreq2./mampfreq1));
            
            allsubj(2, st, nc1, :, :) = mampfreq2;
            allsubj(3, st, nc1, :, :) = mampfreq2./mampfreq1;
        end
    end
    fr_cols = 'brgm';
    
    freq_labels = {'fundamental', 'intermodulation', 'ratio of inter/fund'};
    tx = mean(starttime, 2);
    for ff = 1:3
        figure('Name', freq_labels{ff}), hold on
        %             subplot(1,3 , ff), hold on
        for nc1 = 1:3
            shadedErrorBar(tx, squeeze(diff(mm(ff, :, nc1, :), [], 4)), ...
                ss(ff, :, nc1)/sqrt(nv), fr_cols(nc1), 1)
        end
        ax = axis; ax(2) = tx(end);
        ax(3:4) = [-2 2];
        axis(ax)
        title(freq_labels{ff})
        xlabel('center of 0.8 sec window (s)'), ylabel('Difference scores between IOCD 0.2 and 0.0')
    end
    
end


end