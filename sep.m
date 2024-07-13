%% 0. seperate 条件
% 42 hz总是在左边，44 hz总是在右边。
% trigger：10代表左边先开始闪，20代表右边先开始闪。
% PO3:45; PO4:46

%% 1. SSVEP topography，左右半球各选一个电极,先计算ERP，再做fft
load("E:\Data\Kongqing\Exp1\single\sep.mat");
cfg = [];
cfg.demean  = 'yes';
cfg.baselinewindow = [-0.2, 0];
for i=1:length(eegdata)
    eegdata{1,i} = ft_preprocessing(cfg, eegdata{1,i});
end

% 分开左边先闪和右边先闪，各自计算erp
subj=1:21; 
cfg = [];
for i = 1:length(subj)
    cfg.trials  = find(eegdata{1, subj(i)}.trialinfo == 10);
    erp_l{1, i}    = ft_timelockanalysis(cfg,eegdata{1, subj(i)});
    cfg.trials  = find(eegdata{1, subj(i)}.trialinfo == 20);
    erp_r{1, i}    = ft_timelockanalysis(cfg,eegdata{1, subj(i)});
    i
end

% grandaverage
cfg = [];
erp_l_ave = ft_timelockgrandaverage(cfg,erp_l{1,1},erp_l{1,2},erp_l{1,3},erp_l{1,4},erp_l{1,5},erp_l{1,6},...
    erp_l{1,7},erp_l{1,8},erp_l{1,9},erp_l{1,10},erp_l{1,11},erp_l{1,12},erp_l{1,13},erp_l{1,14},erp_l{1,15},...
    erp_l{1,16},erp_l{1,17},erp_l{1,18},erp_l{1,19},erp_l{1,20},erp_l{1,21});
erp_r_ave = ft_timelockgrandaverage(cfg,erp_r{1,1},erp_r{1,2},erp_r{1,3},erp_r{1,4},erp_r{1,5},erp_r{1,6},...
    erp_r{1,7},erp_r{1,8},erp_r{1,9},erp_r{1,10},erp_r{1,11},erp_r{1,12},erp_r{1,13},erp_r{1,14},erp_r{1,15},...
    erp_r{1,16},erp_r{1,17},erp_r{1,18},erp_r{1,19},erp_r{1,20},erp_r{1,21});
% plot
cfg = [];
cfg.showlabels = 'yes';
cfg.fontsize = 6;
cfg.layout = 'easycapM1.mat';
% cfg.ylim = [-3e-13 3e-13];
ft_multiplotER(cfg, erp_l{1, 2}, erp_r{1, 2});
% ft_multiplotER(cfg, erp_l_ave, erp_r_ave);

cfg = [];
cfg.xlim = [-0.2 1.7];
% cfg.ylim = [-1e-13 3e-13];
cfg.channel = 'PO3';
ft_singleplotER(cfg, erp_l_ave, erp_r_ave);

cfg = [];
cfg.xlim = [-0.2 1.7];
% cfg.ylim = [-1e-13 3e-13];
cfg.channel = 'PO4';
ft_singleplotER(cfg, erp_l_ave, erp_r_ave);

% 导出原始数据
po3_erp_l=erp_l_ave.avg(45,:);
po3_erp_r=erp_r_ave.avg(45,:);
po4_erp_l=erp_l_ave.avg(46,:);
po4_erp_r=erp_r_ave.avg(46,:);
t=erp_l_ave.time;
save('erp_single.mat','po3_erp_l','po3_erp_r','po4_erp_l','po4_erp_r','t');

% 分开基线前和基线后
cfg = [];
cfg.toilim = [-0.2, 0];
for i=1:length(eegdata)
    eegdata_bs{1,i} = ft_redefinetrial(cfg, eegdata{1,i});
end
cfg = [];
cfg.toilim = [0, 1.7];
for i=1:length(eegdata)
    eegdata_er{1,i} = ft_redefinetrial(cfg, eegdata{1,i});
end

% 分开左边先闪和右边先闪，各自计算erp
subj=1:21;
cfg = [];
for i = 1:length(subj)
    cfg.trials  = find(eegdata_bs{1, subj(i)}.trialinfo == 10);
    erp_bs_l{1, i}    = ft_timelockanalysis(cfg,eegdata_bs{1, subj(i)});
    cfg.trials  = find(eegdata_bs{1, subj(i)}.trialinfo == 20);
    erp_bs_r{1, i}    = ft_timelockanalysis(cfg,eegdata_bs{1, subj(i)});
    cfg.trials  = find(eegdata_er{1, subj(i)}.trialinfo == 10);
    erp_er_l{1, i}    = ft_timelockanalysis(cfg,eegdata_er{1, subj(i)});
    cfg.trials  = find(eegdata_er{1, subj(i)}.trialinfo == 20);
    erp_er_r{1, i}    = ft_timelockanalysis(cfg,eegdata_er{1, subj(i)});    
    i
end

% FFT analysis
cfg = [];
cfg.output    = 'pow';
cfg.channel   = 'all';
cfg.method    = 'mtmfft';
cfg.pad       = 10;
cfg.foi       = 36:1:48;
cfg.taper     = 'hanning';
for i=1:length(subj)
    erp_bs_l_fft{1, i} = ft_freqanalysis(cfg, erp_bs_l{1, i});
    erp_bs_r_fft{1, i} = ft_freqanalysis(cfg, erp_bs_r{1, i});
    erp_er_l_fft{1, i} = ft_freqanalysis(cfg, erp_er_l{1, i});
    erp_er_r_fft{1, i} = ft_freqanalysis(cfg, erp_er_r{1, i});    
    i
end

% convert to zscore
erp_er_l_fft_z=erp_er_l_fft;
erp_er_r_fft_z=erp_er_r_fft;
for i = 1:length(subj)
    erp_er_l_fft_z{1, i}.powspctrm = zscore(erp_er_l_fft{1, i}.powspctrm,0,2);
    erp_er_r_fft_z{1, i}.powspctrm = zscore(erp_er_r_fft{1, i}.powspctrm,0,2);
end

erp_er_l_fft_z_ave = erp_er_l_fft_z{1, 1};
erp_er_r_fft_z_ave = erp_er_r_fft_z{1, 1};
spct_er_l = 0;
spct_er_r = 0;
for i = 1:length(subj)
    spct_er_l = spct_er_l + erp_er_l_fft_z{1, i}.powspctrm;
    spct_er_r = spct_er_r + erp_er_r_fft_z{1, i}.powspctrm;   
end
erp_er_l_fft_z_ave.powspctrm = spct_er_l/length(subj);
erp_er_r_fft_z_ave.powspctrm = spct_er_r/length(subj);

% plot
cfg = [];
cfg.parameter    = 'powspctrm';
cfg.showlabels   = 'yes';
cfg.layout       = 'easycapM1.mat';
figure; 
ft_multiplotER(cfg, erp_er_l_fft_z_ave);

figure; 
ft_multiplotER(cfg, erp_er_r_fft_z_ave);

% topo plot
cfg = [];
cfg.parameter    = 'powspctrm';
cfg.xlim         = [42, 42];
cfg.zlim         = [-0.4, 2.4];
cfg.colormap     = brewermap([],'*RdBu');
cfg.colorbar     = 'yes';
cfg.layout       = 'easycapM1.mat';

figure; 
ft_topoplotER(cfg, erp_er_l_fft_z_ave);
figure; 
ft_topoplotER(cfg, erp_er_r_fft_z_ave);

cfg.xlim         = [44, 44];
cfg.zlim         = [-0.4, 2];
figure; 
ft_topoplotER(cfg, erp_er_l_fft_z_ave);
figure; 
ft_topoplotER(cfg, erp_er_r_fft_z_ave);

% 导出频谱
for i = 1:length(subj)
    spec_po3_l(i,:)=erp_er_l_fft_z{1, i}.powspctrm(45,:);
    spec_po3_r(i,:)=erp_er_r_fft_z{1, i}.powspctrm(45,:);
    spec_po4_l(i,:)=erp_er_l_fft_z{1, i}.powspctrm(46,:);
    spec_po4_r(i,:)=erp_er_r_fft_z{1, i}.powspctrm(46,:);    
end
f=erp_er_l_fft_z_ave.freq;
save('spect_single.mat','spec_po3_l','spec_po3_r','spec_po4_l','spec_po4_r','f');

% 提取PO3和PO4的数值来统计
% 44 hz
for i = 1:length(subj)
    fft_po3_bs_l(i)=erp_bs_l_fft{1, i}.powspctrm(45,9);
    fft_po3_bs_r(i)=erp_bs_r_fft{1, i}.powspctrm(45,9);
    fft_po3_er_l(i)=erp_er_l_fft{1, i}.powspctrm(45,9);
    fft_po3_er_r(i)=erp_er_r_fft{1, i}.powspctrm(45,9);
end
fft_po3_44_bs=(fft_po3_bs_l+fft_po3_bs_r)/2;
fft_po3_44_er=(fft_po3_er_l+fft_po3_er_r)/2;
% 42 hz
for i = 1:length(subj)
    fft_po4_bs_l(i)=erp_bs_l_fft{1, i}.powspctrm(46,7);
    fft_po4_bs_r(i)=erp_bs_r_fft{1, i}.powspctrm(46,7);
    fft_po4_er_l(i)=erp_er_l_fft{1, i}.powspctrm(46,7);
    fft_po4_er_r(i)=erp_er_r_fft{1, i}.powspctrm(46,7);
end
fft_po4_42_bs=(fft_po4_bs_l+fft_po4_bs_r)/2;
fft_po4_42_er=(fft_po4_er_l+fft_po4_er_r)/2;
save('spect_value_single.mat','fft_po3_44_bs','fft_po3_44_er','fft_po4_42_bs','fft_po4_42_er');

%% 2. 抽取选择的电极上SSVEP的amplitude envelop
subj=1:21; 
% filter
cfg = [];
cfg.bpfilter    = 'yes';
for i = 1:length(subj)
    cfg.bpfreq     = [41.5, 42.5];
    cfg.trials  = find(eegdata_er{1, subj(i)}.trialinfo == 10);
    eeg_42_l{1, i}    = ft_preprocessing(cfg,eegdata_er{1, subj(i)});
    cfg.trials  = find(eegdata_er{1, subj(i)}.trialinfo == 20);
    eeg_42_r{1, i}    = ft_preprocessing(cfg,eegdata_er{1, subj(i)});    
    i
end
for i = 1:length(subj)
    cfg.bpfreq     = [43.5, 44.5];
    cfg.trials  = find(eegdata_er{1, subj(i)}.trialinfo == 10);
    eeg_44_l{1, i}    = ft_preprocessing(cfg,eegdata_er{1, subj(i)});
    cfg.trials  = find(eegdata_er{1, subj(i)}.trialinfo == 20);
    eeg_44_r{1, i}    = ft_preprocessing(cfg,eegdata_er{1, subj(i)});    
    i
end

% envelope
for i = 1:length(subj)
    for j = 1:160
        for m = 1:65
            [u1,l1]=envelope(eeg_42_l{1, i}.trial{1, j}(m,:),13,'peak');
            eeg_42_l{1, i}.trial{1, j}(m,:)=u1;
            [u2,l2]=envelope(eeg_42_r{1, i}.trial{1, j}(m,:),13,'peak');
            eeg_42_r{1, i}.trial{1, j}(m,:)=u2;  
            [u3,l3]=envelope(eeg_44_l{1, i}.trial{1, j}(m,:),13,'peak');
            eeg_44_l{1, i}.trial{1, j}(m,:)=u3;  
            [u4,l4]=envelope(eeg_44_r{1, i}.trial{1, j}(m,:),13,'peak');
            eeg_44_r{1, i}.trial{1, j}(m,:)=u4; 
        end
    end
    i
end

% time-lock analysis
cfg = [];
for i = 1:length(subj)
    hilbert_42_l{1, i}    = ft_timelockanalysis(cfg,eeg_42_l{1, i});
    hilbert_42_r{1, i}    = ft_timelockanalysis(cfg,eeg_42_r{1, i});
    hilbert_44_l{1, i}    = ft_timelockanalysis(cfg,eeg_44_l{1, i});
    hilbert_44_r{1, i}    = ft_timelockanalysis(cfg,eeg_44_r{1, i});    
    i
end

% grandaverage
cfg = [];
hilbert_42_l_ave = ft_timelockgrandaverage(cfg,hilbert_42_l{1,1},hilbert_42_l{1,2},hilbert_42_l{1,3},hilbert_42_l{1,4},hilbert_42_l{1,5},hilbert_42_l{1,6},...
    hilbert_42_l{1,7},hilbert_42_l{1,8},hilbert_42_l{1,9},hilbert_42_l{1,10},hilbert_42_l{1,11},hilbert_42_l{1,12},hilbert_42_l{1,13},hilbert_42_l{1,14},hilbert_42_l{1,15},...
    hilbert_42_l{1,16},hilbert_42_l{1,17},hilbert_42_l{1,18},hilbert_42_l{1,19},hilbert_42_l{1,20},hilbert_42_l{1,21});
hilbert_42_r_ave = ft_timelockgrandaverage(cfg,hilbert_42_r{1,1},hilbert_42_r{1,2},hilbert_42_r{1,3},hilbert_42_r{1,4},hilbert_42_r{1,5},hilbert_42_r{1,6},...
    hilbert_42_r{1,7},hilbert_42_r{1,8},hilbert_42_r{1,9},hilbert_42_r{1,10},hilbert_42_r{1,11},hilbert_42_r{1,12},hilbert_42_r{1,13},hilbert_42_r{1,14},hilbert_42_r{1,15},...
    hilbert_42_r{1,16},hilbert_42_r{1,17},hilbert_42_r{1,18},hilbert_42_r{1,19},hilbert_42_r{1,20},hilbert_42_r{1,21});
hilbert_44_l_ave = ft_timelockgrandaverage(cfg,hilbert_44_l{1,1},hilbert_44_l{1,2},hilbert_44_l{1,3},hilbert_44_l{1,4},hilbert_44_l{1,5},hilbert_44_l{1,6},...
    hilbert_44_l{1,7},hilbert_44_l{1,8},hilbert_44_l{1,9},hilbert_44_l{1,10},hilbert_44_l{1,11},hilbert_44_l{1,12},hilbert_44_l{1,13},hilbert_44_l{1,14},hilbert_44_l{1,15},...
    hilbert_44_l{1,16},hilbert_44_l{1,17},hilbert_44_l{1,18},hilbert_44_l{1,19},hilbert_44_l{1,20},hilbert_44_l{1,21});
hilbert_44_r_ave = ft_timelockgrandaverage(cfg,hilbert_44_r{1,1},hilbert_44_r{1,2},hilbert_44_r{1,3},hilbert_44_r{1,4},hilbert_44_r{1,5},hilbert_44_r{1,6},...
    hilbert_44_r{1,7},hilbert_44_r{1,8},hilbert_44_r{1,9},hilbert_44_r{1,10},hilbert_44_r{1,11},hilbert_44_r{1,12},hilbert_44_r{1,13},hilbert_44_r{1,14},hilbert_44_r{1,15},...
    hilbert_44_r{1,16},hilbert_44_r{1,17},hilbert_44_r{1,18},hilbert_44_r{1,19},hilbert_44_r{1,20},hilbert_44_r{1,21});

% plot
cfg = [];
cfg.showlabels = 'yes';
cfg.fontsize = 6;
cfg.layout = 'easycapM1.mat';
% cfg.ylim = [-3e-13 3e-13];
ft_multiplotER(cfg, hilbert_42_l{1, 1}, hilbert_42_r{1, 1});

cfg = [];
cfg.xlim = [0 1.7];
% cfg.ylim = [-1e-13 3e-13];
cfg.channel = 'PO3';
ft_singleplotER(cfg, hilbert_42_l_ave, hilbert_42_r_ave);

cfg = [];
cfg.xlim = [0 1.7];
% cfg.ylim = [-1e-13 3e-13];
cfg.channel = 'PO4';
ft_singleplotER(cfg, hilbert_42_l_ave, hilbert_42_r_ave);

%% 3. 两个条件的envelop之间的关系，考察它们之间的相位差
% FFT analysis
cfg = [];
cfg.output    = 'fooof_peaks';
cfg.channel   = 'all';
cfg.method    = 'mtmfft';
cfg.pad       = 10;
cfg.foi       = 2:0.2:10;
cfg.taper     = 'hanning';
for i=1:length(subj)
    eeg_42_l_fft{1, i} = ft_freqanalysis(cfg, eeg_42_l{1, i});
    eeg_42_r_fft{1, i} = ft_freqanalysis(cfg, eeg_42_r{1, i});
    eeg_44_l_fft{1, i} = ft_freqanalysis(cfg, eeg_44_l{1, i});
    eeg_44_r_fft{1, i} = ft_freqanalysis(cfg, eeg_44_r{1, i});    
    i
end
freq=2:0.2:10;
% grand average
eeg_42_l_fft_ave = eeg_42_l_fft{1, 1};
eeg_42_r_fft_ave = eeg_42_r_fft{1, 1};
eeg_44_l_fft_ave = eeg_44_l_fft{1, 1};
eeg_44_r_fft_ave = eeg_44_r_fft{1, 1};
for i = 1:length(subj)
    spct_42_l(i,:,:) = eeg_42_l_fft{1, i}.powspctrm;
    spct_42_r(i,:,:) = eeg_42_r_fft{1, i}.powspctrm;
    spct_44_l(i,:,:) = eeg_44_l_fft{1, i}.powspctrm;
    spct_44_r(i,:,:) = eeg_44_r_fft{1, i}.powspctrm;    
end

eeg_42_l_fft_ave.powspctrm = squeeze(mean(spct_42_l));
eeg_42_r_fft_ave.powspctrm = squeeze(mean(spct_42_r));
eeg_44_l_fft_ave.powspctrm = squeeze(mean(spct_44_l));
eeg_44_r_fft_ave.powspctrm = squeeze(mean(spct_44_r));

% plot
cfg = [];
cfg.parameter    = 'powspctrm';
cfg.showlabels   = 'yes';
cfg.layout       = 'easycapM1.mat';
figure; 
ft_multiplotER(cfg, eeg_42_l_fft_ave, eeg_42_r_fft_ave);

figure; 
ft_multiplotER(cfg, eeg_44_l_fft_ave, eeg_44_r_fft_ave);

% topo plot
eeg_42_l_fft_temp=eeg_42_l_fft_ave;
eeg_42_l_fft_temp.powspctrm=(eeg_42_l_fft_ave.powspctrm+eeg_42_r_fft_ave.powspctrm)/2;
eeg_44_l_fft_temp=eeg_44_l_fft_ave;
eeg_44_l_fft_temp.powspctrm=(eeg_44_l_fft_ave.powspctrm+eeg_44_r_fft_ave.powspctrm)/2;

cfg = [];
cfg.parameter    = 'powspctrm';
cfg.xlim         = [4.4, 4.4];
cfg.zlim         = [6.4, 7.4];
cfg.colormap     = brewermap([],'Reds');
cfg.colorbar     = 'yes';
cfg.layout       = 'easycapM1.mat';

ft_topoplotER(cfg, eeg_42_l_fft_temp);

cfg.xlim         = [6.2, 6.2];
cfg.zlim         = [12.5, 14];
ft_topoplotER(cfg, eeg_44_l_fft_temp);

% permutation test
% FFT analysis
cfg = [];
cfg.output    = 'fooof_peaks';
cfg.channel   = 'all';
cfg.method    = 'mtmfft';
cfg.pad       = 10;
cfg.foi       = 2:0.2:10;
cfg.taper     = 'hanning';

eeg_42_l_temp=eeg_42_l;
eeg_42_r_temp=eeg_42_r;
eeg_44_l_temp=eeg_44_l;
eeg_44_r_temp=eeg_44_r;
parob=parpool;
for n=1:1000
    % shuffle
    parfor i=1:length(subj)
        for j=1:160
            for m=1:65
                eeg_42_l_temp{1, i}.trial{1, j}(m,:)=eeg_42_l{1, i}.trial{1, j}(m,randperm(850));
                eeg_42_r_temp{1, i}.trial{1, j}(m,:)=eeg_42_r{1, i}.trial{1, j}(m,randperm(850));
                eeg_44_l_temp{1, i}.trial{1, j}(m,:)=eeg_44_l{1, i}.trial{1, j}(m,randperm(850));
                eeg_44_r_temp{1, i}.trial{1, j}(m,:)=eeg_44_r{1, i}.trial{1, j}(m,randperm(850));
            end
        end
        eeg_42_l_fft_temp{1, i} = ft_freqanalysis(cfg, eeg_42_l_temp{1, i});
        eeg_42_r_fft_temp{1, i} = ft_freqanalysis(cfg, eeg_42_r_temp{1, i});
        eeg_44_l_fft_temp{1, i} = ft_freqanalysis(cfg, eeg_44_l_temp{1, i});
        eeg_44_r_fft_temp{1, i} = ft_freqanalysis(cfg, eeg_44_r_temp{1, i});
    end

    % grand average
    spct_42_l = 0;
    spct_42_r = 0;
    spct_44_l = 0;
    spct_44_r = 0;
    for i = 1:length(subj)
        spct_42_l = spct_42_l + eeg_42_l_fft_temp{1, i}.powspctrm;
        spct_42_r = spct_42_r + eeg_42_r_fft_temp{1, i}.powspctrm;
        spct_44_l = spct_44_l + eeg_44_l_fft_temp{1, i}.powspctrm;
        spct_44_r = spct_44_r + eeg_44_r_fft_temp{1, i}.powspctrm;
    end

    eeg_42_l_shuf(n,:,:) = spct_42_l/length(subj);
    eeg_42_r_shuf(n,:,:) = spct_42_r/length(subj);
    eeg_44_l_shuf(n,:,:) = spct_44_l/length(subj);
    eeg_44_r_shuf(n,:,:) = spct_44_r/length(subj);
    n
end
delete(parob);
save('env_fft_shuf_single.mat','eeg_42_l_shuf','eeg_42_r_shuf','eeg_44_l_shuf','eeg_44_r_shuf');

% 取95%的阈限
load('env_fft_shuf_single.mat');
eeg_42_l_thre = squeeze(prctile(eeg_42_l_shuf,95,1));
eeg_42_r_thre = squeeze(prctile(eeg_42_r_shuf,95,1));
eeg_44_l_thre = squeeze(prctile(eeg_44_l_shuf,95,1));
eeg_44_r_thre = squeeze(prctile(eeg_44_r_shuf,95,1));

% 导出数据
spct_42_l_po4=squeeze(spct_42_l(:,46,:));
spct_42_r_po4=squeeze(spct_42_r(:,46,:));
spct_44_l_po3=squeeze(spct_44_l(:,45,:));
spct_44_r_po3=squeeze(spct_44_r(:,45,:));
spct_42_l_po4_thre=eeg_42_l_thre(56,:);
spct_42_r_po4_thre=eeg_42_r_thre(56,:);
spct_44_l_po3_thre=eeg_44_l_thre(55,:);
spct_44_r_po3_thre=eeg_44_r_thre(55,:);
save('env_fft_single.mat','freq','spct_42_l_po4','spct_42_r_po4','spct_44_l_po3','spct_44_r_po3',...
    'spct_42_l_po4_thre','spct_42_r_po4_thre','spct_44_l_po3_thre','spct_44_r_po3_thre');

% phase analysis
cfg = [];
cfg.output    = 'fourier';
cfg.channel   = 'all';
cfg.method    = 'mtmfft';
cfg.pad       = 10;
cfg.foi       = 2:0.2:10;
cfg.taper     = 'hanning';
for i=1:length(subj)
    eeg_42_l_phase{1, i} = ft_freqanalysis(cfg, eeg_42_l{1, i});
    eeg_42_r_phase{1, i} = ft_freqanalysis(cfg, eeg_42_r{1, i});
    eeg_44_l_phase{1, i} = ft_freqanalysis(cfg, eeg_44_l{1, i});
    eeg_44_r_phase{1, i} = ft_freqanalysis(cfg, eeg_44_r{1, i});    
    i
end

for i=1:length(subj)
    po4_42_l_phase(i) = circ_mean(angle(eeg_42_l_phase{1, i}.fourierspctrm(:,46,17)));
    po4_42_r_phase(i) = circ_mean(angle(eeg_42_r_phase{1, i}.fourierspctrm(:,46,17)));
    po3_44_l_phase(i) = circ_mean(angle(eeg_44_l_phase{1, i}.fourierspctrm(:,45,17)));
    po3_44_r_phase(i) = circ_mean(angle(eeg_44_r_phase{1, i}.fourierspctrm(:,45,17))); 
    i
end

phase_l_diff =  circ_dist(po4_42_l_phase, po3_44_l_phase);
phase_r_diff =  circ_dist(po3_44_r_phase, po4_42_r_phase);

% plot
figure;
subplot(131);
circ_plot(po4_42_l_phase','hist',[],12,false,true,'linewidth',2,'color','r');
subplot(132);
circ_plot(po3_44_l_phase','hist',[],12,false,true,'linewidth',2,'color','r');
subplot(133);
circ_plot(phase_l_diff','hist',[],12,false,true,'linewidth',2,'color','r');

circ_rad2ang(circ_mean(po4_42_l_phase,[],2));
[pval z] = circ_rtest(po4_42_l_phase); % non-uniform
circ_rad2ang(circ_mean(po3_44_l_phase,[],2));
[pval z] = circ_rtest(po3_44_l_phase); % non-uniform
circ_rad2ang(circ_mean(phase_l_diff,[],2));
[pval z] = circ_rtest(phase_l_diff); % non-uniform

figure;
subplot(131);
circ_plot(po3_44_r_phase','hist',[],12,false,true,'linewidth',2,'color','r');
subplot(132);
circ_plot(po4_42_r_phase','hist',[],12,false,true,'linewidth',2,'color','r');
subplot(133);
circ_plot(phase_r_diff','hist',[],12,false,true,'linewidth',2,'color','r');

circ_rad2ang(circ_mean(po3_44_r_phase,[],2));
[pval z] = circ_rtest(po3_44_r_phase); % non-uniform
circ_rad2ang(circ_mean(po4_42_r_phase,[],2));
[pval z] = circ_rtest(po4_42_r_phase); % non-uniform
circ_rad2ang(circ_mean(phase_r_diff,[],2));
[pval z] = circ_rtest(phase_r_diff); % non-uniform

[pval z] = circ_vtest(phase_l_diff,circ_ang2rad(0)); % difference to 0
[pval z] = circ_vtest(phase_r_diff,circ_ang2rad(0)); % difference to 0

phase_diff_single=circ_mean([phase_l_diff; phase_r_diff],[],1);
save('phase_diff_single.mat','phase_diff_single');

%% 4. Phase-amplitude coupling, SSVEP的amplitude与某个地方脑电的phase有coupling。表现为SSVEP envelope与脑电phase的coherence
% 将提取EEG的频谱，然后跟ssvep的envelope计算coherence
% FFT on eeg
cfg = [];
cfg.output    = 'fourier';
cfg.channel   = 'all';
cfg.method    = 'mtmfft';
cfg.pad       = 10;
cfg.foi       = 2:0.2:10;
cfg.taper     = 'hanning';
for i=1:length(subj)
    cfg.trials  = find(eegdata_er{1, subj(i)}.trialinfo == 10);
    eeg_fft_l{1, i} = ft_freqanalysis(cfg, eegdata_er{1, subj(i)});
    cfg.trials  = find(eegdata_er{1, subj(i)}.trialinfo == 20);
    eeg_fft_r{1, i} = ft_freqanalysis(cfg, eegdata_er{1, subj(i)}); 
    i
end

% 选择电极，替换IO。计算IO与其他电极的coherence
% 左边先开始闪，PO4（46）替换
for i=1:length(subj)
    eeg_fft_l{1, i}.fourierspctrm(:,20,:)=eeg_42_l_phase{1, i}.fourierspctrm(:,46,:);
end

cfg = [];
cfg.method  = 'coh';
cfg.channelcmb = {'IO' 'Fp1'; 'IO' 'Fp2'; 'IO' 'F3'; 'IO' 'F4'; 'IO' 'C3'; 'IO' 'C4';...
    'IO' 'P3'; 'IO' 'P4'; 'IO' 'O1'; 'IO' 'O2'; 'IO' 'F7'; 'IO' 'F8'; 'IO' 'T7';...
    'IO' 'T8'; 'IO' 'P7'; 'IO' 'P8'; 'IO' 'Fz'; 'IO' 'Cz'; 'IO' 'Pz'; 'IO' 'FC1';...
    'IO' 'FC2'; 'IO' 'CP1'; 'IO' 'CP2'; 'IO' 'FC5'; 'IO' 'FC6'; 'IO' 'CP5'; 'IO' 'CP6';...
    'IO' 'FT9'; 'IO' 'FT10'; 'IO' 'TP9';'IO' 'TP10';'IO' 'F1';'IO' 'F2';'IO' 'C1';...
    'IO' 'C2';'IO' 'P1';'IO' 'P2';'IO' 'AF3';'IO' 'AF4';'IO' 'FC3';'IO' 'FC4';...
    'IO' 'CP3';'IO' 'CP4';'IO' 'PO3';'IO' 'PO4';'IO' 'F5';'IO' 'F6';'IO' 'C5';...
    'IO' 'C6';'IO' 'P5';'IO' 'P6';'IO' 'AF7';'IO' 'AF8';'IO' 'FT7';'IO' 'FT8';...
    'IO' 'TP7';'IO' 'TP8';'IO' 'PO7';'IO' 'PO8';'IO' 'Fpz';'IO' 'CPz';'IO' 'POz';...
    'IO' 'Oz';'IO' 'FCz'};
for i=1:length(subj)
     coh_l_po4{1, i}= ft_connectivityanalysis(cfg, eeg_fft_l{1, i});
end

% 左边先开始闪，PO3（45）
for i=1:length(subj)
    eeg_fft_l{1, i}.fourierspctrm(:,20,:)=eeg_44_l_phase{1, i}.fourierspctrm(:,45,:);
end

cfg = [];
cfg.method  = 'coh';
cfg.channelcmb = {'IO' 'Fp1'; 'IO' 'Fp2'; 'IO' 'F3'; 'IO' 'F4'; 'IO' 'C3'; 'IO' 'C4';...
    'IO' 'P3'; 'IO' 'P4'; 'IO' 'O1'; 'IO' 'O2'; 'IO' 'F7'; 'IO' 'F8'; 'IO' 'T7';...
    'IO' 'T8'; 'IO' 'P7'; 'IO' 'P8'; 'IO' 'Fz'; 'IO' 'Cz'; 'IO' 'Pz'; 'IO' 'FC1';...
    'IO' 'FC2'; 'IO' 'CP1'; 'IO' 'CP2'; 'IO' 'FC5'; 'IO' 'FC6'; 'IO' 'CP5'; 'IO' 'CP6';...
    'IO' 'FT9'; 'IO' 'FT10'; 'IO' 'TP9';'IO' 'TP10';'IO' 'F1';'IO' 'F2';'IO' 'C1';...
    'IO' 'C2';'IO' 'P1';'IO' 'P2';'IO' 'AF3';'IO' 'AF4';'IO' 'FC3';'IO' 'FC4';...
    'IO' 'CP3';'IO' 'CP4';'IO' 'PO3';'IO' 'PO4';'IO' 'F5';'IO' 'F6';'IO' 'C5';...
    'IO' 'C6';'IO' 'P5';'IO' 'P6';'IO' 'AF7';'IO' 'AF8';'IO' 'FT7';'IO' 'FT8';...
    'IO' 'TP7';'IO' 'TP8';'IO' 'PO7';'IO' 'PO8';'IO' 'Fpz';'IO' 'CPz';'IO' 'POz';...
    'IO' 'Oz';'IO' 'FCz'};
for i=1:length(subj)
     coh_l_po3{1, i}= ft_connectivityanalysis(cfg, eeg_fft_l{1, i});
end

% 右边先开始闪，PO3（45）
for i=1:length(subj)
    eeg_fft_r{1, i}.fourierspctrm(:,20,:)=eeg_44_r_phase{1, i}.fourierspctrm(:,45,:);
end

cfg = [];
cfg.method  = 'coh';
cfg.channelcmb = {'IO' 'Fp1'; 'IO' 'Fp2'; 'IO' 'F3'; 'IO' 'F4'; 'IO' 'C3'; 'IO' 'C4';...
    'IO' 'P3'; 'IO' 'P4'; 'IO' 'O1'; 'IO' 'O2'; 'IO' 'F7'; 'IO' 'F8'; 'IO' 'T7';...
    'IO' 'T8'; 'IO' 'P7'; 'IO' 'P8'; 'IO' 'Fz'; 'IO' 'Cz'; 'IO' 'Pz'; 'IO' 'FC1';...
    'IO' 'FC2'; 'IO' 'CP1'; 'IO' 'CP2'; 'IO' 'FC5'; 'IO' 'FC6'; 'IO' 'CP5'; 'IO' 'CP6';...
    'IO' 'FT9'; 'IO' 'FT10'; 'IO' 'TP9';'IO' 'TP10';'IO' 'F1';'IO' 'F2';'IO' 'C1';...
    'IO' 'C2';'IO' 'P1';'IO' 'P2';'IO' 'AF3';'IO' 'AF4';'IO' 'FC3';'IO' 'FC4';...
    'IO' 'CP3';'IO' 'CP4';'IO' 'PO3';'IO' 'PO4';'IO' 'F5';'IO' 'F6';'IO' 'C5';...
    'IO' 'C6';'IO' 'P5';'IO' 'P6';'IO' 'AF7';'IO' 'AF8';'IO' 'FT7';'IO' 'FT8';...
    'IO' 'TP7';'IO' 'TP8';'IO' 'PO7';'IO' 'PO8';'IO' 'Fpz';'IO' 'CPz';'IO' 'POz';...
    'IO' 'Oz';'IO' 'FCz'};
for i=1:length(subj)
     coh_r_po3{1, i}= ft_connectivityanalysis(cfg, eeg_fft_r{1, i});
end

% 右边先开始闪，PO4（46）
for i=1:length(subj)
    eeg_fft_r{1, i}.fourierspctrm(:,20,:)=eeg_42_r_phase{1, i}.fourierspctrm(:,46,:);
end

cfg = [];
cfg.method  = 'coh';
cfg.channelcmb = {'IO' 'Fp1'; 'IO' 'Fp2'; 'IO' 'F3'; 'IO' 'F4'; 'IO' 'C3'; 'IO' 'C4';...
    'IO' 'P3'; 'IO' 'P4'; 'IO' 'O1'; 'IO' 'O2'; 'IO' 'F7'; 'IO' 'F8'; 'IO' 'T7';...
    'IO' 'T8'; 'IO' 'P7'; 'IO' 'P8'; 'IO' 'Fz'; 'IO' 'Cz'; 'IO' 'Pz'; 'IO' 'FC1';...
    'IO' 'FC2'; 'IO' 'CP1'; 'IO' 'CP2'; 'IO' 'FC5'; 'IO' 'FC6'; 'IO' 'CP5'; 'IO' 'CP6';...
    'IO' 'FT9'; 'IO' 'FT10'; 'IO' 'TP9';'IO' 'TP10';'IO' 'F1';'IO' 'F2';'IO' 'C1';...
    'IO' 'C2';'IO' 'P1';'IO' 'P2';'IO' 'AF3';'IO' 'AF4';'IO' 'FC3';'IO' 'FC4';...
    'IO' 'CP3';'IO' 'CP4';'IO' 'PO3';'IO' 'PO4';'IO' 'F5';'IO' 'F6';'IO' 'C5';...
    'IO' 'C6';'IO' 'P5';'IO' 'P6';'IO' 'AF7';'IO' 'AF8';'IO' 'FT7';'IO' 'FT8';...
    'IO' 'TP7';'IO' 'TP8';'IO' 'PO7';'IO' 'PO8';'IO' 'Fpz';'IO' 'CPz';'IO' 'POz';...
    'IO' 'Oz';'IO' 'FCz'};
for i=1:length(subj)
     coh_r_po4{1, i}= ft_connectivityanalysis(cfg, eeg_fft_r{1, i});
end

% the coh value of 4 hz
for i=1:length(subj)
    coh_l_po4_val(:,i)=coh_l_po4{1, i}.cohspctrm(:,13);
    coh_l_po3_val(:,i)=coh_l_po3{1, i}.cohspctrm(:,22);
    coh_r_po4_val(:,i)=coh_r_po4{1, i}.cohspctrm(:,13);
    coh_r_po3_val(:,i)=coh_r_po3{1, i}.cohspctrm(:,22);    
end

% topo plot
chann_label=eegdata{1, 1}.label;
chann_label(20)=[];

cfg = [];
cfg.layout = 'easycapM1.mat';
cfg.xlim = 'maxmin';  % time limitation
% cfg.zlim = [0.5 1];
% cfg.style  = 'straight';
cfg.colorbar = 'yes';
cfg.marker = 'on';
cfg.comment = 'xlim';
cfg.commentpos = 'lefttop';

coh_l_po4_plot = {}; 
coh_l_po4_plot.avg = mean(coh_l_po4_val,2); % data definition, channel * time
coh_l_po4_plot.time = 1; % time definition
coh_l_po4_plot.dimord = 'chan_time';
coh_l_po4_plot.label = chann_label;

figure;
ft_topoplotER(cfg, coh_l_po4_plot);
title('coh_l_po4');
axis tight;

coh_l_po3_plot = {}; 
coh_l_po3_plot.avg = mean(coh_l_po3_val,2); % data definition, channel * time
coh_l_po3_plot.time = 1; % time definition
coh_l_po3_plot.dimord = 'chan_time';
coh_l_po3_plot.label = chann_label;

figure;
ft_topoplotER(cfg, coh_l_po3_plot);
title('coh_l_po3');
axis tight;

coh_r_po3_plot = {}; 
coh_r_po3_plot.avg = mean(coh_r_po3_val,2); % data definition, channel * time
coh_r_po3_plot.time = 1; % time definition
coh_r_po3_plot.dimord = 'chan_time';
coh_r_po3_plot.label = chann_label;

figure;
ft_topoplotER(cfg, coh_r_po3_plot);
title('coh_r_po3');
axis tight;

coh_r_po4_plot = {}; 
coh_r_po4_plot.avg = mean(coh_r_po4_val,2); % data definition, channel * time
coh_r_po4_plot.time = 1; % time definition
coh_r_po4_plot.dimord = 'chan_time';
coh_r_po4_plot.label = chann_label;

figure;
ft_topoplotER(cfg, coh_r_po4_plot);
title('coh_r_po3');
axis tight;

% permutation test
eegdata_er_sel=eegdata_er(subj);
parob=parpool;
for n=1:1000
    % shuffle
    eegdata_er_temp=eegdata_er_sel;
    parfor i=1:length(subj)
        for j=1:320
            for m=1:65
                eegdata_er_temp{1, i}.trial{1, j}(m,:)=eegdata_er_sel{1, i}.trial{1, j}(m,randperm(850));
            end
        end
    end

    % FFT on eeg
    cfg = [];
    cfg.output    = 'fourier';
    cfg.channel   = 'all';
    cfg.method    = 'mtmfft';
    cfg.pad       = 10;
    cfg.foi       = 2:0.2:10;
    cfg.taper     = 'hanning';
    for i=1:length(subj)
        cfg.trials  = find(eegdata_er_temp{1, i}.trialinfo == 10);
        eeg_fft_l_temp{1, i} = ft_freqanalysis(cfg, eegdata_er_temp{1, i});
        cfg.trials  = find(eegdata_er_temp{1, i}.trialinfo == 20);
        eeg_fft_r_temp{1, i} = ft_freqanalysis(cfg, eegdata_er_temp{1, i});
    end

    % 选择电极，替换IO。计算IO与其他电极的coherence
    % 左边先开始闪，PO4（46）替换
    for i=1:length(subj)
        eeg_fft_l_temp{1, i}.fourierspctrm(:,20,:)=eeg_42_l_phase{1, i}.fourierspctrm(:,46,:);
    end

    cfg = [];
    cfg.method  = 'coh';
    cfg.channelcmb = {'IO' 'Fp1'; 'IO' 'Fp2'; 'IO' 'F3'; 'IO' 'F4'; 'IO' 'C3'; 'IO' 'C4';...
        'IO' 'P3'; 'IO' 'P4'; 'IO' 'O1'; 'IO' 'O2'; 'IO' 'F7'; 'IO' 'F8'; 'IO' 'T7';...
        'IO' 'T8'; 'IO' 'P7'; 'IO' 'P8'; 'IO' 'Fz'; 'IO' 'Cz'; 'IO' 'Pz'; 'IO' 'FC1';...
        'IO' 'FC2'; 'IO' 'CP1'; 'IO' 'CP2'; 'IO' 'FC5'; 'IO' 'FC6'; 'IO' 'CP5'; 'IO' 'CP6';...
        'IO' 'FT9'; 'IO' 'FT10'; 'IO' 'TP9';'IO' 'TP10';'IO' 'F1';'IO' 'F2';'IO' 'C1';...
        'IO' 'C2';'IO' 'P1';'IO' 'P2';'IO' 'AF3';'IO' 'AF4';'IO' 'FC3';'IO' 'FC4';...
        'IO' 'CP3';'IO' 'CP4';'IO' 'PO3';'IO' 'PO4';'IO' 'F5';'IO' 'F6';'IO' 'C5';...
        'IO' 'C6';'IO' 'P5';'IO' 'P6';'IO' 'AF7';'IO' 'AF8';'IO' 'FT7';'IO' 'FT8';...
        'IO' 'TP7';'IO' 'TP8';'IO' 'PO7';'IO' 'PO8';'IO' 'Fpz';'IO' 'CPz';'IO' 'POz';...
        'IO' 'Oz';'IO' 'FCz'};
    for i=1:length(subj)
        coh_l_po4_temp{1, i}= ft_connectivityanalysis(cfg, eeg_fft_l_temp{1, i});
        coh_l_po4_temp_val(i,:,:)=coh_l_po4_temp{1, i}.cohspctrm;
    end

    % 左边先开始闪，PO3（45）
    for i=1:length(subj)
        eeg_fft_l_temp{1, i}.fourierspctrm(:,20,:)=eeg_44_l_phase{1, i}.fourierspctrm(:,45,:);
    end

    cfg = [];
    cfg.method  = 'coh';
    cfg.channelcmb = {'IO' 'Fp1'; 'IO' 'Fp2'; 'IO' 'F3'; 'IO' 'F4'; 'IO' 'C3'; 'IO' 'C4';...
        'IO' 'P3'; 'IO' 'P4'; 'IO' 'O1'; 'IO' 'O2'; 'IO' 'F7'; 'IO' 'F8'; 'IO' 'T7';...
        'IO' 'T8'; 'IO' 'P7'; 'IO' 'P8'; 'IO' 'Fz'; 'IO' 'Cz'; 'IO' 'Pz'; 'IO' 'FC1';...
        'IO' 'FC2'; 'IO' 'CP1'; 'IO' 'CP2'; 'IO' 'FC5'; 'IO' 'FC6'; 'IO' 'CP5'; 'IO' 'CP6';...
        'IO' 'FT9'; 'IO' 'FT10'; 'IO' 'TP9';'IO' 'TP10';'IO' 'F1';'IO' 'F2';'IO' 'C1';...
        'IO' 'C2';'IO' 'P1';'IO' 'P2';'IO' 'AF3';'IO' 'AF4';'IO' 'FC3';'IO' 'FC4';...
        'IO' 'CP3';'IO' 'CP4';'IO' 'PO3';'IO' 'PO4';'IO' 'F5';'IO' 'F6';'IO' 'C5';...
        'IO' 'C6';'IO' 'P5';'IO' 'P6';'IO' 'AF7';'IO' 'AF8';'IO' 'FT7';'IO' 'FT8';...
        'IO' 'TP7';'IO' 'TP8';'IO' 'PO7';'IO' 'PO8';'IO' 'Fpz';'IO' 'CPz';'IO' 'POz';...
        'IO' 'Oz';'IO' 'FCz'};
    for i=1:length(subj)
        coh_l_po3_temp{1, i}= ft_connectivityanalysis(cfg, eeg_fft_l_temp{1, i});
        coh_l_po3_temp_val(i,:,:)=coh_l_po3_temp{1, i}.cohspctrm;
    end

    % 右边先开始闪，PO3（45）
    for i=1:length(subj)
        eeg_fft_r_temp{1, i}.fourierspctrm(:,20,:)=eeg_44_r_phase{1, i}.fourierspctrm(:,45,:);
    end

    cfg = [];
    cfg.method  = 'coh';
    cfg.channelcmb = {'IO' 'Fp1'; 'IO' 'Fp2'; 'IO' 'F3'; 'IO' 'F4'; 'IO' 'C3'; 'IO' 'C4';...
        'IO' 'P3'; 'IO' 'P4'; 'IO' 'O1'; 'IO' 'O2'; 'IO' 'F7'; 'IO' 'F8'; 'IO' 'T7';...
        'IO' 'T8'; 'IO' 'P7'; 'IO' 'P8'; 'IO' 'Fz'; 'IO' 'Cz'; 'IO' 'Pz'; 'IO' 'FC1';...
        'IO' 'FC2'; 'IO' 'CP1'; 'IO' 'CP2'; 'IO' 'FC5'; 'IO' 'FC6'; 'IO' 'CP5'; 'IO' 'CP6';...
        'IO' 'FT9'; 'IO' 'FT10'; 'IO' 'TP9';'IO' 'TP10';'IO' 'F1';'IO' 'F2';'IO' 'C1';...
        'IO' 'C2';'IO' 'P1';'IO' 'P2';'IO' 'AF3';'IO' 'AF4';'IO' 'FC3';'IO' 'FC4';...
        'IO' 'CP3';'IO' 'CP4';'IO' 'PO3';'IO' 'PO4';'IO' 'F5';'IO' 'F6';'IO' 'C5';...
        'IO' 'C6';'IO' 'P5';'IO' 'P6';'IO' 'AF7';'IO' 'AF8';'IO' 'FT7';'IO' 'FT8';...
        'IO' 'TP7';'IO' 'TP8';'IO' 'PO7';'IO' 'PO8';'IO' 'Fpz';'IO' 'CPz';'IO' 'POz';...
        'IO' 'Oz';'IO' 'FCz'};
    for i=1:length(subj)
        coh_r_po3_temp{1, i}= ft_connectivityanalysis(cfg, eeg_fft_r_temp{1, i});
        coh_r_po3_temp_val(i,:,:)=coh_r_po3_temp{1, i}.cohspctrm;
    end

    % 右边先开始闪，PO4（46）
    for i=1:length(subj)
        eeg_fft_r_temp{1, i}.fourierspctrm(:,20,:)=eeg_42_r_phase{1, i}.fourierspctrm(:,46,:);
    end

    cfg = [];
    cfg.method  = 'coh';
    cfg.channelcmb = {'IO' 'Fp1'; 'IO' 'Fp2'; 'IO' 'F3'; 'IO' 'F4'; 'IO' 'C3'; 'IO' 'C4';...
        'IO' 'P3'; 'IO' 'P4'; 'IO' 'O1'; 'IO' 'O2'; 'IO' 'F7'; 'IO' 'F8'; 'IO' 'T7';...
        'IO' 'T8'; 'IO' 'P7'; 'IO' 'P8'; 'IO' 'Fz'; 'IO' 'Cz'; 'IO' 'Pz'; 'IO' 'FC1';...
        'IO' 'FC2'; 'IO' 'CP1'; 'IO' 'CP2'; 'IO' 'FC5'; 'IO' 'FC6'; 'IO' 'CP5'; 'IO' 'CP6';...
        'IO' 'FT9'; 'IO' 'FT10'; 'IO' 'TP9';'IO' 'TP10';'IO' 'F1';'IO' 'F2';'IO' 'C1';...
        'IO' 'C2';'IO' 'P1';'IO' 'P2';'IO' 'AF3';'IO' 'AF4';'IO' 'FC3';'IO' 'FC4';...
        'IO' 'CP3';'IO' 'CP4';'IO' 'PO3';'IO' 'PO4';'IO' 'F5';'IO' 'F6';'IO' 'C5';...
        'IO' 'C6';'IO' 'P5';'IO' 'P6';'IO' 'AF7';'IO' 'AF8';'IO' 'FT7';'IO' 'FT8';...
        'IO' 'TP7';'IO' 'TP8';'IO' 'PO7';'IO' 'PO8';'IO' 'Fpz';'IO' 'CPz';'IO' 'POz';...
        'IO' 'Oz';'IO' 'FCz'};
    for i=1:length(subj)
        coh_r_po4_temp{1, i}= ft_connectivityanalysis(cfg, eeg_fft_r_temp{1, i});
        coh_r_po4_temp_val(i,:,:)=coh_r_po4_temp{1, i}.cohspctrm;
    end

    coh_l_po4_shuf(n,:,:)=mean(coh_l_po4_temp_val,1);
    coh_l_po3_shuf(n,:,:)=mean(coh_l_po3_temp_val,1);
    coh_r_po3_shuf(n,:,:)=mean(coh_r_po3_temp_val,1);
    coh_r_po4_shuf(n,:,:)=mean(coh_r_po4_temp_val,1);
    n
end
delete(parob);
save('coh_shuf_single.mat','coh_l_po4_shuf','coh_l_po3_shuf','coh_r_po3_shuf','coh_r_po4_shuf');

% 取百分位数
load('coh_shuf_single.mat');
coh_l_po3_shuf_temp=squeeze(coh_l_po3_shuf(:,:,22));
coh_l_po4_shuf_temp=squeeze(coh_l_po4_shuf(:,:,13));
coh_r_po3_shuf_temp=squeeze(coh_r_po3_shuf(:,:,22));
coh_r_po4_shuf_temp=squeeze(coh_r_po4_shuf(:,:,13));

coh_l_po3_val_m=mean(coh_l_po3_val,2);
coh_l_po4_val_m=mean(coh_l_po4_val,2);
coh_r_po3_val_m=mean(coh_r_po3_val,2);
coh_r_po4_val_m=mean(coh_r_po4_val,2);

for i=1:64
    coh_l_po3_p(i) = invprctile(coh_l_po3_shuf_temp(:,i),coh_l_po3_val_m(i));
    coh_l_po4_p(i) = invprctile(coh_l_po4_shuf_temp(:,i),coh_l_po4_val_m(i));
    coh_r_po3_p(i) = invprctile(coh_r_po3_shuf_temp(:,i),coh_r_po3_val_m(i));
    coh_r_po4_p(i) = invprctile(coh_r_po4_shuf_temp(:,i),coh_r_po4_val_m(i));    
end

coh_l_p=(coh_l_po3_p+coh_l_po4_p)/2;
coh_r_p=(coh_r_po3_p+coh_r_po4_p)/2;

% topo plot
chann_label=eegdata{1, 1}.label;
chann_label(20)=[];

cfg = [];
cfg.layout = 'easycapM1.mat';
cfg.xlim = 'maxmin';  % time limitation
cfg.zlim = [50 100];
cfg.style  = 'straight';
cfg.colormap     = brewermap([],'Reds');
cfg.colorbar = 'yes';
% cfg.marker = 'on';
% cfg.comment = 'xlim';
% cfg.commentpos = 'lefttop';

coh_l_plot = {}; 
coh_l_plot.avg = coh_l_p'; % data definition, channel * time
coh_l_plot.time = 1; % time definition
coh_l_plot.dimord = 'chan_time';
coh_l_plot.label = chann_label;

cfg.highlight          = 'on';
cfg.highlightchannel   =  find(coh_l_p'>95);
cfg.highlightsymbol    = '*';
cfg.highlightcolor     = [1 1 0];
% cfg.highlightsize      = 6;

ft_topoplotER(cfg, coh_l_plot);
title('coh_l');
axis tight;

cfg = [];
cfg.layout = 'easycapM1.mat';
cfg.xlim = 'maxmin';  % time limitation
cfg.zlim = [50 100];
cfg.style  = 'straight';
cfg.colormap     = brewermap([],'Reds');
cfg.colorbar = 'yes';
% cfg.marker = 'on';
% cfg.comment = 'xlim';
% cfg.commentpos = 'lefttop';

coh_r_plot = {}; 
coh_r_plot.avg = coh_r_p'; % data definition, channel * time
coh_r_plot.time = 1; % time definition
coh_r_plot.dimord = 'chan_time';
coh_r_plot.label = chann_label;

cfg.highlight          = 'on';
cfg.highlightchannel   =  find(coh_r_p'>95);
cfg.highlightsymbol    = '*';
cfg.highlightcolor     = [1 1 0];

ft_topoplotER(cfg, coh_r_plot);
title('coh_r');
axis tight;

save('coh_single.mat','coh_l_p','coh_r_p');

%% 5. behaviral correlates，两个SSVEP的相位差异与行为之间的相关。circular-linear相关
phase_diff=circ_mean([phase_l_diff; phase_r_diff]);

load('D:\OneDrive - hznu.edu.cn\Matlab_workspace\Attention and oscillation\kongqing\Data\beha\exp1\single_acc.mat');
load('D:\OneDrive - hznu.edu.cn\Matlab_workspace\Attention and oscillation\kongqing\Data\beha\exp1\single_rt.mat');

% acc
for i=1:length(subj)
    acc_sub(i)=mean(single_acc{1, subj(i)});  
end

% RT
for i=1:length(subj)
    for j=1:320
        if single_acc{1, subj(i)}(1,j) == 0
            single_rt{1, subj(i)}(1,j) = NaN;
        end
    end
    rt_sub(i)=mean(single_rt{1, subj(i)},"omitnan");  
end

% correlation
[rho pval] = circ_corrcl(phase_diff', rt_sub');
[rho pval] = circ_corrcl(abs(phase_diff'), rt_sub');

% plot
scatter(phase_diff,rt_sub);
scatter(abs(phase_diff),rt_sub);

save('beha_single.mat','phase_diff','rt_sub');
