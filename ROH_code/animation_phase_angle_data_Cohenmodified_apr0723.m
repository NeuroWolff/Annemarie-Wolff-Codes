%% Animation of phase angle data, timeseries and polar plots, modified from uANTS_synch.m from Cohen 0 to hero

%% exact code of Cohen
clear; clc;
% load in V1 mouse data
load(fullfile('C:/Users/awolf/Documents/programming_courses/Cohen_courses/Neural_signal_processing_and_analysis_0_to_hero/9_synchronization/','v1_laminar.mat'));

% setup connectivity parameters
% The goal here is to explore the mechanism of phase synchronization.
% extract phase angle time series from two channels from trial 1 at 8 Hz.

% channels for connectivity
chan1idx = 1;
chan2idx = 8;

% create a complex Morlet wavelet (don't forget to plot the wavelet!!!)
wave_freq = 8;
time      = -1.5:1/srate:1.5;
s         = 8/(2*pi*wave_freq);
wavelet   = exp(2*1i*pi*wave_freq.*time) .* exp(-time.^2./(2*s^2));
half_wavN = (length(time)-1)/2;

% FFT parameters
nWave = length(time); % time vector for wavelet
nData = size(csd, 2); % time vector for data, only 1 trial
nConv = nWave + nData - 1;

% FFT of wavelet (check nfft)
waveletX = fft(wavelet, nConv);
waveletX = waveletX ./ max(waveletX);

% initialize output time-frequency data
phase_data = zeros(2, length(timevec));
real_data  = zeros(2, length(timevec));

% analytic signal of channel 1
dataX = fft(csd(chan1idx, :, 1), nConv); % take FFT
as = ifft(waveletX .* dataX, nConv); % take inverse FFT
as = as(half_wavN+1:end-half_wavN); % cut off wings of convolution

% collect real and phase data
phase_data(1,:) = angle(as); % extract phase angles
real_data(1,:)  = real(as);  % extract the real part (projection onto real axis)

% analytic signal of channel 2
dataX = fft(csd(chan2idx, :, 1), nConv);
as = ifft(waveletX.*dataX,nConv);
as = as(half_wavN+1:end-half_wavN);

% collect real and phase data
phase_data(2,:) = angle(as);
real_data(2,:)  = real(as);

% setup figure and define plot handles
% note: This cell is just setting up the figure for the following cell. 
%       You can run it and move on.

% open and name figure
figure(1), clf
set(gcf,'NumberTitle','off','Name','Movie magic minimizes the magic.');

% draw the filtered signals
subplot(221)
filterplotH1 = plot(timevec(1),real_data(1,1),'b');
hold on
filterplotH2 = plot(timevec(1),real_data(2,1),'m');
set(gca,'xlim',[timevec(1) timevec(end)],'ylim',[min(real_data(:)) max(real_data(:))])
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
title([ 'Filtered signal at ' num2str(wave_freq) ' Hz' ])

% draw the phase angle time series
subplot(222)
phaseanglesH1 = plot(timevec(1),phase_data(1,1),'b');
hold on
phaseanglesH2 = plot(timevec(1),phase_data(2,1),'m');
set(gca,'xlim',[timevec(1) timevec(end)],'ylim',[-pi pi]*1.1)
xlabel('Time (ms)')
ylabel('Phase angle (radian)')
title('Phase angle time series')

% draw phase angles in polar space
subplot(223)
polar2chanH1 = polar([zeros(1,1) phase_data(1,1)]',repmat([0 1],1,1)','b');
hold on
polar2chanH2 = polar([zeros(1,1) phase_data(2,1)]',repmat([0 1],1,1)','m');
title('Phase angles from two channels')

% draw phase angle differences in polar space
subplot(224)
polarAngleDiffH = polar([zeros(1,1) phase_data(2,1)-phase_data(1,1)]',repmat([0 1],1,1)','k');
title('Phase angle differences from two channels')

% ANIMATE: now update plots at each timestep
for pointi=1:5:length(timevec)
    
    % update filtered signals
    set(filterplotH1,'XData',timevec(1:pointi),'YData',real_data(1,1:pointi))
    set(filterplotH2,'XData',timevec(1:pointi),'YData',real_data(2,1:pointi))
    
    % update cartesian plot of phase angles
    set(phaseanglesH1,'XData',timevec(1:pointi),'YData',phase_data(1,1:pointi))
    set(phaseanglesH2,'XData',timevec(1:pointi),'YData',phase_data(2,1:pointi))
    
    subplot(223), cla
    polar([zeros(1,pointi) phase_data(1,1:pointi)]',repmat([0 1],1,pointi)','b');
    hold on
    polar([zeros(1,pointi) phase_data(2,1:pointi)]',repmat([0 1],1,pointi)','m');
    
    subplot(224), cla
    polar([zeros(1,pointi) phase_data(2,1:pointi)-phase_data(1,1:pointi)]',repmat([0 1],1,pointi)','k');
    
    drawnow
end

%% edit code for my data
clear; clc;
% load data
data_max_file = load(fullfile('E:/EEGdata/Knott_data/Knott_P300_data/','oddball_CON_16_filtpt1_80_noflatnoisychans_cleanline_avgREF_ICAMARA.mat'));
data_min_file = load(fullfile('E:/EEGdata/Knott_data/Knott_P300_data/','oddball_SCZ_01_filtpt1_80_noflatnoisychans_cleanline_avgREF_ICAMARA.mat'));

% specify electrodes and stims
eles = {'Fp1','Fpz','Fp2'}; stims = {'S  1','S  2','S  3'};
idx_eles = ismember({data_max_file.EEG.chanlocs.labels}',eles);
srate = data_max_file.EEG.srate;

% Knott oddball precision index from 3-5Hz and 250-350ms
wave_freq = 4; % Hz
polar_time = 0.3 * srate; % datapoints for time offset

% specify number of datapoints to keep
num_mins = 0.5;
num_pts = num_mins * 60 * data_max_file.EEG.srate;

% specify stim onsets for max file
lats = cell2mat({data_max_file.EEG.event.latency}');
stims_max = lats(ismember({data_max_file.EEG.event.type}',stims));
idx_stim_onset_max = stims_max(1);
stims_max = stims_max - stims_max(1) + 1; 
stims_max(stims_max > num_pts - polar_time) = [];
clear lats
% specify stim onsets for min file
lats = cell2mat({data_min_file.EEG.event.latency}');
stims_min = lats(ismember({data_min_file.EEG.event.type}',stims));
idx_stim_onset_min = stims_min(1); 
stims_min = stims_min - stims_min(1) + 1;
stims_min(stims_min > num_pts - polar_time) = [];
clear lats

% isolate data for wavelet analysis
% data_max = nanmean(data_max_file.EEG.data(idx_eles,idx_stim_onset_max:idx_stim_onset_max + num_pts),1);
% data_min = nanmean(data_min_file.EEG.data(idx_eles,idx_stim_onset_min:idx_stim_onset_min + num_pts),1);
data_max = nanmean(data_max_file.EEG.data(idx_eles,:),1);
data_min = nanmean(data_min_file.EEG.data(idx_eles,:),1);

%% analysis for max data
% create a complex Morlet wavelet (don't forget to plot the wavelet!!!)
time      = -1.5:1/srate:1.5;
s         = 8/(2*pi*wave_freq);
wavelet   = exp(2*1i*pi*wave_freq .* time) .* exp(-time .^ 2 ./ (2*s^2));
half_wavN = (length(time)-1)/2;

% FFT parameters
nWave = length(time); % time vector for wavelet
nData = size(data_max, 2); % time vector for data, only 1 trial
nConv = nWave + nData - 1;

% FFT of wavelet (check nfft)
waveletX = fft(wavelet, nConv);
waveletX = waveletX ./ max(waveletX);

% analytic signal of data_max
dataX = fft(data_max, nConv); % take FFT
as = ifft(waveletX .* dataX, nConv); % take inverse FFT
as = as(half_wavN+1:end-half_wavN); % cut off wings of convolution

% collect real and phase data
data_phase_max = angle(as(:,idx_stim_onset_max:idx_stim_onset_max + num_pts)); % extract phase angles
data_real_max  = real(as(:,idx_stim_onset_max:idx_stim_onset_max + num_pts));  % extract the real part (projection onto real axis)
clear as

%% analysis for min data
% FFT parameters
nWave = length(time); % time vector for wavelet
nData = size(data_min, 2); % time vector for data, only 1 trial
nConv = nWave + nData - 1;

% FFT of wavelet (check nfft)
waveletX = fft(wavelet, nConv);
waveletX = waveletX ./ max(waveletX);

% analytic signal of data_min
dataX = fft(data_min, nConv);
as = ifft(waveletX.*dataX,nConv);
as = as(half_wavN+1:end-half_wavN);

% collect real and phase data
data_phase_min = angle(as(:,idx_stim_onset_min:idx_stim_onset_min + num_pts)); % extract phase angles
data_real_min  = real(as(:,idx_stim_onset_min:idx_stim_onset_min + num_pts));  % extract the real part (projection onto real axis)
timevec = (data_max_file.EEG.times(idx_stim_onset_max:idx_stim_onset_max + num_pts) - data_max_file.EEG.times(idx_stim_onset_max)) ./ 1000;
clear as dataX waveletX nWave nData nConv as half_wavN s

%% plot data
interval = 1.0 * srate;
col_max = [0, 0.4470, 0.7410]; col_min = [0.3010, 0.7450, 0.9330];
lin_wid = 1.5; font_tit = 13; font_ax = font_tit-2;
% setup figure and define plot handles
clear data_max_file data_min
rows = 3; cols = 2;
figure('Renderer', 'painters', 'Position', [10 10 700 700]), clf
set(gcf,'NumberTitle','off','Name','Movie magic minimizes the magic.');

% draw the filtered signals
ax1=subplot(rows,cols,1);
filterplotH1 = plot(timevec(1),data_real_max(1,1),'Color',col_max,'LineWidth',lin_wid);
set(gca,'ylim',[min(data_real_max(:))-1 max(data_real_max(:))+1],'XTick',[],'YTick',[]);
box off;
xlabel('Time','FontSize',font_ax); ylabel('Voltage','FontSize',font_ax);
title({'LOW VARIABILITY OVER TRIALS';'Signal amplitude at 4Hz'},'FontSize',font_tit,'FontWeight','normal')

% draw the filtered signals
ax2=subplot(rows,cols,2);
filterplotH2 = plot(timevec(1),data_real_min(1,1),'Color',col_min,'LineWidth',lin_wid);
set(gca,'ylim',[min(data_real_min(:))-1 max(data_real_min(:))+1],'XTick',[],'YTick',[])
box off;
xlabel('Time','FontSize',font_ax); ylabel('Voltage','FontSize',font_ax);
title({'HIGH VARIABILITY OVER TRIALS';'Signal amplitude at 4Hz'},'FontSize',font_tit,'FontWeight','normal')

% draw the phase angle time series
ax3=subplot(rows,cols,cols+1);
phaseanglesH1 = plot(timevec(1),data_phase_max(1,1),'Color',col_max,'LineWidth',lin_wid);
set(gca,'ylim',[-pi-1 pi+1]*1.1,'XTick',[],'YTick',[]);
box off;
xlabel('Time','FontSize',font_ax); ylabel('Phase angle','FontSize',font_ax);
title('Phase angle time series at 4Hz','FontSize',font_tit,'FontWeight','normal')

% draw the phase angle time series
ax4=subplot(rows,cols,cols+2);
phaseanglesH2 = plot(timevec(1),data_phase_min(1,1),'Color',col_min,'LineWidth',lin_wid);
set(gca,'ylim',[-pi-1 pi+1]*1.1,'XTick',[],'YTick',[]);
box off;
xlabel('Time','FontSize',font_ax); ylabel('Phase angle','FontSize',font_ax);
title('Phase angle time series at 4Hz','FontSize',font_tit,'FontWeight','normal')

% draw phase angles in polar space
polar_plot = 1;
subplot(rows,cols,2*cols+1);
% hold all;
polar2chanH1 = polarplot([zeros(1,polar_plot) simdata_polar(1,1:polar_plot)]',repmat([0 1],1,polar_plot)','Color',col_max,'LineWidth',lin_wid);
% polar2chanH1 = polarplot([zeros(1,1) data_phase_max(1,1)]',repmat([0 1],1,1)','w');
set(gca,'Rlim',[0 1],'RTick',[0 1],'RTickLabel',[],'ThetaTick', [0 90 180 270],'ThetaTickLabels',{'0','\pi/2','\pi','3\pi/2'},'FontSize',font_tit);
title('Phase angles in individual trials','FontSize',font_tit,'FontWeight','normal');

% draw phase angles in polar space
subplot(rows,cols,2*cols+2);
% hold on;
polar2chanH2 = polarplot([zeros(1,polar_plot) simdata_polar(1,1:polar_plot)]',repmat([0 1],1,polar_plot)','Color',col_min,'LineWidth',lin_wid);
% polar2chanH2 = polarplot([zeros(1,1) data_phase_min(1,1)]',repmat([0 1],1,1)','w');
set(gca,'Rlim',[0 1],'RTick',[0 1],'RTickLabel',[],'ThetaTick', [0 90 180 270],'ThetaTickLabels',{'0','\pi/2','\pi','3\pi/2'},'FontSize',font_tit);
title('Phase angles in individual trials','FontSize',font_tit,'FontWeight','normal');

% simulate data for polar plots
num_figs = 2;
circle_proportions = linspace(0.54,.92,num_figs);
num_trials = 50;
% simulate data
for simi = 1:num_figs
    % generate phase-angle distribution
    simdata_polar(simi,:) = rand(1,num_trials) * (2*pi) * circle_proportions(simi);
    % compute ITPC and preferred phase-angle
    itpc_polar(simi,:) = abs(mean(exp(1i*simdata_polar(simi,:))));
    prefAngle_polar(simi,:) = angle(mean(exp(1i*simdata_polar(simi,:))));
end

% setup movie
cd('C:/Users/awolf/Downloads/');
vid_filename = 'ITPC_amp_phaseangles_animation_apr1223.avi';
mov = VideoWriter(vid_filename,'Motion JPEG AVI'); % create video writer object
open(mov);

% ANIMATE: now update plots at each timestep
step = 10;
total_pnts = numel(1:step:length(timevec)-interval);
for pointi = 1:step:length(timevec)-interval
    % update filtered signals
    if pointi >= interval
        set(filterplotH1,'XData',timevec(pointi-interval:pointi+interval),'YData',data_real_max(1,pointi-interval:pointi+interval)); set(ax1,'xlim',[timevec(pointi-interval) timevec(pointi+interval)]);
        set(filterplotH2,'XData',timevec(pointi-interval:pointi+interval),'YData',data_real_min(1,pointi-interval:pointi+interval)); set(ax2,'xlim',[timevec(pointi-interval) timevec(pointi+interval)]);
    else
        set(filterplotH1,'XData',timevec(1:pointi+interval),'YData',data_real_max(1,1:pointi+interval)); set(ax1,'xlim',[timevec(1) timevec(pointi+interval)]);
        set(filterplotH2,'XData',timevec(1:pointi+interval),'YData',data_real_min(1,1:pointi+interval)); set(ax2,'xlim',[timevec(1) timevec(pointi+interval)]);
    end

    % update cartesian plot of phase angles
    if pointi >= interval
        set(phaseanglesH1,'XData',timevec(pointi-interval:pointi+interval),'YData',data_phase_max(1,pointi-interval:pointi+interval)); set(ax3,'xlim',[timevec(pointi-interval) timevec(pointi+interval)]);
        set(phaseanglesH2,'XData',timevec(pointi-interval:pointi+interval),'YData',data_phase_min(1,pointi-interval:pointi+interval)); set(ax4,'xlim',[timevec(pointi-interval) timevec(pointi+interval)]);
    else
        set(phaseanglesH1,'XData',timevec(1:pointi+interval),'YData',data_phase_max(1,1:pointi+interval)); set(ax3,'xlim',[timevec(1) timevec(pointi+interval)]);
        set(phaseanglesH2,'XData',timevec(1:pointi+interval),'YData',data_phase_min(1,1:pointi+interval)); set(ax4,'xlim',[timevec(1) timevec(pointi+interval)]);
    end

    % plot polar plots of phase data
    if mod(pointi,300) == 1
        % MAX
        subplot(rows,cols,2*cols+1);
        hold all;
        polarplot([zeros(1,polar_plot) simdata_polar(1,1:polar_plot)]',repmat([0 1],1,polar_plot)','Color',col_max,'LineWidth',lin_wid);
        % MIN
        subplot(rows,cols,2*cols+2);
        hold all;
        polarplot([zeros(1,polar_plot) simdata_polar(2,1:polar_plot)]',repmat([0 1],1,polar_plot)','Color',col_min,'LineWidth',lin_wid);
        % COUNT
        polar_plot = polar_plot + 1;
    end

    % update movie frame
    writeVideo(mov,getframe(gcf));

    % animate
    pause(.005); % drawnow;

end

% MAX
subplot(rows,cols,2*cols+1);
hold on;
polarplot([0 prefAngle_polar(1,:)],[0 itpc_polar(1,:)],'r','LineWidth',lin_wid*3);
annotation('textarrow',[.2 .295],[.15 .195],'String',['ITPC = ' num2str(round(itpc_polar(1,:),2))],'FontSize',font_tit,'FontWeight','bold','LineWidth',lin_wid*1.5);
% MIN
subplot(rows,cols,2*cols+2);
hold on;
polarplot([0 prefAngle_polar(2,:)],[0 itpc_polar(2,:)],'r','LineWidth',lin_wid*3);
annotation('textarrow',[.86 .735],[.15 .195],'String',['ITPC = ' num2str(round(itpc_polar(2,:),2))],'FontSize',font_tit,'FontWeight','bold','LineWidth',lin_wid*1.5);

% update movie frame
writeVideo(mov,getframe(gcf));
% write out movie
close(mov);

