% Load data structures
clear; clc;
load(fullfile...
    ('E:\EEGdata\Knott_data\Knott_P300_data\',...
    'knott_oddball_ITPC_eleclusters_dec3122.mat'));
disp('Data structures loaded.');

%% Create Simulated data
% parameters
clear EEG
num_parts = 50;
EEG.srate  = 500; % sampling rate in Hz
EEG.times   = -1:1/EEG.srate:1.2; % -1s to 1.2s
EEG.pnts   = length(EEG.times);
EEG.trials = 50;

% continuous sinewave frequencies and corresponding amplitudes
freq_range = [1 80];
frex = logspace(freq_range(1), freq_range(2), 2*diff(freq_range));
amps = logspace(10, 1, length(frex));

% phase offsets
num_intervals = 100;
phase_offsets = linspace(0, 2, num_intervals);
idx.before_onset = find(EEG.times==0)-1; % determine which datapoint is just before stim onset

%%%%% WHITE NOISE
noise_white_gain = 1.0;

%%%%% PINK NOISE
exponen_decay = 50; % exponential decay parameter
noiseamp = 1.0; % noise gain

% create Gaussian taper for transient oscillation
peaktime = .175; % seconds
width    = .15; % FWHM in seconds
sinefreq = linspace(4,8,5); % for sine wave
gaus = exp(-(EEG.times - peaktime) .^2 / (2*width^2)); % Gaussian in the time domain
clear EEG.data

%%%%%%%%% MAKE sine wave
% Loop over phase offsets
for offseti = 1:length(phase_offsets)
    % Lopp over participants
    for parti = 1:max(num_parts) 
        % Loop over trials
        for triali = 1:EEG.trials

            %%% create transient oscillation
            sinewave = zeros(1, EEG.pnts);
            % Loop over Gaussian kernel peak frequencies
            for si=1:length(sinefreq)
                sinewave = sinewave + sin(2*pi*sinefreq(si)*EEG.times + phase_offsets(offseti)*pi*rand);
            end

            %%% Make Pink Noise
            amplitud_spect = rand(1,floor(EEG.pnts/2)-1) .* exp(-(1:floor(EEG.pnts/2)-1)/exponen_decay);
            amplitud_spect = [amplitud_spect(1) amplitud_spect 0 amplitud_spect(:,end:-1:1)];
            % Fourier coefficients
            fourier_coeffs = amplitud_spect .* exp(1i*2*pi*rand(size(amplitud_spect)));
            % inverse Fourier transform to create the noise
            mat_noise_pink = real(ifft(fourier_coeffs)) * EEG.pnts;
            % make pink noise
            noise_pink = noiseamp*[mat_noise_pink mat_noise_pink(end)];

            %%% Make White Noise
            noise_white = noise_white_gain*randn(1, EEG.pnts);

            % Combine all into one signal
            EEG.data(offseti, parti, :, triali) = (sinewave .* gaus) + noise_pink + noise_white;
        end
    end
    disp(['Offset ' num2str(offseti) ' is done.']);
end

% save in data structure
simData.data_50t = EEG.data;
simData.params = {'phase_offsets', 'participants', 'timepoints', 'trials'};
simData.srate  = 500; % sampling rate in Hz
simData.times   = (-1:1/simData.srate:1.2)*1000; % -1s to 1.2s
simData.pnts   = length(simData.times);
disp('Simulated Data generation finished.');

% Run timefreq analysis
epoch_interval = [EEG.times(1)*1000 EEG.times(end)*1000];
base_interv = [-400 0];
% freq_range = [1 80];
nFrex = diff(freq_range)*2;
range_cycles = [3 10];
% initialize matrices
simData.ITPC_50t = zeros(length(phase_offsets), max(num_parts), nFrex, EEG.pnts);
% TF analysis
% for offseti = 1:length(phase_offsets)
%     for parti = 1:max(num_parts)
%         tmp_data = reshape(EEG.data(offseti, parti, :, :), [1, EEG.pnts, EEG.trials]);
%         [~,~,~,simData.ITPC_50t(offseti,parti,:,:),~,~,tmp_frex] = cohen_MorletWaveTF(tmp_data(1, :, :),...
%                                                                                     EEG.srate,...
%                                                                                     epoch_interval,...
%                                                                                     base_interv,...
%                                                                                     freq_range,...
%                                                                                     nFrex,...
%                                                                                     range_cycles, ...
%                                                                                     1);
%         simData.frex = tmp_frex;
%         clear tmp_data
%     end
%     disp(['Offset ' num2str(offseti) ' finished.']);
%     clear elapsed_time
% end
% disp('Time Freq analysis done.');

% extract ITPC data of interest from simulations and real data
% load('ITPCsim_100offsets_50parts_jan0323.mat');
phase_offset_labs = {'None', '\pi/2', '\pi', '3\pi/2', '2\pi'};
phase_offset_vals = linspace(0,2*pi,num_intervals);
hz2get = [4 8];
time2get = [0 122];
idx_hz2get_sim = dsearchn(simData.frex',hz2get');
idx_time2get_sim = dsearchn(simData.times',time2get');
hz2get = [3 5];
time2get = [200 300];
idx_hz2get_knott = dsearchn(mat_ITPC.frex',hz2get');
idx_time2get_knott = dsearchn(mat_ITPC.times',time2get');

% extract data
extracted_simITPC = mean(...
                    mean(...
                    simData.ITPC_50t(:,:,idx_hz2get_sim(1):idx_hz2get_sim(2),idx_time2get_sim(1):idx_time2get_sim(2))...
                    ,3)...
                    ,4);

extracted_KnottITPC = mean(...
                      mean(...
                      mat_ITPC.devi(:,1,idx_hz2get_knott(1):idx_hz2get_knott(2),idx_time2get_knott(1):idx_time2get_knott(2))...
                      ,3)...
                      ,4);

%% plot results
x_ticks = [0 pi/2 pi (3*pi)/2 2*pi];
x_ticks_labs = {'0','\pi/2','\pi','3\pi/2','2\pi'};
figure; set(gcf,'Position',[5 5 900 700]);
hold on;
scatter(phase_offset_vals,extracted_simITPC,'xk');
scatter(phase_offset_vals,mean(extracted_simITPC,2),'^r', 'filled');
set(gca,'xlim',[phase_offset_vals(1)-.25 phase_offset_vals(end)+.25],'XTick',x_ticks,'XTickLabel',x_ticks_labs);
xlabel('Phase offset','FontSize',12);
ylabel('ITPC','FontSize',12);
% save(['ITPCsim_' num2str(num_intervals) 'offsets_' num2str(num_parts) 'parts_jan0323.mat'],'simData');

%% fit model to simulated data in order to predict phase offsets for real data
clear data_for_modelfit_x data_for_modelfit_y
data_for_modelfit_x = []; data_for_modelfit_y = [];
for intervi = 1:size(extracted_simITPC,1)
    data_for_modelfit_x = [data_for_modelfit_x; intervi*ones([size(extracted_simITPC,2) 1])];
    data_for_modelfit_y = [data_for_modelfit_y; squeeze(extracted_simITPC(intervi,:))'];
end

% specify parameters of fit
% for getting the betas and the delta
fit_poly1 = 1; fit_poly2 = 3;
% for getting the goodness of fit
fit_gof_1 = 'poly1'; fit_gof_2 = 'poly3';
lin_wid = 2;

% fit data to curve
% to get the betas and the delta
% fit 1
[fit1_polynomials,Structure_fit1] = polyfit(data_for_modelfit_x,data_for_modelfit_y,fit_poly1);
[yvals_fit1,delta_fit1] = polyval(fit1_polynomials,data_for_modelfit_x,Structure_fit1);
% [fit1_polynomials,Structure_fit1] = polyfit(data_for_modelfit_y,data_for_modelfit_x,fit_poly1);
% [yvals_fit1,delta_fit1] = polyval(fit1_polynomials,data_for_modelfit_y,Structure_fit1);
% fit 2
[fit2_polynomials,Structure_fit2] = polyfit(data_for_modelfit_x,data_for_modelfit_y,fit_poly2);
[yvals_fit2,delta_fit2] = polyval(fit2_polynomials,data_for_modelfit_x,Structure_fit2);
% [fit2_polynomials,Structure_fit2] = polyfit(data_for_modelfit_y,data_for_modelfit_x,fit_poly2);
% [yvals_fit2,delta_fit2] = polyval(fit2_polynomials,data_for_modelfit_y,Structure_fit2);
clear Structure_*

% to get the goodness of fit
% fit 1
[fitobject_1,gof_1] = fit(data_for_modelfit_x,data_for_modelfit_y,fit_gof_1);
% [fitobject_1,gof_1] = fit(data_for_modelfit_y,data_for_modelfit_x,fit_gof_1);
% fit 2
[fitobject_2,gof_2] = fit(data_for_modelfit_x,data_for_modelfit_y,fit_gof_2);
% [fitobject_2,gof_2] = fit(data_for_modelfit_y,data_for_modelfit_x,fit_gof_2);

% plot the results
hold on;
scatter(data_for_modelfit_x,data_for_modelfit_y,'kx');
plot(data_for_modelfit_x,yvals_fit1,'r-','LineWidth',lin_wid);
plot(data_for_modelfit_x,yvals_fit1+2*delta_fit1,'r:','LineWidth',lin_wid);
plot(data_for_modelfit_x,yvals_fit1-2*delta_fit1,'r:','LineWidth',lin_wid);
plot(data_for_modelfit_x,yvals_fit2,'r-','LineWidth',lin_wid);
plot(data_for_modelfit_x,yvals_fit2+2*delta_fit2,'r:','LineWidth',lin_wid);
plot(data_for_modelfit_x,yvals_fit2-2*delta_fit2,'r:','LineWidth',lin_wid);

%% simulate data for polar plots
num_figs = 9;
circle_proportions = linspace(0.05,1,num_figs);
num_trials = 50;

% simulate data
for simi = 1:num_figs
    % generate phase-angle distribution
    simdata_polar(simi,:) = rand(1,num_trials) * (2*pi) * circle_proportions(simi);
    % compute ITPC and preferred phase-angle
    itpc_polar(simi,:) = abs(mean(exp(1i*simdata_polar(simi,:))));
    prefAngle_polar(simi,:) = angle(mean(exp(1i*simdata_polar(simi,:))));
    disp(['Simulation ' num2str(simi) ' is done.']);
end

polar_tits = {'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'};
font_tit = 14;
rows = 7; cols = 9;
lin_wid = 3;
itpc_intervals = round(linspace(1,100,num_figs));
sims_itpc_data = squeeze(mean(simData.ITPC_50t(itpc_intervals,:,:,:),2));
y_lim = [2 30]; x_lim = [-250 400]; c_lim = [.1 .82];
x_ticks = [0 pi/4 pi/2 (3*pi)/4 pi (5*pi)/4 (3*pi)/2 (7*pi)/4 2*pi];

% plot
% figure; set(gcf,'Position',[5 5 900 700]);
figure; set(gcf,'Position',[5 5 700 700]);
for ploti = 1:num_figs
    subplot(rows,cols,ploti);
    polarplot([zeros(1,num_trials); simdata_polar(ploti,:)],[zeros(1,num_trials); ones(1,num_trials)],'k');
    hold on;
    polarplot([0 prefAngle_polar(ploti,:)],[0 itpc_polar(ploti,:)],'r','LineWidth',lin_wid);
    set(gca,'rlim',[0 1],'RTick',[0 1],'RTickLabel',[],'ThetaTick', [0 90 180 270],'ThetaTickLabels',[]);
    title(polar_tits{ploti},'FontSize',font_tit,'FontWeight','normal');
end

for ploti = 1:num_figs
    subplot(rows,cols,cols+ploti);
    contourf(simData.times,simData.frex,squeeze(sims_itpc_data(ploti,:,:)),40,'LineStyle','none');
    set(gca,'xlim',x_lim,'ylim',y_lim,'YScale','log','clim',c_lim,'YTick',[2 4 8 13 30],'colormap',parula,'XTickLabels',[],'FontSize',font_tit-3);
    vline(0,':w');
    if ploti ~= 1
        set(gca,'YTickLabel',[]);
    end
    if ploti == 1
        ylabel('Freq (Hz)','FontSize',font_tit-2);
    end
end

subplot(rows,cols,2*cols+1:cols*rows);
hold on;
scatter(phase_offset_vals,extracted_simITPC,'xk');
set(gca,'xlim',[phase_offset_vals(1)-.25 phase_offset_vals(end)+.25],'XTick',x_ticks,'XTickLabel',polar_tits,'FontSize',font_tit-2);
xlabel('Simulation phase offset','FontSize',font_tit-2);
ylabel('ITPC','FontSize',font_tit-2);

%% plot phase offset and ITPC results from sim, with two models
rows=2;cols=1;
y_lim = [0 .8];
color_model = [0.4660, 0.6740, 0.1880];
color_dots = [0.5, 0.5, 0.5];
size_dots = 10;
font_tit = 15;
[equ_fit1,~] = polyfit(data_for_modelfit_y,data_for_modelfit_x,fit_poly1);
[~,equ_gof1] = fit(data_for_modelfit_y,data_for_modelfit_x,fit_gof_1);
[equ_fit2,~] = polyfit(data_for_modelfit_y,data_for_modelfit_x,fit_poly2);
[~,equ_gof2] = fit(data_for_modelfit_y,data_for_modelfit_x,fit_gof_2);

% plot
figure; set(gcf,'Position',[5 5 550 800]);
subplot(rows,cols,1);
hold on;
scatter(data_for_modelfit_x,data_for_modelfit_y,size_dots,color_dots,'.');
plot(data_for_modelfit_x,yvals_fit1,'b-','LineWidth',lin_wid,'Color',color_model);
plot(data_for_modelfit_x,yvals_fit1+2*delta_fit1,'b--','LineWidth',lin_wid,'Color',color_model);
plot(data_for_modelfit_x,yvals_fit1-2*delta_fit1,'b--','LineWidth',lin_wid,'Color',color_model);
set(gca,...
    'xlim',[data_for_modelfit_x(1)-1 data_for_modelfit_x(end)+1],...
    'XTick',round(linspace(data_for_modelfit_x(1),data_for_modelfit_x(end),num_figs)),...
    'XTickLabel',polar_tits,...
    'ylim',y_lim,...
    'FontSize',font_tit-2);
title(['Polynomial = ' num2str(fit_poly1) ', {\itR}^2 = ' num2str(equ_gof1.rsquare) ', RMSE = ' num2str(equ_gof1.rmse)],'FontSize',font_tit,'FontWeight','normal');
xl = xlim; yl = ylim;
xt = 0.015 * (xl(2)-xl(1)) + xl(1); yt = 0.95 * (yl(2)-yl(1)) + yl(1);
caption = {['y = ' num2str(equ_fit1(1)) 'x + ' num2str(equ_fit1(2))]};
text(xt, yt, caption, 'FontSize', font_tit-1, 'Color', 'k', 'FontWeight', 'normal','FontAngle','italic');
ylabel('ITPC','FontSize',font_tit-2);

subplot(rows,cols,2);
hold on;
scatter(data_for_modelfit_x,data_for_modelfit_y,size_dots,color_dots,'.');
plot(data_for_modelfit_x,yvals_fit2,'b-','LineWidth',lin_wid,'Color',color_model);
plot(data_for_modelfit_x,yvals_fit2+2*delta_fit2,'b--','LineWidth',lin_wid,'Color',color_model);
plot(data_for_modelfit_x,yvals_fit2-2*delta_fit2,'b--','LineWidth',lin_wid,'Color',color_model);
set(gca,...
    'xlim',[data_for_modelfit_x(1)-1 data_for_modelfit_x(end)+1],...
    'XTick',round(linspace(data_for_modelfit_x(1),data_for_modelfit_x(end),num_figs)),...
    'XTickLabel',polar_tits,...
    'ylim',y_lim,...
    'FontSize',font_tit-2);
title(['Polynomial = ' num2str(fit_poly2) ', {\itR}^2 = ' num2str(equ_gof2.rsquare) ', RMSE = ' num2str(equ_gof2.rmse)],'FontSize',font_tit,'FontWeight','normal');
xl = xlim; yl = ylim;
xt = 0.015 * (xl(2)-xl(1)) + xl(1); yt = 0.95 * (yl(2)-yl(1)) + yl(1);
caption = {['y = ' num2str(equ_fit2(1)) 'x^3 ' num2str(equ_fit2(2)) 'x^2 + ', num2str(equ_fit2(3)) 'x + ' num2str(equ_fit2(4))]};
text(xt, yt, caption, 'FontSize', font_tit-1, 'Color', 'k', 'FontWeight', 'normal','FontAngle','italic');
ylabel('ITPC','FontSize',font_tit-2);
xlabel('Simulation phase offsets','FontSize',font_tit-2);

%% plot Knott ITPC data and related offsets
data_itpc_knott_devi = squeeze(mean(mean(mat_ITPC.devi(:,:,idx_hz2get_knott(1):idx_hz2get_knott(2),idx_time2get_knott(1):idx_time2get_knott(2)),4),3));
for clusti = 1:size(data_itpc_knott_devi,2)
    for parti = 1:size(data_itpc_knott_devi,1)
        tmp_data = data_itpc_knott_devi(parti,clusti);
        predicted_vals(parti,clusti) = equ_fit2(1)*(tmp_data.^3) +...
                                        equ_fit2(2)*(tmp_data.^2) +...
                                        equ_fit2(3)*(tmp_data) +...
                                        equ_fit2(4);
        clear tmp_data
    end
end
idx_CON = ismember(mat_demos.diagnosis,'CON');
idx_SCZ = ismember(mat_demos.diagnosis,'SCZ');
idx_MDD = ismember(mat_demos.diagnosis,'MDD');
rows=2;cols=1;
col_scz = [0, 0.4470, 0.7410]; alp = .3; lin_wid = 2;
font_tit = 14;
polar_tits = {'0','\pi/2','\pi','3\pi/2','2\pi'};

figure; set(gcf,'Position',[5 5 300 600]);
subplot(rows,cols,1);
tmp = [[data_itpc_knott_devi(idx_CON,1);NaN(2,1)], [data_itpc_knott_devi(idx_SCZ,1);NaN(3,1)], data_itpc_knott_devi(idx_MDD,1)];
CategoricalScatterplot(tmp,...
                    'Marker','x',...
                    'color',col_scz,...
                    'Labels',{'CON','SCZ','MDD'},...
                    'BoxColor',col_scz,...
                    'BoxAlpha',alp,...
                    'MedianLineWidth',lin_wid,...
                    'WhiskerLineWidth',lin_wid-.5,...
                    'WhiskerLineStyle','-');
set(gca,'xlim',[0.5 3.5],...
        'XTick',[1 2 3],...
        'XTickLabels',{'CON','SCZ','MDD'},...
        'ylim',[0 .75],...
        'FontSize',font_tit-2);
ylabel('ITPC','FontSize',font_tit-2);
if ranksum(tmp(:,1),tmp(:,2)) < 0.05
    sigstar({[1,2]},0.05);
end
if ranksum(tmp(:,1),tmp(:,3)) > 0.05
    sigstar({[1,3]},NaN);
end
clear tmp
title('ITPC in oddball deviants','FontSize',font_tit,'FontWeight','normal');

subplot(rows,cols,2);
tmp = [[predicted_vals(idx_CON,1);NaN(2,1)], [predicted_vals(idx_SCZ,1);NaN(3,1)], predicted_vals(idx_MDD,1)];
CategoricalScatterplot(tmp,...
                    'Marker','x',...
                    'color',col_scz,...
                    'Labels',{'CON','SCZ','MDD'},...
                    'BoxColor',col_scz,...
                    'BoxAlpha',alp,...
                    'MedianLineWidth',lin_wid,...
                    'WhiskerLineWidth',lin_wid-.5,...
                    'WhiskerLineStyle','-');
set(gca,'xlim',[0.5 3.5],...
        'XTick',[1 2 3],...
        'XTickLabels',{'CON','SCZ','MDD'},...
        'YTick',round(linspace(data_for_modelfit_x(1),data_for_modelfit_x(end),numel(polar_tits))),...
        'YTickLabel',polar_tits,...
        'FontSize',font_tit-2);
ylabel('Phase offsets','FontSize',font_tit-2);
if ranksum(tmp(:,1),tmp(:,2)) < 0.05
    sigstar({[1,2]},0.05);
end
if ranksum(tmp(:,1),tmp(:,3)) > 0.05
    sigstar({[1,3]},NaN);
end
clear tmp
title('Phase offset in oddball deviants','FontSize',font_tit,'FontWeight','normal');

