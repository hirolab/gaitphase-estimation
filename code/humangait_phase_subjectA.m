clear
clc
close all


r1 = [0.091305262300364;-0.006017738044598;0.005762894304378] ;
Tw1 = [-0.035766217056836,0.224271496587231,-0.986855483131145;-0.799666833943993,0.575575075329483,0.154996430729911;0.598636436256773,0.788519931772248,0.164842131192030];
Ta1 = [-0.015533192510838,0.227452625653973,-0.938506470036921;-0.803033431853574,0.576353036226075,0.151062173829026;0.595436706713049,0.791674638072338,0.170923221853200];
Tm1 = [-0.038025666259945,-0.010640685987398,-0.999220108137944;-0.834348161843248,0.550629611294890,0.025887757609184;0.549924716242397,0.834681859732771,-0.029816094626279];
abias1 = [0.202534697017711,-0.335720616841910,0.380142380236680];
wbias1 =  [-0.006191581347783,-0.011518655999633,0.003261243232253];
mbias1= [-34.537674676025050,-86.271648812013200,-73.025554005125580];


r2 = [-0.012042146113055;-0.020933453618862;-0.040687037093892] ;
Tw2 = [0.026609109556093,1.007641001818010,-0.101724300408771;0.995034910208084,-0.016939743254783,0.033188356892398;0.028616620945591,-0.096837535901401,-1.007774892110355];
Ta2 =  [0.028666848796954,0.998610258511563,-0.105271165939392;0.963169622218637,-0.023470355157086,0.027453469222973;0.039544193902019,-0.108601850486066,-1.024785151729591];
Tm2 = [-0.422815826573027,-0.906116877366499,-0.013378391194486;0.794900729246127,-0.363749878168159,-0.485611837558173;0.435154693693339,-0.215958863407443,0.874072172004361] ;
abias2 = [0.261940921219439,-0.323127819335090,0.413641352368841]; 
wbias2 =  [-0.009163141961407,-6.665360925654811e-04,-0.002722902654606] ;
mbias2 =  [-47.437304022632200,-12.820750403930155,-71.031656346552720];

rfa =  [-0.004004689514154;-0.006964576615733;0.009818576708270];
rsa =  [3.538934950412793e-04;-0.212277798170128;0.001200695872606];



%% Load walking data
Fs = 400;

% outdoor ===============================
load('\\nas01.itap.purdue.edu\puhome\Desktop\Gait_phase_paper\data\subjectA\subjectA_outdoor walk 40min.mat')

% pattern =================================================
load('\\nas01.itap.purdue.edu\puhome\Desktop\Gait_phase_paper\data\subjectA\subjectA_gait_pattern_model.mat')


xcnt = 40;
coeffcnnt = 100;


t = trial;
t_start_train = 1*60;
t_end_train = 20;
t_start_test =20*60;
t_end_test = 40;


time = (t_start_train : 1/Fs : t_end_train*60)';

dlen = length(time);

wf = interp1(t.itime, (t.w1 - wbias1) * inv(Tw1)', time);
af = interp1(t.itime, (t.a1 - abias1) * inv(Ta1)', time);
mf = interp1(t.itime, (t.m1 - mbias1) * Tm1, time);

ws = interp1(t.itime, (t.w2 - wbias2) * inv(Tw2)', time);
as = interp1(t.itime, (t.a2 - abias2) * inv(Ta2)', time);
ms = interp1(t.itime, (t.m2 - mbias2) * Tm2, time);

% Label stance
warning('Change me if the subject changes - is_Stance')
wftol = .55;
aftol = 1;
pitchtol = 30 * pi/180;

yf = af ./ sqrt(sum(af.^2, 2));

is_stance = (sqrt(sum(wf.^2, 2)) < wftol) & ...
    (abs(sqrt(sum(af.^2, 2)) - 9.81) < aftol) & ...
    (abs(yf(:,2)) > cos(pitchtol));

figure, 
h1 = subplot(221); plot(time, wf, time, (~is_stance)*5), grid on, title('\omega_f')
h2 = subplot(222); plot(time, af, time, (~is_stance)*5), grid on, title('a_f')
h3 = subplot(223); plot(time, ws), grid on, title('\omega_s')
h4 = subplot(224); plot(time, as), grid on, title('a_s')
linkaxes([h1, h2, h3, h4], 'x')


%% create features
x = zscore([af, wf]);

%% Identifying HS and TO using loaded pattern


patterns_hs = gait_pattern_model.patterns_hs;
patterns_to = gait_pattern_model.patterns_to;

lag = linspace(-0.5, 0.5, size(patterns_hs,1))';

figure, 
h1 = subplot(221); hleg = plot_errbar(lag, mean(patterns_hs(:,1:3,:), 3), std(patterns_hs(:,1:3,:), [], 3)); title('HS pattern 1:3'), grid on
h2 = subplot(222); plot_errbar(lag, mean(patterns_to(:,1:3,:), 3), std(patterns_to(:,1:3,:), [], 3)); title('TO pattern 1:3'), grid on, xlabel('lag [s]')
h3 = subplot(223); hleg = plot_errbar(lag, mean(patterns_hs(:,4:6,:), 3), std(patterns_hs(:,4:6,:), [], 3)); title('HS pattern 4:6'), grid on
h4 = subplot(224); plot_errbar(lag, mean(patterns_to(:,4:6,:), 3), std(patterns_to(:,4:6,:), [], 3)); title('TO pattern 4:6'), grid on, xlabel('lag [s]')
legend(hleg, 'x', 'y', 'z')

ncha = size(patterns_hs);
dlen = size(x, 1);

cor_hs = zeros(dlen, ncha(3));
cor_to = zeros(dlen, ncha(3));
for i = 1 : ncha(3)
    aa = xcorr2(x, patterns_hs(:,:,i));
    cor_hs(:,i) = sum(aa(floor(ncha(1)/2) + (1:dlen),floor(ncha(2)/2) + (1:ncha(2))), 2);
    
    aa = xcorr2(x, patterns_to(:,:,i));
    cor_to(:,i) = sum(aa(floor(ncha(1)/2) + (1:dlen),floor(ncha(2)/2) + (1:ncha(2))), 2);
end

cor_hs = mean(cor_hs, 2);
cor_to = mean(cor_to, 2);

[~, hs_event] = findpeaks(cor_hs, 'MinPeakDistance', 0.8*Fs, 'MinPeakHeight', prctile(cor_hs, 75));
[~, to_event] = findpeaks(cor_to, 'MinPeakDistance', 0.8*Fs, 'MinPeakHeight', prctile(cor_to, 75));

hs_event([1:3,end-2:end]) = [];
to_event([1:3,end-2:end]) = [];

figure
h1 = subplot(211); plot(time, af, time(hs_event), af(hs_event), 'rs', ...
    time(to_event), af(to_event), 'bo')
h2 = subplot(212); plot(time, cor_hs, 'r', time, cor_to, 'b'), grid on
legend('HS confidence', 'TO confidence')
linkaxes([h1, h2], 'x')

%% Select valid stance phases

slen_min = 0.5;
slen_max = 1;

glen_min = 0.8;
glen_max = 1.5;

count = 0;
walk = repmat(struct('idx', 0, 'slen', 0, 'glen', 0), [length(hs_event)-1, 1]);

for i = 1 : length(hs_event)-1
    
    i_to = find(to_event > hs_event(i) + slen_min*Fs, 1, 'first');
    i_hs2 = find(hs_event > hs_event(i) + glen_min*Fs, 1, 'first');
    
    if ~isempty(i_to) && ~ isempty(i_hs2) && ...
       (to_event(i_to) - hs_event(i) < slen_max*Fs) && ...
       (hs_event(i_hs2) - hs_event(i) < glen_max*Fs)
        
        count = count + 1;
        walk(count).idx = hs_event(i);
        walk(count).slen = (to_event(i_to) - hs_event(i)) / Fs;
        walk(count).glen = (hs_event(i_hs2) - hs_event(i)) / Fs;
    end
end
walk(count + 1 : end) = [];

%% Create predictor and label matrices
tlen = 1.0;  % estimation window
wrows = (round(-tlen*Fs) : 0)';
wlen = length(wrows);

features = zscore([as, ws]);
nx = size(features, 2);

N = round(glen_max * Fs * length(walk));
prows = 0;
X = zeros(N, wlen*nx);
Y = zeros(N, 2);

for i = 1 : length(walk)
    rows = walk(i).idx + (0 : round(walk(i).glen*Fs));
    
    phase = linspace(0, 2*pi, length(rows))'; % ground truth calculated using interpolation based on detected heel strikes
    
    xx = permute(reshape(features(rows + wrows,:), wlen, [], nx), [1, 3, 2]);
    
    prows = prows(end) + (1:length(rows));
    X(prows,:) = reshape(xx, [], length(rows))';
    Y(prows,:) = [sin(phase), cos(phase)]; 
    
    if mod(i, round(length(walk)/10)) == 0
        fprintf('%d%%, ', round(100*i/length(walk)))
    end
end

X(prows(end)+1:end,:) = [];
Y(prows(end)+1:end,:) = [];
%% train model 


[coeff,score,latent,~,] = pca(X(1:xcnt:end,:));
Xavg = mean(X);

Xdown = (X - Xavg) * coeff(:,1:coeffcnnt);
fit_phase = Xdown \ Y;
Yh = Xdown * fit_phase;
Yh = Yh ./ sqrt(sum(Yh.^2, 2));

phase_error = acos(sum(Y .* Yh, 2));


phase_error1 = phase_error;


%% test trained model 

time = (t_start_test : 1/Fs : t_end_test*60)';

dlen = length(time);

wf = interp1(t.itime, (t.w1 - wbias1) * inv(Tw1)', time);
af = interp1(t.itime, (t.a1 - abias1) * inv(Ta1)', time);
mf = interp1(t.itime, (t.m1 - mbias1) * Tm1, time);

ws = interp1(t.itime, (t.w2 - wbias2) * inv(Tw2)', time);
as = interp1(t.itime, (t.a2 - abias2) * inv(Ta2)', time);
ms = interp1(t.itime, (t.m2 - mbias2) * Tm2, time);

% Label stance
warning('Change me if the subject changes - is_Stance')
wftol = .55;
aftol = 1;
pitchtol = 30 * pi/180;

yf = af ./ sqrt(sum(af.^2, 2));

is_stance = (sqrt(sum(wf.^2, 2)) < wftol) & ...
    (abs(sqrt(sum(af.^2, 2)) - 9.81) < aftol) & ...
    (abs(yf(:,2)) > cos(pitchtol));

figure, 
h1 = subplot(221); plot(time, wf, time, (~is_stance)*5), grid on, title('\omega_f')
h2 = subplot(222); plot(time, af, time, (~is_stance)*5), grid on, title('a_f')
h3 = subplot(223); plot(time, ws), grid on, title('\omega_s')
h4 = subplot(224); plot(time, as), grid on, title('a_s')
linkaxes([h1, h2, h3, h4], 'x')


%% create features
x = zscore([af, wf]);

%% Identifying HS and TO using loaded pattern


patterns_hs = gait_pattern_model.patterns_hs;
patterns_to = gait_pattern_model.patterns_to;

lag = linspace(-0.5, 0.5, size(patterns_hs,1))';

figure, 
h1 = subplot(221); hleg = plot_errbar(lag, mean(patterns_hs(:,1:3,:), 3), std(patterns_hs(:,1:3,:), [], 3)); title('HS pattern 1:3'), grid on
h2 = subplot(222); plot_errbar(lag, mean(patterns_to(:,1:3,:), 3), std(patterns_to(:,1:3,:), [], 3)); title('TO pattern 1:3'), grid on, xlabel('lag [s]')
h3 = subplot(223); hleg = plot_errbar(lag, mean(patterns_hs(:,4:6,:), 3), std(patterns_hs(:,4:6,:), [], 3)); title('HS pattern 4:6'), grid on
h4 = subplot(224); plot_errbar(lag, mean(patterns_to(:,4:6,:), 3), std(patterns_to(:,4:6,:), [], 3)); title('TO pattern 4:6'), grid on, xlabel('lag [s]')
legend(hleg, 'x', 'y', 'z')

ncha = size(patterns_hs);
dlen = size(x, 1);

cor_hs = zeros(dlen, ncha(3));
cor_to = zeros(dlen, ncha(3));
for i = 1 : ncha(3)
    aa = xcorr2(x, patterns_hs(:,:,i));
    cor_hs(:,i) = sum(aa(floor(ncha(1)/2) + (1:dlen),floor(ncha(2)/2) + (1:ncha(2))), 2);
    
    aa = xcorr2(x, patterns_to(:,:,i));
    cor_to(:,i) = sum(aa(floor(ncha(1)/2) + (1:dlen),floor(ncha(2)/2) + (1:ncha(2))), 2);
end

cor_hs = mean(cor_hs, 2);
cor_to = mean(cor_to, 2);

[~, hs_event] = findpeaks(cor_hs, 'MinPeakDistance', 0.8*Fs, 'MinPeakHeight', prctile(cor_hs, 75));
[~, to_event] = findpeaks(cor_to, 'MinPeakDistance', 0.8*Fs, 'MinPeakHeight', prctile(cor_to, 75));

hs_event([1:3,end-2:end]) = [];
to_event([1:3,end-2:end]) = [];

figure
h1 = subplot(211); plot(time, af, time(hs_event), af(hs_event), 'rs', ...
    time(to_event), af(to_event), 'bo')
h2 = subplot(212); plot(time, cor_hs, 'r', time, cor_to, 'b'), grid on
legend('HS confidence', 'TO confidence')
linkaxes([h1, h2], 'x')

%% Select valid stance phases

slen_min = 0.5;
slen_max = 1;

glen_min = 0.8;
glen_max = 1.5;

count = 0;
walk = repmat(struct('idx', 0, 'slen', 0, 'glen', 0), [length(hs_event)-1, 1]);

for i = 1 : length(hs_event)-1
    
    i_to = find(to_event > hs_event(i) + slen_min*Fs, 1, 'first');
    i_hs2 = find(hs_event > hs_event(i) + glen_min*Fs, 1, 'first');
    
    if ~isempty(i_to) && ~ isempty(i_hs2) && ...
       (to_event(i_to) - hs_event(i) < slen_max*Fs) && ...
       (hs_event(i_hs2) - hs_event(i) < glen_max*Fs)
        
        count = count + 1;
        walk(count).idx = hs_event(i);
        walk(count).slen = (to_event(i_to) - hs_event(i)) / Fs;
        walk(count).glen = (hs_event(i_hs2) - hs_event(i)) / Fs;
    end
end
walk(count + 1 : end) = [];

%% Create predictor and label matrices
tlen = 1.0;  % estimation window
wrows = (round(-tlen*Fs) : 0)';
wlen = length(wrows);

features = zscore([as, ws]);
nx = size(features, 2);

N = round(glen_max * Fs * length(walk));
prows = 0;
X = zeros(N, wlen*nx);
Y = zeros(N, 2);

for i = 1 : length(walk)
    rows = walk(i).idx + (0 : round(walk(i).glen*Fs));
    
    phase = linspace(0, 2*pi, length(rows))'; % ground truth calculated using interpolation based on detected heel strikes
    
    xx = permute(reshape(features(rows + wrows,:), wlen, [], nx), [1, 3, 2]);
    
    prows = prows(end) + (1:length(rows));
    X(prows,:) = reshape(xx, [], length(rows))';
    Y(prows,:) = [sin(phase), cos(phase)]; 
    
    if mod(i, round(length(walk)/10)) == 0
        fprintf('%d%%, ', round(100*i/length(walk)))
    end
end

X(prows(end)+1:end,:) = [];
Y(prows(end)+1:end,:) = [];
%% test model 

[coeff,score,latent,~,] = pca(X(1:xcnt:end,:));
Xavg = mean(X);
Xdown = (X - Xavg) * coeff(:,1:coeffcnnt);

%% 
Yh = Xdown * fit_phase;
Yh = Yh ./ sqrt(sum(Yh.^2, 2));

phase_error = acos(sum(Y .* Yh, 2));


figure, set(gcf, 'position', [681   690   403   289])
histogram(phase_error * 100 / (2*pi), 0:0.5:15, 'Normalization', 'count')
xlabel('Gait Phase Error [%]'), xlim([0, 10]), grid on
title('Histogram of residuals for testing')
ylabel('Counts')
set(gca,'box', 'off', 'tickdir', 'out', 'fontName', 'times new roman')
accuracy_err = mean(phase_error * 100 / (2*pi));

% close all
figure, set(gcf, 'position', [681   472   741   507])
subplot(121)
histogram(phase_error1 * 100 / (2*pi), 0:0.5:15, 'Normalization', 'count')
xlabel({'Gait Phase Error [%]';'(a)'}), xlim([0, 10]), grid on
% title('Histogram of residuals for training model - $\overline{e}_{\gamma}$')
ylabel('Counts')
set(gca,'box', 'off', 'tickdir', 'out', 'fontName', 'Helvetica')
 
h = title(['\bf Histogram  of  residuals  for  training  model- ($\bf\overline{e}_{\gamma} =  $',' ',num2str(mean(phase_error1 * 100) / (2*pi),3), ' \%)']);
set(h,'interpreter','Latex','FontSize',12,'FontName','Helvetica')

subplot(122)
histogram(phase_error * 100 / (2*pi), 0:0.5:15, 'Normalization', 'count')
xlabel({'Gait Phase Error [%]';'(b)'}), xlim([0, 10]), grid on
% xlabel({'Time [s]';'(b)'});

% title('Histogram of residuals for training model - $\overline{e}_{\gamma}$')
ylabel('Counts')
set(gca,'box', 'off', 'tickdir', 'out', 'fontName', 'Helvetica')
 
h = title(['\bf Histogram  of  residuals  for  testing  - ($\bf\overline{e}_{\gamma} =  $',' ',num2str(mean(phase_error * 100) / (2*pi),3), ' \%)']);
set(h,'interpreter','Latex','FontSize',12,'FontName','Helvetica')
return



