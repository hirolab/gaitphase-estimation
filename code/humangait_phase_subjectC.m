clear
clc
close all

r1 = [0.106518614215597;-0.0122284407286391;-0.000489674783200498];
Tw1 = [0.0190572790438554,-0.0545767194973006,-0.987141475961665;-0.818362894671529,0.581691358222422,-0.0755358656280325;0.516104786014129,0.794224031585750,-0.0441359452847375];
Ta1 = [0.0545548033438833,-0.0186956599204625,-1.11002005698209;-0.775734525784356,0.560856078014906,-0.0749197910089614;0.603477236118774,0.852064592621460,-0.0342370807324579];
Tm1 = quat2rotm(rotm2quat(Tw1 + Ta1));
abias1 = [-0.176027327106184,-0.563608249596711,0.217801989188737];
wbias1 = [-0.0167705082480907,-0.00215315618462306,0.00716710917440805];
mbias1 = [-20.1431866611659,-29.5314117701198,-81.7907600600365];

r2 = [0.0488994991459101;-0.0594425988846504;0.00755524297601401];
Tw2 = [0.114390765820317,0.135936910621720,0.998075336342434;-0.0215800304773111,-0.997666804733703,0.141703219077139;0.995839396095769,-0.0215156897407333,-0.0920046667481200];
Ta2 = [0.0973519875794760,0.105244916984141,0.998748638676139;-0.0241755650949058,-0.993503775752579,0.107609287678423;1.00148751474769,-0.0330445948445313,-0.103535894101711];
Tm2 = quat2rotm(rotm2quat(Tw2 + Ta2));
abias2 = [0.397151592074530,-0.157171147331356,0.697967428263995];
wbias2 = [0.00466370114304344,0.00250120673548035,0.00800608047758022];
mbias2 = [34.7924929530826,4.66658042753315,-52.8492739106742];

rfa = [-0.0115010448196404;0.00114083513064818;-0.000676346203193720];
rsa = [-0.0258406707164383;-0.206089578304228;-0.00286186105272751];

%% Load walking data
Fs = 400;

% outdoor ===============================

trial = load('\\nas01.itap.purdue.edu\puhome\Desktop\Gait_phase_paper\data\subjectC\subjectC_outdoor_walk_40min.mat');
t = trial.trial;
% pattern =================================================
load('\\nas01.itap.purdue.edu\puhome\Desktop\Gait_phase_paper\data\subjectC\gait_pattern_model.mat')




xcnt = 40;
coeffcnnt = 100;

% t = trial;
t_start_train =2*60;
t_end_train = 20.6;
t_start_test = 20*60;
t_end_test = 35.6;



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

patterns_hs = gait_pattern_model.patterns_hs;
patterns_to = gait_pattern_model.patterns_to;

% lag = linspace(-0.5, 0.5, size(patterns_hs,1))';
% figure, 
% h1 = subplot(211); hleg = plot_errbar(lag, mean(patterns_hs, 3), std(patterns_hs, [], 3)); title('HS pattern'), grid on
% h2 = subplot(212); plot_errbar(lag, mean(patterns_to, 3), std(patterns_to, [], 3)); title('TO pattern'), grid on, xlabel('lag [s]')
% legend(hleg, 'a_x', 'a_y', 'a_z', '\omega_x', '\omega_y', '\omega_z', 'NumColumns', 2)

lag = linspace(-0.5, 0.5, size(patterns_hs,1))';
figure, 
h1 = subplot(221); hleg = plot_errbar(lag, mean(patterns_hs(:,1:3,:), 3), std(patterns_hs(:,1:3,:), [], 3)); title('HS pattern  a^F_t'), grid on, xlabel('lag [s]'),ylabel('a [^{m}/_{s^2}]')
h2 = subplot(222); plot_errbar(lag, mean(patterns_to(:,1:3,:), 3), std(patterns_to(:,1:3,:), [], 3)); title('TO pattern  a^F_t'), grid on, xlabel('lag [s]'),ylabel('a [^{m}/_{s^2}]')
h3 = subplot(223); hleg = plot_errbar(lag, mean(patterns_hs(:,4:6,:), 3), std(patterns_hs(:,4:6,:), [], 3)); title('HS pattern  w^F_t'), grid on, xlabel('lag [s]'),ylabel('\omega [^{rad}/_{s}]')
h4 = subplot(224); plot_errbar(lag, mean(patterns_to(:,4:6,:), 3), std(patterns_to(:,4:6,:), [], 3)); title('TO pattern  w^F_t'), grid on, xlabel('lag [s]'),ylabel('\omega [^{rad}/_{s}]')
legend(hleg, 'x', 'y', 'z')


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

% 
% figure, set(gcf, 'position', [681   690   403   289])
% histogram(phase_error * 100 / (2*pi), 0:0.5:15, 'Normalization', 'count')
% xlabel('Gait Phase Error [%]'), xlim([0, 10]), grid on
% title('Histogram of residuals for training model')
% ylabel('Counts')
% set(gca,'box', 'off', 'tickdir', 'out', 'fontName', 'times new roman')
 
phase_error1 = phase_error;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% test trained model 

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


close all
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

