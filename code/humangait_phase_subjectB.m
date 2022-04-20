clear
clc
close all


r1 =  [0.099771490453043;0.001644277807084;-0.006819983125572];
Tw1 =[-0.188267794043455,-4.056656466959235e-04,-0.994832398620836;-0.799165389550038,0.570143135250826,0.151281373770168;0.565066521243301,0.818354996253341,-0.107707847557452] ;
Ta1 = [-0.173183493232567,0.010971855927482,-0.953441374111618;-0.801459466255059,0.572371604163400,0.148583141634766;0.563745060051088,0.825043045482814,-0.106720318759009] ;
Tm1 = [-0.038025666259945,-0.010640685987398,-0.999220108137944;-0.834348161843248,0.550629611294890,0.025887757609184;0.549924716242397,0.834681859732771,-0.029816094626279];
abias1 = [0.163034995364795,-0.319514697595614,0.384060644222848];
wbias1 = [-0.008524988350446,-0.014798800703154,0.004656544103471] ;
mbias1=[21.518394065515950,-78.216555392047010,-60.343380178362935] ;


r2 = [0.003141194702994;0.008742936179178;0.065303726601225] ;
Tw2 = [0.938513427049527,0.132427648925817,-0.352683776488376;-0.046500606347537,0.977474715290924,0.222231026412359;0.367361442322140,-0.189679755521854,0.909760791886977];
Ta2 = [0.917733069425707,0.127366495797751,-0.347154114005872;-0.046176405120172,0.933720786617610,0.211447930337766;0.381866785991804,-0.198655846266996,0.930203146901214] ;
Tm2 = [-0.422815826573027,-0.906116877366499,-0.013378391194486;0.794900729246127,-0.363749878168159,-0.485611837558173;0.435154693693339,-0.215958863407443,0.874072172004361] ;
abias2 = [0.339348048657677,-0.440280788147191,0.333386025322913] ; 
wbias2 = [-0.009835975316450,-0.002256048936303,-0.001425393269962] ;
mbias2 = [-16.814117394059450,-59.390506146596180,-81.411902621121290] ;

rfa = [0.004676125645037;-0.003701602901284;0.004218619980954] ;
rsa = [8.377455150248608e-04;-0.206418742789493;6.496951887002833e-05] ;

%% Load walking data
Fs = 400;

% outdoor ===============================
load('\\nas01.itap.purdue.edu\puhome\Desktop\Gait_phase_paper\data\subjectB\subjectB_outdoor walk 46min.mat')
% pattern =================================================
load('\\nas01.itap.purdue.edu\puhome\Desktop\Gait_phase_paper\data\subjectB\subjectB_gait_pattern_model.mat')

xcnt = 40;
coeffcnnt = 100;


t = trial;
t_start_train = 5*60;
t_end_train = 25;
t_start_test = 22*60;
t_end_test = 37.5;


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


figure, set(gcf, 'position', [681   690   403   289])
histogram(phase_error * 100 / (2*pi), 0:0.5:15, 'Normalization', 'count')
xlabel('Gait Phase Error [%]'), xlim([0, 10]), grid on
title('Histogram of residuals for training model')
ylabel('Counts')
set(gca,'box', 'off', 'tickdir', 'out', 'fontName', 'times new roman')

phase_error1 = phase_error;
figure, set(gcf, 'position', [681   472   741   507])
histogram(phase_error * 100 / (2*pi), 0:0.5:15, 'Normalization', 'count')
xlabel('Gait Phase Error [%]'), xlim([0, 10]), grid on
% title('Histogram of residuals for training model - $\overline{e}_{\gamma}$')
ylabel('Counts')
set(gca,'box', 'off', 'tickdir', 'out', 'fontName', 'Helvetica')
 
h = title(['\bf Histogram  of  residuals  for  training  model - ($\bf\overline{e}_{\gamma} =  $',' ',num2str(mean(phase_error * 100) / (2*pi),3), ' \%)']);
set(h,'interpreter','Latex','FontSize',12,'FontName','Helvetica')
 

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
