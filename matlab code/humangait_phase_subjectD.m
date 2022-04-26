clear
clc
close all

r1 = [0.095026572520900;-3.945808440737073e-04;-0.012215536575702];
Tw1 = [0.071678660788978,0.077991073961893,-1.011845839167415;-0.841257915321978,0.545824935665565,-0.022535977381557;0.534292226467099,0.834695938230340,0.099677944622817];
Ta1 = [0.083024044472027,0.079245218750652,-0.976260362755186;-0.862815564353740,0.551876139545587,-0.018974253862459;0.562318899193943,0.850249885305265,0.099597937159941];
Tm1 = [0.083024044472027,0.079245218750652,-0.976260362755186;-0.862815564353740,0.551876139545587,-0.018974253862459;0.562318899193943,0.850249885305265,0.099597937159941];
abias1 = [0.239327370166102,-0.487187326511065,0.373706818480228];
wbias1 = [-0.004364524050998,-0.013253968537931,0.006623539902981];
mbias1=[0,0,0];


r2 = [0.006947997609772;0.046648126316578;-0.043206192549542];
Tw2 =[-0.047910149029392,1.017401862435155,0.027186834466638;0.749963900591967,0.039227227320776,0.679285821406879;0.680586284438009,0.042417450889460,-0.750127207749930];
Ta2 = [-0.053597568015329,0.995989037638687,0.003731795514654;0.695952627563320,0.020606177312225,0.656792529342617;0.707965990888661,0.042333129485294,-0.769656806401236];
Tm2 = [-0.422815826573027,-0.906116877366499,-0.013378391194486;0.794900729246127,-0.363749878168159,-0.485611837558173;0.435154693693339,-0.215958863407443,0.874072172004361];
abias2 = [0.332285039071802,-0.217533095148346,0.493587513826306];
wbias2 = [-0.006498627843007,-0.001556319920726,-9.239718203433528e-04];
mbias2 = [0,0,0];

rfa = [0.007285503939759;-8.741828959805633e-04;-0.003253437310353];
rsa = [0.003086842952211;-0.199674937857264;0.001418152320710];


%% Load walking data
Fs = 400;

% outdoor ===============================
load('C:\Users\solimana\Documents\gaitphase-estimation\data\subjectD\subjectD_outdoor_walk_40min.mat')

% pattern =================================================
load('C:\Users\solimana\Documents\gaitphase-estimation\data\subjectD\subjectD_gait_pattern_model.mat')

xcnt = 40;
coeffcnnt = 100;


t = trial;
t_start_train =1*60;
t_end_train = 19.6;
t_start_test = 21*60;
t_end_test = 39.6;


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


figure, 
sgtitle('Total Outdoor walk data','fontsize',10,'fontweight','bold')
h5 = subplot(313);plot(time, wf), grid on, title('Segment of outdoor walk (\omega_f)'),xlim([613 615]);
ylabel('\omega [^{rad}/_{s}]')
xlabel('Time [s]');
grid on

h1 = subplot(321); plot(time, wf), grid on, title('\omega_f')
ylabel('\omega [^{rad}/_{s}]')
% xlabel('Time [s]');
grid on
h2 = subplot(322); plot(time, af), grid on, title('a_f')
ylabel('a [^{m}/_{s^2}]')
% xlabel('Time [s]');
grid on

h3 = subplot(323); plot(time, ws), grid on, title('\omega_s')
ylabel('\omega [^{rad}/_{s}]')
% xlabel('Time [s]');
grid on
h4 = subplot(324); plot(time, as), grid on, title('a_s')
ylabel('a [^{m}/_{s^2}]')
ylabel('a [^{m}/_{s^2}]')
% xlabel('Time [s]');
grid on

% xlabel('Time [s]');
grid on
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

% figure, 
% h1 = subplot(221); hleg = plot_errbar(lag, mean(patterns_hs(:,1:3,:), 3), std(patterns_hs(:,1:3,:), [], 3)); title('HS pattern  a^F_t'), grid on, xlabel('lag [s]'),ylabel('a [^{m}/_{s^2}]')
% h2 = subplot(222); plot_errbar(lag, mean(patterns_to(:,1:3,:), 3), std(patterns_to(:,1:3,:), [], 3)); title('TO pattern  a^F_t'), grid on, xlabel('lag [s]'),ylabel('a [^{m}/_{s^2}]')
% h3 = subplot(223); hleg = plot_errbar(lag, mean(patterns_hs(:,4:6,:), 3), std(patterns_hs(:,4:6,:), [], 3)); title('HS pattern  w^F_t'), grid on, xlabel('lag [s]'),ylabel('\omega [^{m}/_{s^2}]')
% h4 = subplot(224); plot_errbar(lag, mean(patterns_to(:,4:6,:), 3), std(patterns_to(:,4:6,:), [], 3)); title('TO pattern  w^F_t'), grid on, xlabel('lag [s]'),ylabel('\omega [^{m}/_{s^2}]')
% legend(hleg, 'x', 'y', 'z')

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



figure
h1 = subplot(221);plot(time, af),hold on
ylabel('a [^{m}/_{s^2}]')
xlabel('Time [s]');
grid on
plot(time(hs_event), af(hs_event,1), 'ks', ...
   time(to_event), af(to_event,1), 'ko','LineWidth',1.5);title ('a^F_t during outdoor gait'),xlim([987.5 990]);
legend('','','','HS','TO','location','best','Fontsize',9)

h2 = subplot(223); plot(time, cor_hs, 'r', time, cor_to, 'b'), grid on,xlim([987.5 990]);
legend('HS confidence', 'TO confidence','location','best','Fontsize',9),ylabel('Correlation confidence');
xlabel('Time [s]');

h3 = subplot(222);plot(time, wf),hold on
ylabel('\omega [^{rad}/_{s}]')
xlabel('Time [s]');
grid on
plot(time(hs_event), wf(hs_event,1), 'ks', ...
   time(to_event), wf(to_event,1), 'ko','LineWidth',1.5);title ('\omega ^F_t during outdoor gait'),xlim([800 802.5]);
% legend('','','','HS','TO','location','best','Fontsize',9)

h4 = subplot(224); plot(time, cor_hs, 'r', time, cor_to, 'b'), grid on,xlim([800 802.5]);
% legend('HS confidence', 'TO confidence','location','best','Fontsize',9),
ylabel('Correlation confidence');
xlabel('Time [s]');

linkaxes([h1, h2], 'x')
linkaxes([h3, h4], 'x')




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
close all
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
ts = linspace(0,1000/Fs,1000)';
figure
subplot(211),plot(ts,Y(2000:3000-1,:),'LineWidth',1.2),xlim([0 1000/400]);
legend('sin(\theta)' ,'cos(\theta)','Location','best')
ylabel('\gamma')
xlabel({'Time [s]';'(a)'});
title ('\gamma during outdoor gait')
grid on

subplot(223),plot(ts,features(2000:3000-1,1:3),'LineWidth',1.2),xlim([0 1000/400]);
legend('x' ,'y','z','Location','best')
ylabel('a [^{m}/_{s^2}]')
xlabel({'Time [s]';'(b)'});
title ('Z-score of a ^S_t during outdoor gait')
grid on

subplot(224),plot(ts,features(2000:3000-1,4:6),'LineWidth',1.2),xlim([0 1000/400]);
% legend('x' ,'y','z','Location','best')
ylabel('\omega [^{rad}/_{s}]')
xlabel({'Time [s]';'(c)'});
title ('Z-score of  \omega ^S_t during outdoor gait')
grid on

%% train model 


[coeff,score,latent,~,] = pca(X(1:xcnt:end,:));
Xavg = mean(X);

Xdown = (X - Xavg) * coeff(:,1:coeffcnnt);

% Xdown = (X - Xavg) ;

fit_phase = Xdown \ Y;
Yh = Xdown * fit_phase;
% Yh1 = Yh;
% size(Yh)
Yh = Yh ./ sqrt(sum(Yh.^2, 2));
% size(Yh)

phase_error = acos(sum(Y .* Yh, 2));

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


% figure, set(gcf, 'position', [681   690   403   289])
% histogram(phase_error * 100 / (2*pi), 0:0.5:15, 'Normalization', 'count')
% xlabel('Gait Phase Error [%]'), xlim([0, 10]), grid on
% title('Histogram of residuals for testing')
% ylabel('Counts')
% set(gca,'box', 'off', 'tickdir', 'out', 'fontName', 'times new roman')
% accuracy_err = mean(phase_error * 100 / (2*pi))

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


figure, set(gcf, 'position', [681   472   741   507])
subplot(121)
histogram(phase_error1 * 100 / (2*pi), 0:0.5:15, 'Normalization', 'count')
xlabel({'Gait Phase Error [%]';'(a)'},'Fontsize',8), xlim([0, 7]), grid on
% title('Histogram of residuals for training model - $\overline{e}_{\gamma}$')
ylabel('Counts','Fontsize',8)
set(gca,'box', 'off', 'tickdir', 'out', 'fontName', 'Helvetica')
 
h = title(['\bf Training  dataset - ($\bf\overline{e}_{\gamma} =  $',' ',num2str(mean(phase_error1 * 100) / (2*pi),3), ' \%)']);
set(h,'interpreter','Latex','FontSize',12,'FontName','Helvetica')

subplot(122)
histogram(phase_error * 100 / (2*pi), 0:0.5:15, 'Normalization', 'count')
xlabel({'Gait Phase Error [%]';'(b)'},'Fontsize',8), xlim([0, 7]), grid on
% xlabel({'Time [s]';'(b)'});

% title('Histogram of residuals for training model - $\overline{e}_{\gamma}$')
ylabel('Counts','Fontsize',8)
set(gca,'box', 'off', 'tickdir', 'out', 'fontName', 'Helvetica')
 
h = title(['\bf Validation dataset  - ($\bf\overline{e}_{\gamma} =  $',' ',num2str(mean(phase_error * 100) / (2*pi),3), ' \%)']);
set(h,'interpreter','Latex','FontSize',12,'FontName','Helvetica')
return

