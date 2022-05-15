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



load('/Users/ahmysoli/Documents/GitHub/gaitphase-estimation/matlab code/subjectD_indoor_walk.mat')


t = trial;

t_tmp_i = (5<t.itime_tsync)&(t.itime_tsync<490);

t_tmp = t.itime_tsync(t_tmp_i);
time = (t_tmp(1) : 1/Fs :t_tmp(end))';

dlen = length(time);
wf = interp1(t_tmp, (t.w1(t_tmp_i,:) - wbias1) * inv(Tw1)', time);
af = interp1(t_tmp, (t.a1(t_tmp_i,:) - abias1) * inv(Ta1)', time);
mf = interp1(t_tmp, (t.m1(t_tmp_i,:) - mbias1) * Tm1, time);

ws = interp1(t_tmp, (t.w2(t_tmp_i,:) - wbias2) * inv(Tw2)', time);
as = interp1(t_tmp, (t.a2(t_tmp_i,:) - abias2) * inv(Ta2)', time);
ms = interp1(t_tmp, (t.m2(t_tmp_i,:) - mbias2) * Tm2, time);

% Label stance
warning('Change me if the subject changes - is_Stance')

wftol = 0.55;

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

% if isfield(t, 'mtime')  % has ground-truth
% 
%     missing_opt = {'spline', 'EndValues', 'nearest'};
%     fquat = quatnormalize(fillmissing(interp1(t.mtime, QTMParser.fill_missing_quat(t.fquat, missing_opt{:}), time), 'nearest'));
%     squat = quatnormalize(fillmissing(interp1(t.mtime, QTMParser.fill_missing_quat(t.squat, missing_opt{:}), time), 'nearest'));
%     s1 = fillmissing(interp1(t.mtime, fillmissing(t.ftrans + quatrotate(quatinv(t.fquat), r1'), missing_opt{:}), time), 'nearest');
%     s2 = fillmissing(interp1(t.mtime, fillmissing(t.strans + quatrotate(quatinv(t.squat), r2'), missing_opt{:}), time), 'nearest');
%     v1 = deriv_sgolay(s1, Fs, [3, 9]);
%     v2 = deriv_sgolay(s2, Fs, [3, 9]);
%     
% end

%%
x = zscore([af, wf]);
%% identify stance events using ground truth
pforce = interp1(t.mtime, t.pforce, time);

ftresh = 50;
hs_event = find(pforce(2:end,2) > ftresh & pforce(1:end-1,2) < ftresh);
to_event = find(pforce(2:end,2) < ftresh & pforce(1:end-1,2) > ftresh);

figure
% h1 = subplot(211); plot(time, af, time(hs_event), af(hs_event,1), 'ks',...
%  time(to_event), af(to_event), 'ko')
h1 = subplot(211); plot(time, af)
h2 = subplot(212); plot(time, pforce, time(hs_event), pforce(hs_event), 'ks', ...
    time(to_event), pforce(to_event), 'ko')
linkaxes([h1, h2], 'x')

figure
% h1 = subplot(211); plot(time, af, time(hs_event), af(hs_event,1), 'ks',...
%  time(to_event), af(to_event), 'ko')
h1 = subplot(223); plot(time, af); title ('Total a^F_t during indoor gait') 
ylabel('a [^{m}/_{s^2}]')
xlabel({'Time [s]';'(b)'});
grid on
xlim([0 300])
ylim([-50 50])

h2 = subplot(221); plot(time, pforce),hold on
plot(time(hs_event),pforce(hs_event), 'ks', ...
    time(to_event), pforce(to_event), 'ko','LineWidth',0.8); title('Total F_N during indoor gait')
xlabel({'Time [s]';'(a)'});
xlim([0 300])
ylim([-100 900])
linkaxes([h1, h2], 'x')

ylabel('Force [N]')
grid on

h3 = subplot(224); plot(time, af),hold on
ylabel('a [^{m}/_{s^2}]')
xlabel({'Time [s]';'(d)'});
grid on

plot(time(hs_event), af(hs_event,1), 'ks', ...
    time(to_event), af(to_event,1), 'ko','LineWidth',1.5);title ('Segmented a^F_t during indoor gait');xlim([135 141]);
ylim([-30 40])


h4 = subplot(222); plot(time, pforce),hold on
plot(time(hs_event), pforce(hs_event), 'ks', ...
    time(to_event), pforce(to_event), 'ko','LineWidth',1.5); title('Segmented F_N during indoor gait');ylim([-100 750]);
xlim([135 141]);
xlabel({'Time [s]';'(c)'});

legend('','','','HS','TO','location','best','Fontsize',9)
linkaxes([h3, h4], 'x')
ylabel('Force [N]')
grid on



trows = (round(-0.5*Fs) : round(0.5*Fs))';
% size(x(hs_event' + trows,:))
patterns_hs = permute(reshape(x(hs_event' + trows,:), length(trows), length(hs_event), []), [1, 3, 2]);
patterns_to = permute(reshape(x(to_event' + trows,:), length(trows), length(to_event), []), [1, 3, 2]);

return
%% Identifying HS and TO && saving patttern
gait_pattern_model = struct('patterns_hs', patterns_hs, 'patterns_to', patterns_to, ...
    'readme', 'Input: zscore([vaf, aaf, wfh]) around 401 sample window, 400 Hz');

uisave('gait_pattern_model','gait_pattern_model')


    