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



load('\\nas01.itap.purdue.edu\puhome\Desktop\Gait_phase_paper\data\subjectA\subjectA_indoor_walk.mat')


t = trial(2);

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

if isfield(t, 'mtime')  % has ground-truth

    missing_opt = {'spline', 'EndValues', 'nearest'};
    fquat = quatnormalize(fillmissing(interp1(t.mtime, QTMParser.fill_missing_quat(t.fquat, missing_opt{:}), time), 'nearest'));
    squat = quatnormalize(fillmissing(interp1(t.mtime, QTMParser.fill_missing_quat(t.squat, missing_opt{:}), time), 'nearest'));
    s1 = fillmissing(interp1(t.mtime, fillmissing(t.ftrans + quatrotate(quatinv(t.fquat), r1'), missing_opt{:}), time), 'nearest');
    s2 = fillmissing(interp1(t.mtime, fillmissing(t.strans + quatrotate(quatinv(t.squat), r2'), missing_opt{:}), time), 'nearest');
    v1 = deriv_sgolay(s1, Fs, [3, 9]);
    v2 = deriv_sgolay(s2, Fs, [3, 9]);
    
end

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
xlim([0 500])

h2 = subplot(221); plot(time, pforce),hold on
plot(time(hs_event),pforce(hs_event), 'ks', ...
    time(to_event), pforce(to_event), 'ko','LineWidth',0.8); title('Total F_N during indoor gait')
xlabel({'Time [s]';'(a)'});
xlim([0 500])

linkaxes([h1, h2], 'x')

ylabel('Force [N]')
grid on

h3 = subplot(224); plot(time, af),hold on
ylabel('a [^{m}/_{s^2}]')
xlabel({'Time [s]';'(d)'});
grid on

plot(time(hs_event), af(hs_event,1), 'ks', ...
    time(to_event), af(to_event,1), 'ko','LineWidth',1.5);title ('Segmented a^F_t during indoor gait');xlim([224 232]);

h4 = subplot(222); plot(time, pforce),hold on
plot(time(hs_event), pforce(hs_event), 'ks', ...
    time(to_event), pforce(to_event), 'ko','LineWidth',1.5); title('Segmented F_N during indoor gait')
xlim([224 232]); 
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


    