close all;
clear all;
clc;

load('stall_duration_and_resumption_mark_DNAP_100nM_Dct_helicase_45nM_231208_ch5.mat')
load('stall_duration_and_resumption_mark_DNAP_100nM_Dct_helicase_45nM_231208_ch6.mat')
load('stall_duration_and_resumption_mark_DNAP_100nM_Dct_helicase_45nM_240227_ch4.mat')
load('stall_duration_and_resumption_mark_DNAP_100nM_Dct_helicase_45nM_240227_ch5.mat')
load('stall_duration_and_resumption_mark_DNAP_100nM_Dct_helicase_45nM_240227_ch6.mat')

load('stall_duration_and_resumption_mark_DNAP_100nM_Dct_helicase_45nM_240320_ch1.mat')
load('stall_duration_and_resumption_mark_DNAP_100nM_Dct_helicase_45nM_240320_ch2.mat')


%%
stall_title = '100nM DNAP + 45nM Dct helicase, unwind -100 turns, stall duration';

% choose stall bining
stall_duration = [stall_duration_DNAP_100nM_Dct_helicase_45nM_231208_ch5, stall_duration_DNAP_100nM_Dct_helicase_45nM_231208_ch6, stall_duration_DNAP_100nM_Dct_helicase_45nM_240227_ch4, stall_duration_DNAP_100nM_Dct_helicase_45nM_240227_ch5, stall_duration_DNAP_100nM_Dct_helicase_45nM_240227_ch6, stall_duration_DNAP_100nM_Dct_helicase_45nM_240320_ch1, stall_duration_DNAP_100nM_Dct_helicase_45nM_240320_ch2];
resumption_recovery_duration = [resum_recovery_duration_DNAP_100nM_Dct_helicase_45nM_231208_ch5, resum_recovery_duration_DNAP_100nM_Dct_helicase_45nM_231208_ch6, resum_recovery_duration_DNAP_100nM_Dct_helicase_45nM_240227_ch4, resum_recovery_duration_DNAP_100nM_Dct_helicase_45nM_240227_ch5, resum_recovery_duration_DNAP_100nM_Dct_helicase_45nM_240227_ch6, resum_recov_duration_DNAP_100nM_Dct_helicase_45nM_240320_ch1, resum_recov_duration_DNAP_100nM_Dct_helicase_45nM_240320_ch2];

% sort the stall duration and adjust the resumption recovery as well
[stall_duration_sorted, sortIdx] = sort(stall_duration);
resumption_recovery_duration_selected = resumption_recovery_duration(sortIdx);

stall_duration_bin = stall_duration_sorted(64:96); % 
resumption_recovery_bin = resumption_recovery_duration_selected(64:96); % 

%% plot
fig_handle = figure;
pos = [200 200 1200 550];
set(fig_handle, 'Pos', pos);
h = gca;
h.XAxis.Visible = 'off';
box(h,'off');
set(h,'FontSize',13,'LineWidth',1.5);
            
            
subplot (1,2,1)
h1 = histfit(stall_duration_bin,5)
h1(1).FaceColor = [.8 .8 1];
pd1 = fitdist(stall_duration_bin','Normal')
text(220, 4, ['N=' num2str(length(stall_duration_bin))], 'FontSize',14) 
text(220, 3, ['Mean:' num2str(pd1.mu) '+-' num2str(pd1.sigma) 's'], 'FontSize',14) 
xlabel('stall duration (s)')
ylabel('Counts')
ax = gca;
ax.FontSize = 13;
title(stall_title);


% plot accumulative figure of resumption
resumption_recovery_duration_sort = sort(resumption_recovery_bin(resumption_recovery_bin<10000)); % only retain the resumed traces
recovery_fraction_temp = (1:length(stall_duration_bin))/length(stall_duration_bin);  % fraction that resumed
recovery_fraction = recovery_fraction_temp(1:length(resumption_recovery_duration_sort));
recovery_duration = [resumption_recovery_duration_sort,300];
duration_ind = find(resumption_recovery_duration_sort<300);
recovery_fraction_300s_ind = duration_ind(end);
recovery_fraction_300s_end = [recovery_fraction, recovery_fraction(recovery_fraction_300s_ind)];

subplot (1,2,2)
plot(recovery_duration, recovery_fraction_300s_end, '-o')
xlabel('resumption recovery time (s)', 'FontSize',14);
ylabel('Fraction recovered',  'FontSize',14);
title('cumulative plot of resumption',  'FontSize',14);
xlim ([0  300]); 
ylim ([0  1]); 


recovery_duration_DNAP_100nM_Dct_helicase_240s_stall = recovery_duration';
recovery_fraction_300s_end_DNAP_100nMDct_thelicase_240s_stall = recovery_fraction_300s_end';
save("resumption_recovery_100nM_DNAP_Dct_helicase_240s_stall.mat","recovery_duration_DNAP_100nM_Dct_helicase_240s_stall","recovery_fraction_300s_end_DNAP_100nMDct_thelicase_240s_stall")

