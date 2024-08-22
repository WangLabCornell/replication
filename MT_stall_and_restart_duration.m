%% load trace data to be analyzed
data_path = "/data/";
filename = "Trace_099_240131_171528.txt";
trace_name = extractBefore(filename, ".txt");
output_path = "/results/" + trace_name + "/";
mkdir(output_path);
output_filename = output_path + trace_name;


Trace_table = readtable(data_path + filename);
Time = Trace_table.Time_s_; % load Time
Magnet_turns = Trace_table.Magnet_turns; % load Magnet turns
Extension = Trace_table.Extension_nm_; % load Extension

% get Tasklist index: 1 is wait for stall, 2 is unwinding -100 turns, 3 is checking restart
Task1 = find(Magnet_turns == -30); % task 1 is wait for stall
Task2 = find(Magnet_turns < -30 & Magnet_turns > -130); % task 2 is unwind -100 turns
Task3 = find(Magnet_turns == -130); % task 3 is wait for restart
Timespot_index = [Task1(1); Task2(1); Task3(1);];

%% calculate the stall duration, judge restart activity, and calculate restart time
% find stall position and stall duration, resumption recovery duration (defined as reaching half extension expected)
dt = Time(2)- Time(1); % time step of data points
idx_1 = find(Extension(Timespot_index(1):Timespot_index(2))< 0.5*Extension(1)); % index for stalling, when extension <1/2
idx_2 = find(Extension(Timespot_index(3):length(Time))< 0.4*Extension(1)); % index for restart, when extension <1/2 of the remaining DNA template length
stall_duration_start = Time(idx_1(1)); % get the indext of stall start point   
stall_duration = Time(Timespot_index(2)) - stall_duration_start; % stall duration is counted from stall start to the unwinding step
% judge restart, and calculate the restart time 
if length(idx_2) == 0 % if extension never drops down to half extension, count as no resumption
   restart_time = [];
else
   restart_time = (idx_2(1)-1)*dt; % restart time if replication resumes
end 

%% plot and print analysis result
fig = figure(1);

subplot (3,1,1)   
plot (Time, Magnet_turns, 'Color', 'k', 'Linewidth', 2);
xlabel('Time (s)');
ylabel('Magnet Turns');
xlim([0 180])

subplot (3,1,[2,3])     
plot (Time, Extension, 'Color', 'k', 'Linewidth', 2);
line([stall_duration_start stall_duration_start], [0 1600], 'Color','red','LineStyle','--', 'Linewidth', 1); % mark the stall start point
line([Time(Timespot_index(2))  Time(Timespot_index(2))], [0 1600], 'Color','red','LineStyle','--', 'Linewidth', 1); % mark the stall start point
text(30, 1800, ['stall duration:' num2str(stall_duration) 's'], 'Color', 'r', 'FontSize',12); % print stall duration on figure
if length(idx_2) == 0 % if extension never drops down to half extension, count as no resumption
   text(100, 1500, ['no restart'], 'FontSize',12)
else
   line([Time(Timespot_index(3)) Time(Timespot_index(3))], [0 1600], 'Color','blue','LineStyle','--', 'Linewidth', 1);
   line([restart_time + Time(Timespot_index(3))  restart_time + Time(Timespot_index(3))], [0 1600], 'Color','blue','LineStyle','--', 'Linewidth', 1);
   text(100, 1500, ['restart time:' num2str(restart_time) 's'], 'Color', 'b', 'FontSize',12) 
end 
xlabel('Time (s)');
ylabel('Extension (nm)');    
xlim([0 180])

% Save plot and data

trace_data.stall_duration = stall_duration;
trace_data.restart_time = restart_time;
save(output_filename + ".mat","trace_data" )

saveas(fig,output_filename + ".png")


