% Script to measure the velocity of the replisome.


traces = dir("/data/Trace_230821_0186_Shuming/*_data.txt");

for k = 1:length(traces)

data_folder = string(traces(k).folder) + "/";
filename = string(traces(k).name);

trace_name = regexp(filename,".*(?=_data.txt)","match");
tracedata = readtable(data_folder+filename);
load(char(data_folder + trace_name + "_trace_data.mat"));
mkdir("/results/Trace_230821_0186_Shuming/")
trace_data.output_file_name = "/results/Trace_230821_0186_Shuming/" + trace_data.output_file_name;

time = tracedata.time;
x_raw = tracedata.time;
y_raw = tracedata.replicated_bp;
torque_raw = tracedata.torque;
dt = x_raw(2)-x_raw(1);


%Trim down the data
%Detect fork regression as a maximum followed by regression of at least 100
%bp.  If found, cut the trace at the maximum.
[y_max,I_max] = max(y_raw);
I_end = length(y_raw);
[y_min,~] = min(y_raw(I_max:end));
if(y_max-y_min > 100)
    I_end = I_max;
end

I_start = find(torque_raw<13,1,"first");
if(I_end <= I_start)
    I_end = I_start+1;
end

y = y_raw(I_start:I_end);
x = x_raw(I_start:I_end);
torque = torque_raw(I_start:I_end);

% Parameters
dwell_thresh = 0.1;
time_cutoff = x(end);
time_cutoff2 = x(end);


[y_smoothed,v] = calculate_velocity_position(x,y,0.5);


% Make figure and define subplots
fig = figure(1);
fig.Position = [10 10 1800 1500];
clf
rows = 3;
cols = 4;
ax_position = subplot(rows,cols,[1 2]);
ax_velocity = subplot(rows,cols,[9 10]);
ax_torque = subplot(rows,cols,[5 6]);
ax_dwell = subplot(rows,cols,3);
ax_velocity_binned = subplot(rows,cols,4);
ax_velocity_hist = subplot(rows,cols,[7,8]);
ax_velocity_regions = subplot(rows,cols,11:12);
linkaxes([ax_position ax_dwell ax_velocity_binned],'y');
linkaxes([ax_position ax_velocity, ax_torque],'x');

% Plot position vs time
subplot(ax_position);
hold on
%plot(x,y, 'LineWidth',2)
plot(x_raw,y_raw, 'LineWidth',2)
%plot(x,y_smoothed, 'LineWidth',2)
ylabel('Fork Position (bp)')
xlabel('Time (s)')
title(trace_data.output_file_name, "Interpreter","none");

% Plot Velocity vs time
subplot(ax_velocity);
plot(x,v, 'LineWidth',2);
ylabel('Velocity (bp/s)');
xlabel('Time (s)');

title('Velicity vs Time');

% Plot torque vs time
subplot(ax_torque);

plot(x_raw, torque_raw, 'LineWidth',2);
ylim([0,22]);
ylabel("Indirect Torque (pN nm)");
xlabel("Time (s)");
title("Torque vs Time");

% Dwell time histogram
subplot(ax_dwell);
bins = 700:3000;

h = histogram(y_smoothed,bins,'Orientation', 'horizontal');
h.BinCounts = h.BinCounts*dt;
xline(dwell_thresh,'--');
xlabel('Dwell Time (s/bp)')
title('Dwell Histogram');

%  Plot the velocity binned by the postion data.
ind = round(y_smoothed);
binned = accumarray(ind,v)./accumarray(ind,ones(size(v)));
subplot(ax_velocity_binned);
plot(binned,1:length(binned), 'LineWidth',2);
xline(0,'--');
xlabel('Velocity (bp/s)')
title('Velicity binned by position');

% Plot the histogram of velocity binned by position
subplot(ax_velocity_hist);
histogram(binned);
title('Binned velocity histogram')
xlabel('Velocity (bp/s)');
ylabel('Count')
average_velocity_binned = mean(binned(binned>10));
legend(sprintf("Mean (v>10) = %0.1f",average_velocity_binned))

pauses = get_pauses(h.BinCounts,bins,dwell_thresh);


subplot(ax_position);
hold on

% Create and index array to select all the regions with pausing.
index = false(size(y_smoothed));

for i = 1:length(pauses)
    min_index = find((y_smoothed > pauses(i).start),1);
    max_index = find((y_smoothed < pauses(i).end),1,'last');
    index = index | ((x > x(min_index)) & (x < x(max_index)));
end
plot(x(index),y_smoothed(index),'.','MarkerSize',10,'Color','Red')
%tracedata.paused = index;
xline(time_cutoff);


% Find pauses less than cutoff:

index_cutoff = find(x < time_cutoff,1,'last');
[pause_regions, numRegions] = bwlabel(index(1:index_cutoff));
pause_durations = zeros(numRegions,1);
pause_start_index = zeros(numRegions,1);
for i = 1:numRegions
    x_region = x(pause_regions == i);
    pause_durations(i) = x_region(end)-x_region(1);
    pause_start_index(i) = find(pause_regions == i,1,'first');
end

pause_duration = mean(pause_durations);

index_cutoff = find(x < time_cutoff2,1,'last');
[active_regions, numRegions] = bwlabel(~index(1:index_cutoff));
active_distances = zeros(numRegions,1);

for i = 1:numRegions
    
    y_region = y_smoothed(active_regions == i);
    active_distances(i) = y_region(end)-y_region(1);

end

active_distance = mean(active_distances);


%Find the procesivity
processivity = max(y)-y(1);
trace_data.velocity.processivity = processivity;

% Find the overall velocity from a linear fit to the position vs time.

durations = [30,60];

subplot(ax_position)
hold on


overall_velocities = zeros(size(durations));
overall_torques = zeros(size(durations));
for i = 1:length(durations)
    time_end = x(1)+durations(i);
    i_end = find(x>time_end,1,"first");
    if isempty(i_end) 
        i_end = length(x);
    end
    overall_x = x(1:i_end);
    overall_y = y(1:i_end);
    fit = polyfit(overall_x,overall_y,1);
    fit_y = polyval(fit,overall_x);
    plot(overall_x,fit_y, '--', 'LineWidth',2, 'Color', 'Black');
    overall_velocities(i) = fit(1);
    overall_torques(i) = mean(torque(1:i_end));
end

trace_data.velocity.durations = durations;
trace_data.velocity.overall_velocities = overall_velocities;
trace_data.velocity.overall_torques = overall_torques;


plot_text = [...
    sprintf('Processivity = %0.1f bp', processivity)...
    newline...
    sprintf('Overall Velocity (%d s) = %0.1f bp/s', durations(1), overall_velocities(1))...
    newline...
    sprintf('Overall Velocity (%d s) = %0.1f bp/s', durations(2), overall_velocities(2))...
    ];
text(10,2500,plot_text);



% Identify the active regions (not paused).
[active_regions, numRegions] = bwlabel(~index);

tstart = x_raw(I_start);

% Add the vertical lines to all the plots.
subplot(ax_position)
yline(y(1));
yline(y(1)+processivity);
xline(tstart)
for i = 1:length(durations)
    xline(x(1)+durations(i),'--');
end
subplot(ax_torque)
xline(tstart)
for i = 1:length(durations)
    xline(x(1)+durations(i),'--');
end
subplot(ax_velocity)
xline(tstart)
for i = 1:length(durations)
    xline(x(1)+durations(i),'--');
end


subplot(ax_torque)
hold on
average_torque = overall_torques(1);
plot(overall_x,ones(size(overall_x))*average_torque, '--', 'LineWidth',2, 'Color', 'Black')

plot_text = [sprintf('Average torque (%d s) = %0.1f pNnm',durations(1),overall_torques(1))...
    newline...
    sprintf('Average torque (%d s) = %0.1f pNnm',durations(2),overall_torques(2))...
    newline...
    sprintf('Average torque (%d s) = %0.1f pNnm',durations(2),overall_torques(2))...
    newline...
    sprintf('Average torque (%d s) = %0.1f pNnm',durations(2),overall_torques(2))...
];
text(10,5,plot_text);
xlim([0 200])

subplot(ax_velocity_regions)
hold on

velocity_regions = [];
for i = 1:numRegions
    region = active_regions == i;
    x_cropped = x(region);
    y_cropped = y_smoothed(region);

    fit = polyfit(x_cropped,y_cropped,1);
    velocity_regions(i).distance = polyval(fit,x_cropped(end)) - polyval(fit,x_cropped(1));
    velocity_regions(i).velocity = fit(1);
    velocity_regions(i).torque = mean(torque(region));

    
end
hold off
plot([velocity_regions.distance],[velocity_regions.velocity],'.','MarkerSize',30);
average_velocity = sum([velocity_regions.distance].*[velocity_regions.velocity])/sum([velocity_regions.distance]);
average_torque = sum([velocity_regions.distance].*[velocity_regions.torque])/sum([velocity_regions.distance]);
title('Velocity of active regions')
ylabel('Velocity of region (bp/s)');
xlabel('Distance of active region (bp)')
legend(sprintf("Weighted Mean = %0.1f bp/s" + newline + "Weighted torque = %0.1f pNnm",average_velocity, average_torque),'Location','Northwest');

trace_data.velocity.pause_free_velocity = average_velocity;
trace_data.velocity.pause_free_torque = average_torque;



set(findall(gcf,'-property','FontSize'),'FontSize',18)
saveas(fig, strcat(trace_data.output_file_name,"_velocity.png") );

save(trace_data.output_file_name+"_trace_data.mat","trace_data");



end


function pauses = get_pauses(dwell_hist,bins,dwell_thresh)
    
    ind = dwell_hist > dwell_thresh;
    
    % Label regions with unique value
    [labeledVector, numRegions] = bwlabel(ind);
    
    pauses = [];
    for i = 1:numRegions
        bins_selected = bins(labeledVector == i);
        pauses(i).start = min(bins_selected);
        pauses(i).end = max(bins_selected);
    end

end