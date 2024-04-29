
close all;
clear all;

files = dir("*.txt");

stall_time = 30;

all_data = cell(length(files),1);

for k = 1:length(files)
    data = readtable(strcat(files(k).folder,"\",files(k).name));
    data.time = data.time - stall_time;
    data.replicated_bp = real(data.replicated_bp);

    resumption_plot_combined(data)
    all_data{k} = data;
end

%Create the average trace.
average_data = table;
for k = 1:width(data)
    data_x = cell(length(all_data));
    data_y = cell(length(all_data));
    for j = 1:length(all_data)
        cur_table = all_data{j};
        data_x{j} = cur_table{:,1};
        data_y{j} = cur_table{:,k};
    end
    average_data{:,k} = average(data_x,data_y)';
end
average_data.Properties.VariableNames = data.Properties.VariableNames;

%Plot the average trace in a new figure.
figure(1);
resumption_plot_combined(average_data,"Color","Black","Linewidth",4);
subplot(2,1,1)
times = [-stall_time, 0 , 12.48];
for t = times
    xline(t,'--',"LineWidth",1.5);
end
subplot(2,1,2)

for t = times
    xline(t,'--',"LineWidth",1.5);
end


figure(2);
resumption_plot_combined(average_data,"Color","Black","Linewidth",4);
subplot(2,1,1)

for t = times
    xline(t,'--',"LineWidth",1.5);
end
subplot(2,1,2)


for t = times
    xline(t,'--',"LineWidth",1.5);
end


writetable(average_data,"average_trace.csv")


function resumption_plot_combined(data,varargin)
    %Parse vargin to see if we specify plot color
    if length(varargin) == 0
        varargin = {"LineWidth", 2};
    end


    %Collect data
    time = data.time;
    turns = data.turns;
    replicated_bp = data.replicated_bp;
    torque = data.torque;

    %Start plot

    %fig = figure(2);
    
    ax_replicated_bp = subplot(2,1,1);
    ax_torque = subplot(2,1,2);

    linkaxes([ax_replicated_bp,ax_torque],'x');

    %First subplot
    subplot(ax_replicated_bp);
    hold on
    plot(time, replicated_bp, varargin{:});

    %xlim([-50, 200]);
    ylim([-250,2000]);
    xlabel('Time (s)');
    ylabel('Fork Position (bp)');

    % Second subplot
    subplot(ax_torque);
    hold on
    plot(time, torque, varargin{:});

    xlim([-150, 150]);
    ylim([0,25]);
    xlabel('Time (s)');
    ylabel('Indirect Torque (pNnm)');


end


function result = average(data_x,data_y)

    % Determine the required size of the matrix
    align_val = 0;
    rows = length(data_x);
    
    max_right = -Inf;
    max_left = -Inf;
    index = zeros(size(data_x));
    
    for k=1:rows
        cur_data = data_x{k};
        index(k) = find(cur_data>=align_val,1,"first");
        left = index(k);
        right = length(cur_data) - index(k);
        if left > max_left
            max_left = left;
        end
        if right > max_right
            max_right = right;
        end
    
    end
    
    cols = max_left + max_right ;
    all_data = nan([rows,cols]);
    
    for k = 1:rows
        cur_data = data_y{k};
        left_index = max_left-index(k) + 1;
        right_index = left_index + length(cur_data) -1;
        all_data(k,left_index:right_index) = cur_data;
    end
    
    result = mean(all_data,"omitnan");
end