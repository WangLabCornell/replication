% Script to load analyzed traces and create and average trace.
% Shift the time = 0 to the maximum fork position to create the average
% trace for the regression plot.

clear;
close all;

files = dir("*stall.txt");
files_trace_data = dir("*_trace_data.mat");

all_data = {};

for k = 1:length(files)
    load(files_trace_data(k).name)
    if(trace_data.regression.has_regression_data)
        all_data{end+1} = readtable(strcat(files(k).folder,"\",files(k).name));
        all_data{end}.time = all_data{end}.time - trace_data.regression.time_max;
        all_data{end}.replicated_bp_shifted = all_data{end}.replicated_bp - trace_data.regression.replicated_bp_max;
    end
end

length(all_data)

%Create the average trace.
average_data = average(all_data);

% Save the average traces
writetable(average_data,"average_trace.txt")

% Plotting
ax_replicated_bp = subplot(2,1,1);
ax_torqueFromForce = subplot(2,1,2);
linkaxes([ax_replicated_bp, ax_torqueFromForce],'x');

subplot(ax_replicated_bp);
hold on
for k = 1:length(all_data)
    plot(all_data{k}.time,all_data{k}.replicated_bp);
end
plot(average_data.time,average_data.replicated_bp,"LineWidth",3,"Color","Black");

SEM = average_data.replicated_bp_SD./sqrt(average_data.N);
plot(average_data.time,average_data.replicated_bp + SEM, "Color", [0.2 0.2 0.2])
plot(average_data.time,average_data.replicated_bp - SEM, "Color", [0.2 0.2 0.2])
ylabel("Fork position (bp)")


subplot(ax_torqueFromForce);
hold on
for k = 1:length(all_data)
    plot(all_data{k}.time,all_data{k}.torque);
end
plot(average_data.time,average_data.torque,"LineWidth",3,"Color","Black");
plot(average_data.time,average_data.torque_mean_plus_SEM, "Color", [0.2 0.2 0.2])
plot(average_data.time,average_data.torque_mean_minus_SEM, "Color", [0.2 0.2 0.2])
xlim([-300 300])
ylabel("Indirect torque (pN nm)")
xlabel("Time (s)")

set(findall(gcf,'-property','FontSize'),'FontSize',12)

saveas(gcf(), "average_trace.png" );




function output_table = average(data_tables)
    % Assume for now the first column in the data table is the one to use
    % as the value to average over.
    % This function does not require data_tables have the same length.
    % Requires data to be sampled at equal intervals
    
    % Determine the required size of the matrix

    % What x value should be used to align all the traces? 
    % This value must appear in all data sets.
    % Code assumes data is ascending.
    align_val = 0;
    align_col = 1;

    % Rows is the number of data sets.
    num_data_tables = length(data_tables);
    
    % Find the length of the output data table
    % Also find the 
    max_right = -Inf;
    max_left = -Inf;
    align_index = zeros(num_data_tables,1);
    
    for k=1:num_data_tables
        cur_table = data_tables{k};
        align_index(k) = find(cur_table{:,align_col}>=align_val,1,"first");
        left = align_index(k);
        right = length(cur_table{:,align_col}) - align_index(k);
        if left > max_left
            max_left = left;
        end
        if right > max_right
            max_right = right;
        end
    end
    
    VariableNames = string(cur_table.Properties.VariableNames);
    output_table_length = max_left + max_right;
    % Initialize output table to an empty table.
    output_table = table;

    % Make align column and 
    all_data = nan([num_data_tables,output_table_length]);
    for k = 1:num_data_tables
        cur_table = data_tables{k};
        left_index = max_left-align_index(k) + 1;
        right_index = left_index + length(cur_table{:,align_col}) -1;
        all_data(k,left_index:right_index) = cur_table{:,align_col};
    end

    
    output_table.(VariableNames(align_col)) = mean(all_data,"omitnan")';
    output_table.N = sum(~isnan(all_data))';

    cols = 1:length(VariableNames);
    for j = cols(cols~=align_col)
        all_data = nan([num_data_tables,output_table_length]);
        for k = 1:num_data_tables
            cur_table = data_tables{k};
            left_index = max_left-align_index(k) + 1;
            right_index = left_index + length(cur_table{:,j}) -1;
            all_data(k,left_index:right_index) = cur_table{:,j};
        end
        [col_std, col_mean] = std(all_data,"omitnan");
        output_table.(VariableNames(j)) = col_mean';
        SD_col_name = VariableNames(j) + "_SD";
        output_table.(SD_col_name) = col_std';
        
        output_table.(VariableNames(j) + "_SEM") = col_std'./sqrt(output_table.N) ;
        output_table.(VariableNames(j) + "_mean_plus_SEM") = output_table.(VariableNames(j)) + output_table.(VariableNames(j) + "_SEM");
        output_table.(VariableNames(j) + "_mean_minus_SEM") = output_table.(VariableNames(j)) - output_table.(VariableNames(j) + "_SEM");
    end


end
