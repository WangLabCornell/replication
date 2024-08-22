% Script to plot replication traces from the AOT
% Calculates the replicated basepairs from the extension.
% Export data for further analysis

clear;
close all;

% Make a list of all the traces in the current directory to analyze

traces = dir("/data/Trace_220805_0080_*");

%% Define which tasks to use
% Define some named taskes to use later in the code:
% These define the task list index for specific tasks.
% All the task list

task_plot_list = 2:5; % List of tasks to plot versus time.
task_torque_list = [3,5]; %Tasks to plot torque
task_HAT_const_height = 1; % Task of the initial constant trap height HAT curve.
task_HAT_const_height_post = 4; % Task to attempt to extract a post-replication contant trap height HAT curve.
task_initial_activity = 3; % Task to look for the replisome intitial activity where the DNA changes chirality from (-) to (+) supercoiled.
task_stalling = 3; % Task to determine stall start and stall torque value
task_velocity = 3; % Task to measure velocity

% Parameters for defining stalling
stall_start_thresh = 13; % stall start threshold, unit: pNnm
stall_window = 6000; % stall torque duration window in datapoints; 6000 is 60s

% Plotting options
linewidth = 1;
movemean_step = 50; % size of the moving window (for showing data)
movemean_torque = 2000; % size of the moving window (for showing torque data)
ylimits_FT_table = [-10 60];
surface_threshold = 300; %threshold for surface drift selection, unit: nm
ylimim_torque = [-20 40]; % range of torque plot

%% For loop
for i = 1:length(traces)
    
    % Import data from all tasks in a processed trace:
    % The resulting variable 'tasks' is a cell arrary of structs
    % tasks{i}.name File name of the task
    % tasks{i}.data Table of the data
    trace_data.name = traces(i).name;
    trace_data.folder = traces(i).folder;
    trace_data.index = i;
    trace_data.notes = "";
    tasks = import_trace(strcat(trace_data.folder,'/',trace_data.name));
    

    %Parse file and folder names to extract date, trace number, and user.
    %For this data set user name and date are in the name
    folder_parts = strsplit(trace_data.folder,"/");
    folder_parts = strsplit(folder_parts{end},"_");
    trace_name_parts = strsplit(trace_data.name,"_");
    trace_data.user = trace_name_parts{4};
    trace_data.date = trace_name_parts{2};
    trace_data.number = trace_name_parts{3}; 
    %Create a unique output filename for this trace.
    output_name = sprintf("%03d_%s_%s_%s",trace_data.index,trace_data.date,trace_data.number,trace_data.user);
    save_directory = "/results/Trace_" + output_name
    mkdir(save_directory);
    trace_data.output_file_name = save_directory + "/" + output_name ;

    time0 = tasks{task_plot_list(1)}.data.Time_s_(1);
   

    % Get the hat curves at constant trap height
    % Fit the initial HAT curve.  Fit with order 1 for quadratic, order 2
    % for quartic
    trace_data.HAT_const_height_order = 1;
    trace_data.HAT_const_height_fit = polyfit(tasks{task_HAT_const_height}.data.Turns.^2, tasks{task_HAT_const_height}.data.Extension_nm_, trace_data.HAT_const_height_order); 
    % Check if a post replication HAT curve exists and fit it.
    trace_data.HAT_post = fit_post_HAT_curve(tasks{task_HAT_const_height_post}.data.Turns,tasks{task_HAT_const_height_post}.data.Extension_nm_);


    % Get the locations of the replication sign flips:
    % Make sure we are at least 10 turns in.
    index_start = find(tasks{task_initial_activity}.data.Turns <= -10, 1, "first");
    [~,index_end] = max(tasks{task_initial_activity}.data.Force_pN_);
    
  
    [M, I] = max(movmean(tasks{task_initial_activity}.data.Extension_nm_(index_start:index_end),100));
    trace_data.fork_position.flip_times = tasks{task_initial_activity}.data.Time_s_(I+index_start -1);
    trace_data.fork_position.extension_max = M;
    

    if(trace_data.HAT_post.exists)
        [M, I] = max(movmean(tasks{task_HAT_const_height_post}.data.Extension_nm_,100));
        trace_data.fork_position.flip_times = [trace_data.fork_position.flip_times, tasks{task_HAT_const_height_post}.data.Time_s_(I)];
    end

    %Calculate the replicated bp
    for task = task_plot_list
        tasks{task}.data.replicated_bp = get_fork_position(tasks{task}.data.Time_s_,tasks{task}.data.Turns,tasks{task}.data.Extension_nm_ ,trace_data);
    end


    %% Start the main plot
    %Get the x limits for the time figures.
    
    xlimits = [0, tasks{task_plot_list(end)}.data.Time_s_(end) - time0];   
    
    f_main = figure(1);
    clf("reset");
    f_main.Position = [10 10 1800 1500];

    % Define subplots and link x-axis:
    ax_turns = subplot(6,4, [1,2,3]);
    ax_extension = subplot (6,4, [5,6,7]);
    ax_force = subplot (6,4, [9,10,11]);
    ax_replicated_bp = subplot (6,4, [13,14,15]);
    ax_torque = subplot (6,4, [17,18,19,21,22,23]);
    
    
   
    ax_hat_cons_height = subplot (6,4,8);

    
    %linkaxes([],'y')
    linkaxes([ax_turns, ax_extension, ax_force, ax_replicated_bp, ax_torque],'x');

    %% plot turns
    subplot(ax_turns);
    title(trace_data.output_file_name, 'Interpreter', 'none','FontSize', 18);
    hold on
    for task = task_plot_list
        plot (tasks{task}.data.Time_s_ - time0, tasks{task}.data.Turns, 'Linewidth', linewidth);
    end
    plot_xlines(trace_data.fork_position.flip_times - time0);
    ylabel('Turns');
    xlim(xlimits);

    %% plot extension
    subplot(ax_extension);
    hold on
    for task = task_plot_list
        plot (tasks{task}.data.Time_s_ - time0, movmean(tasks{task}.data.Extension_nm_,movemean_step), 'Linewidth', linewidth);
    end
    plot_xlines(trace_data.fork_position.flip_times - time0);
    ylabel('Extension (nm)');
    xlim(xlimits);
    
    %% plot force
    subplot(ax_force);
    hold on
    for task = task_plot_list
        plot (tasks{task}.data.Time_s_ -time0, movmean(tasks{task}.data.Force_pN_,movemean_step), 'Linewidth', linewidth);
    end
	plot_xlines(trace_data.fork_position.flip_times - time0);
    ylabel('Force (pN)');
    xlim(xlimits);

    %% plot replicated bp
    subplot(ax_replicated_bp);
    hold on
    for task = task_plot_list
        plot (tasks{task}.data.Time_s_ - time0, tasks{task}.data.replicated_bp, 'Linewidth', linewidth);
    end
	plot_xlines(trace_data.fork_position.flip_times - time0);
    xlim(xlimits);  
    ylim([-200 1500]);
    ylabel('Fork Position (bp)');
    
    %% plot torque 
    subplot (ax_torque)
    hold on
    for j = 1:length(task_plot_list)
        task = task_plot_list(j);
        if(find(task_torque_list == task))
            %Only plot the torque data for tasks selected to plot
            [time_trimmed,toque_trimmed] = trim_torque(tasks{task}.data.Time_s_, tasks{task}.data.Turns, tasks{task}.data.Torque_pNnm_); 
            plot (time_trimmed - time0, movmean(toque_trimmed, movemean_torque), 'Linewidth', 3, 'Color', default_MATLAB_colors(j));           
        end
        
    end
	plot_xlines(trace_data.fork_position.flip_times - time0);
    
    % Measure the stall torque
    [stall_start, indir_stall_torque_XXs, dir_stall_torque_XXs] = stall_start_stall_torque(tasks, task_stalling, movemean_step, movemean_torque, stall_start_thresh, stall_window); %find the stall start position
    if isempty(stall_start)        
    else
    stall_lines_plot(stall_start, time0, stall_window, indir_stall_torque_XXs, dir_stall_torque_XXs)
    end

    grid on
    xlim(xlimits);
    ylim(ylimim_torque);
    xlabel('Time (s)');
    ylabel('Torque (pNnm)');
    yline(-10, '--', 'Color', 'blue', 'Linewidth', 2);
    
    
    %% plot pre- and post-replication hat curves at constant height
    subplot (ax_hat_cons_height)
    hold on
    title('hat curve at constant trap height');
    plot(tasks{task_HAT_const_height}.data.Turns, movmean(tasks{task_HAT_const_height}.data.Extension_nm_,30)); % this is the pre-replication hat curve at constance height

    plot(tasks{task_HAT_const_height}.data.Turns,polyval(trace_data.HAT_const_height_fit,tasks{task_HAT_const_height}.data.Turns.^2),'Color',	lighten_plot_color(default_MATLAB_colors(1),0.5),'LineWidth',2);
    
    xlim([-10 70]);
    ylim([600 950]);
    ylabel('Extension (nm)');
    xlabel('Turns');
    
    saveas(f_main, strcat(trace_data.output_file_name,".png") );
    save(strcat(trace_data.output_file_name,"_trace_data.mat"),"trace_data");

    %% Special plots for further ananlysis
    % Uncomment to do the choosen analysis and export data.
    %resumption_plot(tasks,trace_data,task_plot_list);
    stall_plot(tasks,trace_data,task_plot_list);
    %velocity_plot(tasks,trace_data,task_velocity);
end

%% Functions
% ####################################################################
% ####################################################################
% ####################################################################
% ####################################################################
% ####################################################################
function velocity_plot(tasks,trace_data,task_velocity)
    movemean_step = 50;
    task = task_velocity;
    time = tasks{task}.data.times;
    time = time - time(1);
    replicated_bp = tasks{task}.data.replicated_bp;
    force = movmean(tasks{task}.data.FzpN, movemean_step);
    torqueFromForce = movmean(tasks{task}.data.TorqueFromForce_pNnm_, movemean_step);

    figure(2);
    clf;
    subplot(4,1,1:3);
    plot(time, replicated_bp);
    ylim([-100, 2500]);
    ylabel("Fork position (bp)");
    title(trace_data.output_file_name, 'Interpreter', 'none');
    subplot(4,1,4);
    plot(time, torqueFromForce);

    ylim([0,25]);
    ylabel("Torque" + newline + "(pNnm)");
    xlabel("Time (s)");

    T = table(time, replicated_bp, force, torqueFromForce);
    writetable(T,trace_data.output_file_name+ "_data.txt");
end

function stall_plot(tasks,trace_data,task_plot_list)
        %Collect data
    time = [];
    replicated_bp = [];
    torque = [];
    for task = task_plot_list(1:2)
        time = [time; tasks{task}.data.Time_s_ - trace_data.fork_position.flip_times(1)];
        replicated_bp = [replicated_bp; tasks{task}.data.replicated_bp];
        torque = [torque; movmean(tasks{task}.data.TorqueFromForce_pNnm_, 50)];
    end

    % Set t = 0 to when the torque rises.
    index = find(torque>10,1,'first');
    
    time = time - time(index);

    [replicated_bp_max, index_max] = max(replicated_bp);

    trace_data.regression.replicated_bp_max = replicated_bp_max;
    trace_data.regression.index_max = index_max;
    trace_data.regression.time_max = time(index_max);
    trace_data.regression.has_regression_data = false;

    if time(end) - trace_data.regression.time_max > 60
        trace_data.regression.has_regression_data = true;
    end
     
    fig = figure(2);
    fig.Position = [1510 10 1600 1400];
    linewidth = 2;
    ax_replicated_bp_stall = subplot(2,1,1);
    ax_torque_stall = subplot(2,1,2);
    linkaxes([ax_replicated_bp_stall,ax_torque_stall],'x');


    subplot(ax_replicated_bp_stall);
    %hold on
    plot(time, replicated_bp, 'Linewidth', linewidth);
    xline(0,'--');
    xline(trace_data.regression.time_max,'--',"LineWidth", linewidth);
    xline(trace_data.regression.time_max+60,'--',"LineWidth", linewidth);
    ylim([-100,2500]);
    xlim([-50, 200]);
    xlabel('Time (s)');
    ylabel('Fork Position (bp)');

    subplot(ax_torque_stall);
    %hold on
    plot(time, torque, 'Linewidth', linewidth);
    xline(0,'--');
    xline(trace_data.regression.time_max,'--',"LineWidth", linewidth);
    xline(trace_data.regression.time_max+60,'--',"LineWidth", linewidth);
    ylim([-15,30]);
    xlim([-50, 400]);
    
    
    if(trace_data.regression.has_regression_data)
        text(0,-5,"Included","Color","#77AC30")
    else
        text(0,-5,"Discarded","Color","Red")
    end
    xlabel('Time (s)');
    ylabel('Torque (pNnm)');
    set(findall(gcf,'-property','FontSize'),'FontSize',18)
    %Save output
    saveas(fig, strcat(trace_data.output_file_name,"_stall.png"));
    save(strcat(trace_data.output_file_name,"_trace_data.mat"),"trace_data");

    T = table(time, replicated_bp, torque);
    writetable(T,strcat(trace_data.output_file_name,"_stall.txt"));
end

function resumption_plot_combined(tasks,trace_data,task_plot_list)
    %Collect data
    time = [];
    turns = [];
    replicated_bp = [];
    torque = [];
    task_start_times = [];
    for task = task_plot_list(1:5)
        time = [time; tasks{task}.data.times];
        turns = [turns; tasks{task}.data.Turn];
        task_start_times = [task_start_times; tasks{task}.data.times(1)];
        replicated_bp = [replicated_bp; tasks{task}.data.replicated_bp];
        torque = [torque; movmean(tasks{task}.data.TorqueFromForce_pNnm_, 50)];
    end

    time0 = tasks{task_plot_list(3)}.data.times(1);
    time = time - time0;
    task_start_times = task_start_times - time0
    %Start plot
     
    fig = figure(2);
    linewidth = 2;
    ax_replicated_bp_stall = subplot(2,1,1);
    ax_torque_stall = subplot(2,1,2);
  
    linkaxes([ax_replicated_bp_stall,ax_torque_stall],'x');

    %First subplot
    subplot(ax_replicated_bp_stall);
    
    plot(time, replicated_bp, 'Linewidth', linewidth);
    
    plot_xlines(task_start_times)
    %xlim([-50, 200]);
    ylim([-250,2000]);
    xlabel('Time (s)');
    ylabel('Fork Position (bp)');
    title(strcat(trace_data.output_file_name," (",trace_data.description,")"), 'Interpreter', 'none','FontSize', 18);

    % Second subplot
    subplot(ax_torque_stall);
    
    plot(time, torque, 'Linewidth', linewidth);
    plot_xlines(task_start_times)
    xlim([-50, 200]);
    ylim([0,25]);
    xlabel('Time (s)');
    ylabel('Indirect Torque (pNnm)');
    torque_smoothed = movmean(tasks{task_plot_list(2)}.data.torque,100);
    text(10,7,sprintf("Actual torque at start of stall : %0.1f pN", torque_smoothed(end) ),'FontSize', 12);
    text(10,3,sprintf("Stall duration : %0.0f s", task_start_times(4)- task_start_times(3)),'FontSize', 12);


    %Save output
    saveas(fig, strcat(trace_data.output_file_name,"_resumption.png"));
    save(strcat(trace_data.output_file_name,"_trace_data.mat"),"trace_data");

    T = table(time, turns, replicated_bp, torque);
    writetable(T,strcat(trace_data.output_file_name,"_resumption.txt"));

end

function data = importfile(filename)
    warning off
    data = readtable(filename, 'Delimiter','tab', 'VariableNamingRule','modify');
    warning on
end

function tasks = import_trace(folder)
    files = dir(folder);
    files = files(3:end);
    tasks = cell(size(files));
    for j = 1:length(files)
        tasks{j} = struct('name',files(j).name,'data',importfile(strcat(files(j).folder,'/',files(j).name)));
    end
end

function [time_trimmed, torque_trimmed] = trim_torque(time, turns, torque)
    % Function to trim the beginning of the torque data
    bool =(mod(turns,1) == 0);
    idx = find(bool);
    if ~isempty(idx)
        starting_index = idx(1);
    else
        starting_index = 1;
    end
    torque_trimmed = torque(starting_index:end);
    time_trimmed = time(starting_index:end);
end

function fork_position = get_fork_position(time,turns,extension,trace_data)
    extension_smoothed = movmean(extension,100);
    p = trace_data.HAT_const_height_fit;
    %Replace the offset term in the const_height HAT curve with new offset.
    p(end) = trace_data.fork_position.extension_max;
    replication_turns_raw = get_turns_from_extension(p,extension_smoothed);
    sign = get_sign(time,trace_data.fork_position.flip_times);
    fork_position = (sign.*replication_turns_raw - turns)*10.5;
    fork_position(imag(fork_position) ~= 0 ) = NaN;
    fork_position = real(fork_position);
end

function sign = get_sign(a,t)
    
    index = false(size(a));
    for i  = 1:length(t)
        if(mod(i,2) == 0)
            index = index & (a<t(i));
        else
            index = index | (a>=t(i));
        end
    end
    sign = (index).*2-1;
end

function post_HAT = fit_post_HAT_curve(turn,extension)
    post_HAT.turn = turn;
    post_HAT.extension = extension;
    segment_window = 300; % set a window for fitting the hat curve peak after unwinding
    extension_thresh = 20; % threshold for testing if the extension changed.
    
    extension_smoothed = movmean(extension,30);
    [ext_max, I_max] = max(extension_smoothed);
    [ext_min, ~] = min(extension_smoothed);
    post_HAT.exists = false; % default to no post HAT exists.
    post_HAT.fit = [];
    post_HAT.turn_shift = 0;
    post_HAT.max_extension = ext_max;
    if((ext_max - ext_min) > extension_thresh)
        % If the extension changed.
        if I_max-segment_window>0 && I_max+segment_window<length(extension_smoothed) % There is enough data to select around the peak
            extension_smoothed_select = extension_smoothed(I_max-segment_window:I_max+segment_window); % select a region centering at the hat curve peak
            turn_select = turn(I_max-segment_window:I_max+segment_window);
            post_HAT.fit = polyfit(turn_select, extension_smoothed_select, 2); % fit the 2nd unwinding using a poly function
            %post_HAT.extension_fit = polyval(post_HAT.fit,turn_select); % calculate the "fitted data"
            post_HAT.turn_shift = -post_HAT.fit(2)/(2*post_HAT.fit(1));
            post_HAT.max_extension = polyval(post_HAT.fit, post_HAT.turn_shift);
            post_HAT.turn_fit = turn_select - post_HAT.turn_shift;
            post_HAT.extension_fit = polyval(post_HAT.fit,turn_select);
            post_HAT.exists = true; %Post HAT exists.
        end
    end
end

function plot_color_light = lighten_plot_color(plot_color,percent)
    %Lighten plot color by mixing with white.
    plot_color_light = plot_color + percent*(1-plot_color);
end

function color = default_MATLAB_colors(index)
    if index < 1
        index = 1;
    end
    index = mod(index-1,7)+1;
    plot_colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
    color = plot_colors{index};
end

function x = quadratic(p,y)
    x = real(sqrt((y-p(2))/p(1)));
end

function turns = get_turns_from_extension(p,extension)
    if length(p) == 2
        turns = quadratic(p,extension);
    elseif length(p) == 3
        %model is extension = p1*turns^4 + p2*turns^2 + p3
        
        turns = sqrt((-p(2)-sqrt(p(2).^2 - 4*p(1)*(p(3)-extension)))/(2*p(1)));
    end
end

function plot_xlines(xlines)
    for i= 1:length(xlines)
        xline(xlines(i),'--');
    end
end

function [stall_start, indire_stall_torque_XXs, dir_stall_torque_XXs] = stall_start_stall_torque(tasks, task_stalling, movemean_step, movemean_torque, stall_start_thresh, stall_window)
    % Find the stall torque from the trace data
    torqueFromForce = movmean(tasks{task_stalling}.data.TorqueFromForce_pNnm_, movemean_step); 
    torque   = movmean(tasks{task_stalling}.data.Torque_pNnm_, movemean_torque);
    Time_stall = tasks{task_stalling}.data.Time_s_;
    idx = find(torqueFromForce > stall_start_thresh);
    if isempty(idx)
       stall_start = []; 
       indire_stall_torque_XXs =[];
       dir_stall_torque_XXs =[];
    else    
    stall_start = Time_stall(idx(1));
         if length(torque(idx(1):end))>stall_window
            indire_stall_torque_XXs = max(torqueFromForce(idx(1):idx(1)+stall_window));
            dir_stall_torque_XXs = max(torque(idx(1):idx(1)+stall_window));
         else
            indire_stall_torque_XXs = max(torqueFromForce(idx(1):end));
            dir_stall_torque_XXs = max(torque(idx(1):end));
         end
    end
end

function stall_lines_plot(stall_start, time0, stall_window, indir_stall_torque_XXs, dir_stall_torque_XXs)
    xline(stall_start - time0, '--', 'Color', 'blue', 'Linewidth', 2);
    xline(stall_start - time0 + stall_window/100, '--', 'Color', 'blue', 'Linewidth', 2);
    
    text(stall_start - time0, -10, ['stall start to 60s'], 'FontSize', 14);
    text(stall_start - time0, 35, ['Stall torque within 60s:    ' num2str(dir_stall_torque_XXs)], 'FontSize', 14);
end