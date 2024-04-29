clear;
close all;
clc;


tether_sel_file = 'tether_selection.txt';
key = '1 Wind';
%key = '5 unclear';

stall_threshold = 0.5;
resume_threshold = 0.4;
%%

stall_duration = [];
resumption_recovery_duration = [];
Count_inactive = 0;
Count_excluded = 0;


L0 = 0.266; % contuor length per bp at 0.5pN
pickedcolor = [255, 0, 0]/255;
EndingT = 11;
tweezer2 = 0;
tweezer2 = 2;
resuming_index = 1;

flag = 0; % To indicate if we choose any files at all, even if we choose one file, this would be 1;
file_index = 0; % Number of files chosen;
Maximum_allowed_files = 1;
Maximum_allowed_boxes = 200;
Data_points_number = zeros(Maximum_allowed_files,1);
position_and_force = zeros(4,Maximum_allowed_boxes,Maximum_allowed_files);

essential_columns = {{'DAQ time (s)'}; {'Time (ms)'}; {'tasklist (index)'};{'Magnet Angle (Turns)'};{'Magnet Height (mm)'}};% This is the number of columns essential for our hatcurve;


for ii = 1:Maximum_allowed_files
    position_and_force(2,:,ii) = (1:Maximum_allowed_boxes)-1;
end


while (1) % while not stopped, this loop keep on choosing files
    %% Access data files
    
    Address = 'Y:\Magnetic Tweezers\XJ67';
    if tweezer2 == 2
        Address = 'Y:\Magnetic Tweezers\XJ67';
    end
    
    [filename, pathname] = uigetfile('*.*','Select the data file',Address); % Choose MT data files
    if filename == 0 % If you do not choose a file, then the program proceeds;
        break;
    else flag = 1;
    end
    new_file_name = strrep(filename,'_',' ');
    fileID = fopen ([pathname,filename],'r'); % Open the data file
    file_index = file_index + 1;
    new_folder = [pwd,'\',[filename,' all']]; 
    mkdir([filename,' all']);
    
    new_folder_2 = [pwd,'\',[filename,' selected']]; 
    mkdir([filename,' selected']);
    
    %% Access force calibration files
    pathname_Force = pathname(1:(length(pathname)-length('Data Files\'))); % This is to access the force calibration results
    filename_Force = tether_sel_file; % Find the force calibration file
    fileID_Force = fopen([pathname_Force,filename_Force],'r'); % Open the force calibration file
    
    %% Get Title lines (the first line of each data file should be word descriptions, ditch that)
    titleline = fgetl(fileID);
    essential_columns_position = find_essential_columns_by_name( titleline, essential_columns );
    Box_cell_position  = find_boxes( titleline );
    num_of_box_cells = length(Box_cell_position);
    numofboxes = num_of_box_cells/8;
    end_of_array = Box_cell_position(end);
    Box_Ext_positions = (Box_cell_position(1)+6):8:end_of_array;
    % This is the actually num of boxes, not the maximum box index. If the
    % numofboxes = 10, the maximum box index would be 9.
    
    %% Figure out which box to plot
    % The boxes being plot should have the category of "1 Wind", which is
    % the content of 'key'.
    TCindex = (1:numofboxes) - 1; % All of the box indexes
    return_matrix = find_TC_index (fileID_Force, key);
    plotting_index = return_matrix(1,:); % These are the <box indexes> that are good traces
    Start_force = return_matrix(2,:); % This is the first calibrated force
    Jump_force = return_matrix(3,:); % This is the second calibrated force
    position_and_force(1,plotting_index+1,file_index) = 1;
    position_and_force(3,plotting_index+1,file_index) = Start_force;
    position_and_force(4,plotting_index+1,file_index) = Jump_force;
    
    %% Format Spec for reading data files
    good_box_index = plotting_index+1;
    good_box_ext_position = Box_Ext_positions(good_box_index);
    reading_position = [essential_columns_position, good_box_ext_position];
    format_spec = [];
    for ii = 1:end_of_array
        if isempty(find(reading_position==ii))
            format_spec = [format_spec, '%*f'];
        else
            format_spec = [format_spec, '%f'];
        end
    end
    Data_all_cell = textscan(fileID, format_spec, 'MultipleDelimsAsOne',true, 'Delimiter','[;', 'HeaderLines',0);
    %% Close all the files   
    fclose (fileID);
    fclose (fileID_Force);
    disp (['File ', num2str(file_index), ': ', filename]);
    
    %% A precaution not to exceed the maximum allowed number of files
    if file_index == Maximum_allowed_files
        break;
    end
end



if flag == 1 % This means we actually choose some files
    %% Find Time spots
    time_for_timespot = Data_all_cell{2}/1000; %This is only used for determining the steps; We use the first chosen file for that.
    task_for_timespot = Data_all_cell{3}; % The task number for the first file selected
    magnet_height_for_subplot = Data_all_cell{5}; % The magnetic height for the first file selected, used for plotting subplot;
    magnet_turns_for_subplot = Data_all_cell{4}; % The magnetic turns for the first file selected, used for plotting subplot;
    Timespot_index = timespot_finding (task_for_timespot);% Time spots given by the first file selected;
    Timespot = time_for_timespot(Timespot_index);
    EndingT = length(Timespot);
    
    total_offset = [];  % total offset due to replication
    shifted_turns = [];
%%
    for which_file = 1:file_index
        %Extension_index = (21+tweezer2):8:4000;
        plotting_index = find (position_and_force(1,:,which_file)==1);
        Time = Data_all_cell{2}/1000;

        for boxes = resuming_index:length(plotting_index)
            %for boxes = resuming_index:5

            figure (boxes)
            fig_handle = figure (boxes);
            pos = [200 200 1200 550];
            set(fig_handle, 'Pos', pos);
            %% Plotting subplots: Magnetic turns
            subplot (4,1,1)
            plot (time_for_timespot, magnet_turns_for_subplot,'k','LineWidth',1.5);
            axis([0 inf -inf inf])
            ylabel('M-turns');
            xlim ([Timespot(1) Timespot(end)]);
            ylim ([-200 51]);
            title ([new_file_name(1:(length(filename)-4)),' box ',num2str(plotting_index(boxes)-1)]);
            h = gca;
            h.XAxis.Visible = 'off';
            box(h,'off');
            set(h,'FontSize',13,'LineWidth',1.5);

            ymin = -100;
            ymax = 2200;

            Extension = 1000*Data_all_cell{boxes+length(essential_columns_position)};
            Extension_unfiltered = Extension;
            Extension = movmean(Extension, 5);
            dt = Time(2)- Time(1);

            Fstart = position_and_force(3,plotting_index(boxes),which_file);
            Fjump = position_and_force(4,plotting_index(boxes),which_file);
            Tracename = ['File',num2str(which_file),' Box',num2str(plotting_index(boxes)-1),'; Force: ',num2str(Fstart),' pN ~ ',num2str(Fjump),' pN.'];
            
            
            subplot (4,1,4)
            post_unwind_window = 600; % observe 15s pst replication
            plot (Time(Timespot_index(4):Timespot_index(5)), Extension(Timespot_index(4):Timespot_index(5)), 'Color', [0.4660 0.6740 0.1880], 'DisplayName',Tracename, 'Linewidth', 2);
            hold on
            plot (Time(Timespot_index(5): Timespot_index(5)+post_unwind_window), Extension(Timespot_index(5): Timespot_index(5)+post_unwind_window), 'Color', [0.9290 0.6940 0.1250], 'DisplayName',Tracename, 'Linewidth', 2);
            ylim([0 1800]);
            
            xlabel('Time (s)');
            ylabel('Extension (nm)');
            
            
            
            %% Plotting all the good traces
            subplot (4,1,[2,3])     
            plot (Time, Extension_unfiltered, 'Color', [0.7 0.7 0.7], 'Linewidth', 0.2);
            hold on;
            plot (Time, Extension,'Color', [0.9290 0.6940 0.1250], 'Linewidth', 1);
            xlim ([Timespot(1) Timespot(end)]);
            ylim ([ymin  ymax]); 
            xlabel('Time (s)');
            ylabel('Extension (nm)');
            
            for t = [2,4,6]
                plot (Time(Timespot_index(t):Timespot_index(t+1)), Extension(Timespot_index(t):Timespot_index(t+1)), 'Color', [0.4660 0.6740 0.1880], 'DisplayName',Tracename, 'Linewidth', 2);
                hold on
            end
            
             for t = 8
                plot (Time(Timespot_index(t):Timespot_index(t+1)), Extension(Timespot_index(t):Timespot_index(t+1)), 'Color', [0.8500 0.3250 0.0980], 'DisplayName',Tracename, 'Linewidth', 2);
                hold on
            end           
            

            for i = 2:9
                hold on
                line([Timespot(i) Timespot(i)], [ymin  ymax], 'Color','magenta','LineStyle','--', 'Linewidth', 0.5);   
            end 
            
            % extension offset due to replication, assuming 100 turns
            % replicated
            %Extension_offset = (0.266-0.031)*100*10.5;
            %Peak_hat_no_resumption = Extension(1) - Extension_offset;
            %line([Timespot(4)-20 Timespot(5)+20], [Peak_hat_no_resumption  Peak_hat_no_resumption], 'Color','blue','LineStyle','--', 'Linewidth', 0.5);
            %text(Timespot(5)+20, Peak_hat_no_resumption, ['expected hat peak if replicated 100 turns and no resumption during unwinding'], 'FontSize',10, 'Color','blue');
            
            text(Timespot(2)-30, 2300, ['add -30 turns'], 'FontSize',12, 'Color',[0.4660 0.6740 0.1880]);
            text(Timespot(4)-30, 2300, ['add -100 turns'], 'FontSize',12, 'Color',[0.4660 0.6740 0.1880]);
            text(750, 2300, ['add -80 turns then wind back to 0'], 'FontSize',12, 'Color',[0.8500 0.3250 0.0980]);
            %text(400, 1900, ['add -80 turns then wind back to 0'], 'FontSize',12, 'Color',[0.8500 0.3250 0.0980]);
            
            % find stall position and stall duration, resumption duration(defined as reaching half extension)
            if Extension(Timespot_index(1))<1150 % exclude short tethers
               stall_duration_temp = 100000;
               resumption_recovery_duration_temp = 100000;
               text(10, 2050, ['short tether'], 'FontSize',14) 
               Count_excluded = Count_excluded + 1;
            else
            
            idx_1 = find(Extension(Timespot_index(3):Timespot_index(4))< stall_threshold*Extension(1)); % index for stall extension <1/2
            idx_2 = find(Extension(Timespot_index(5):Timespot_index(6))< resume_threshold*Extension(1)); % index for resumptiom extension <1/2
            
              if length(idx_1) ==0 % exclude inactive traces
                  stall_duration_temp = 100000;
                  resumption_recovery_duration_temp = 100000;
                  text(10, 2050, ['inactive'], 'FontSize',14) 
                  Count_inactive = Count_inactive+1;
              else
                   stall_duration_temp_start = Time(idx_1(1) + Timespot_index(3)) ;            
                   stall_duration_temp = Timespot(4) - stall_duration_temp_start;
                   line([stall_duration_temp_start stall_duration_temp_start], [300 1500], 'Color','red','LineStyle','--', 'Linewidth', 0.2);
                   text(10, 2050, ['Extension:' num2str(Extension(1)) 'nm'], 'FontSize',14) 
                   text(10, 1850, ['stall duration:' num2str(stall_duration_temp) 's'], 'Color', 'r', 'FontSize',14) 
            
                    if length(idx_2) == 0 
                       resumption_recovery_duration_temp = 100000;
                       text(500, 1850, ['no resumption'], 'FontSize',14)
                       %text(200, 1850, ['no resumption'], 'FontSize',14)
                    else
                       resumption_recovery_duration_temp = (idx_2(1)-1)*dt; %duration takes to reach half extension again in resumption step
                       line([resumption_recovery_duration_temp + Time(Timespot_index(5))  resumption_recovery_duration_temp + Time(Timespot_index(5))], [300 1500], 'Color','blue','LineStyle','--', 'Linewidth', 0.2);
                       text(500, 1850, ['resumption recovery duration:' num2str(resumption_recovery_duration_temp) 's'], 'Color', 'b', 'FontSize',14) 
                       %text(200, 1850, ['resumption recovery duration:' num2str(resumption_recovery_duration_temp) 's'], 'Color', 'b', 'FontSize',14) 
                    end 
              end                   
                
            end               
            
            if stall_duration_temp == 100000
            else
              stall_duration = [stall_duration, stall_duration_temp]; 
              resumption_recovery_duration = [resumption_recovery_duration, resumption_recovery_duration_temp];
              saveas(gcf,[new_folder_2,['\ZBox ', num2str(plotting_index(boxes)-1),'.png']]);
            end          
   
            saveas(gcf,[new_folder,['\ZBox ', num2str(plotting_index(boxes)-1),'.png']]);
      
        end
        
        %output = [plotting_index'-1, mov, mov_rate, start_moving_time, breaking_time, start_mov_rec, finish_mov_rec, finish_flowing];
    end
    %% Plotting time devide lines
    
else % There is a possibility that one doesn't choose any files
    disp('You did not choose any files!');
end

Count_total = length(plotting_index); % total TC tethers
stall_duration_DNAP_100nM_Dct_helicase_45nM_240320_ch1 = stall_duration;
resum_recov_duration_DNAP_100nM_Dct_helicase_45nM_240320_ch1 = resumption_recovery_duration;
save('stall_duration_and_resumption_mark_DNAP_100nM_Dct_helicase_45nM_240320_ch1.mat','stall_duration_DNAP_100nM_Dct_helicase_45nM_240320_ch1', 'resum_recov_duration_DNAP_100nM_Dct_helicase_45nM_240320_ch1', 'Count_inactive', 'Count_excluded', 'Count_total');

figure
h1 = histfit(stall_duration,5)
h1(1).FaceColor = [.8 .8 1];
pd1 = fitdist(stall_duration','Normal')
text(10, 3.5, ['N=' num2str(length(stall_duration))], 'FontSize',14) 
text(10, 3, ['Mean:' num2str(pd1.mu) '+-' num2str(pd1.sigma) 's'], 'FontSize',14) 
xlabel('stall duration (s)')
ylabel('Counts')
ax = gca;
ax.FontSize = 13;
title('stall duration');


figure
h2 = histfit(resumption_recovery_duration,5)
h2(1).FaceColor = [.8 .8 1];
pd2 = fitdist(resumption_recovery_duration','Normal')
text(3, 5.5, ['N=' num2str(length(resumption_recovery_duration))], 'FontSize',14) 
text(3, 5, ['Mean:' num2str(pd2.mu) '+-' num2str(pd2.sigma) 's'], 'FontSize',14) 
xlabel('resumption recovery time (s)')
ylabel('Counts')
ax = gca;
ax.FontSize = 13;
title('resumption recovery duration');


% plot accumulative figure of resumption
resumption_recovery_duration_sort = sort(resumption_recovery_duration(resumption_recovery_duration<10000)); % only retain the resumed traces
recovery_fraction_temp = (1:length(stall_duration))/length(stall_duration);  % fraction that resumed
recovery_fraction = recovery_fraction_temp(1:length(resumption_recovery_duration_sort));
recovery_duration = [resumption_recovery_duration_sort,300];
recovery_fraction_300s_end = [recovery_fraction, recovery_fraction(end)];


figure
plot(recovery_duration, recovery_fraction_300s_end, '-o')
xlabel('resumption recovery time (s)', 'FontSize',14);
ylabel('Fraction recovered',  'FontSize',14);
title('cumulative plot of resumption fraction vs recovery duration',  'FontSize',14);
xlim ([0  300]); 
ylim ([0  1]); 


%% Functions
function return_matrix = find_TC_index (fileID_Force, key)
    % This function is to deal with the force calibration file. We want to
    % get which traces are good traces (in the previously catagorized good
    % trace catagory), and what are the two forces (there supposed to be
    % two forces);
    
    titleline = fgetl(fileID_Force); % Get the title line, since the first line is just text
    
    ii = 1;
    box_index = 0;
    while feof(fileID_Force)==0
        string = fgetl(fileID_Force);
        is_TC = strfind(string,key);
        if isempty(is_TC) == 1
        else
            TC_bead_index_matrix (ii) = box_index;
            starting = is_TC+length(key)+1;
            shorter_string = string(starting:end);
            shorter_data = str2num (shorter_string);
            Force (1,ii) = shorter_data(2);
            Force (2,ii) = shorter_data(1);
            ii = ii+1;
        end
        box_index = box_index+1;
    end
    return_position = TC_bead_index_matrix (1:ii-1);
    return_matrix = [return_position; Force(:,1:ii-1)];
end

function index = timespot_finding (task)
    % This is a function to find the time for changing to different task
 
    index = length(task);
    task_temp = max(task)-task;
    while (1)
        mask = sign(task_temp);
        if mask == 0
            break;
        end
        index = [sum(mask);index];
        task_temp = task_temp-1;
        task_temp = task_temp.*mask;
    end
    %timespot = time(index);
end

function [fitresult, gof] = Linear_Fit(xx, yy)

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( xx, yy );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );
end

function [ essential_columns_position ] = find_essential_columns_by_name( title_line, essential_columns )
    % This function returns the column numbers of the columns we are
    % caring, searching by name.
    tabs_pos = find(title_line == 9); % Find were are the "tab" signs are
    tabs_pos = [tabs_pos, (length(title_line)+1)];
    [sizex, sizey] = size(essential_columns);
    essential_columns_position = ones(1,sizex);
        for ii = 1:sizex
            key = essential_columns{ii};
            find_key = strfind(title_line,key);
            essential_columns_position(ii) = find(tabs_pos>=find_key(1),1);
        end
    %essential_columns_position = sort(essential_columns_position);
end

function [ Box_cell_position ] = find_boxes( title_line )
    % This function returns the column numbers of the columns we are
    % caring, searching by name.
    tabs_pos = find(title_line == 9); % Find were are the "tab" signs are
    tabs_pos = [tabs_pos, (length(title_line)+1)];
    key = 'Box';
    find_key = strfind(title_line,key);
    pp_boxes = length(find_key);
    for pp = 1:pp_boxes
        Box_cell_position(pp) = find(tabs_pos>=find_key(pp),1);
    end
    
    %essential_columns_position = sort(essential_columns_position);
end

function [distlist]=corrDist(x,y)
distlist=zeros(length(x)-length(y)+1,1);
for i=1:length(x)-length(y)+1
    subx=x(i:i+length(y)-1);
    dist=sum((subx-y).^2);
    distlist(i)=dist;
end
end
    
        
function Dvec = movingslope(vec,supportlength,modelorder,dt)
% movingslope: estimate local slope for a sequence of points, using a sliding window

if (nargin==0)
  help movingslope
  return
end
if ~isvector(vec)
  error('vec must be a row or column vector')
end
n = length(vec);
% supply defaults
if (nargin<4) || isempty(dt)
  dt = 1;
end
if (nargin<3) || isempty(modelorder)
  modelorder = 1;
end
if (nargin<2) || isempty(supportlength)
  supportlength = 3;
end
% check the parameters for problems
if (length(supportlength)~=1) || (supportlength<=1) || (supportlength>n) || (supportlength~=floor(supportlength))
  error('supportlength must be a scalar integer, >= 2, and no more than length(vec)')
end
if (length(modelorder)~=1) || (modelorder<1) || (modelorder>min(10,supportlength-1)) || (modelorder~=floor(modelorder))
  error('modelorder must be a scalar integer, >= 1, and no more than min(10,supportlength-1)')
end
if (length(dt)~=1) || (dt<0)
  error('dt must be a positive scalar numeric variable')
end
% now build the filter coefficients to estimate the slope
if mod(supportlength,2) == 1
  parity = 1; % odd parity
else
  parity = 0;
end
s = (supportlength-parity)/2;
t = ((-s+1-parity):s)';
coef = getcoef(t,supportlength,modelorder);
% Apply the filter to the entire vector
f = filter(-coef,1,vec);
Dvec = zeros(size(vec));
Dvec(s+(1:(n-supportlength+1))) = f(supportlength:end);
% patch each end
vec = vec(:);
for i = 1:s
  % patch the first few points
  t = (1:supportlength)' - i;
  coef = getcoef(t,supportlength,modelorder);
  
  Dvec(i) = coef*vec(1:supportlength);
  
  % patch the end points
  if i<(s + parity)
    t = (1:supportlength)' - supportlength + i - 1;
    coef = getcoef(t,supportlength,modelorder);
    Dvec(n - i + 1) = coef*vec(n + (0:(supportlength-1)) + 1 - supportlength);
  end
end
% scale by the supplied spacing
Dvec = Dvec/dt;
% all done
end % mainline end
% =========================================================
% subfunction, used to compute the filter coefficients
function coef = getcoef(t,supportlength,modelorder)
% Note: bsxfun would have worked here as well, but some people
% might not yet have that release of matlab.
A = repmat(t,1,modelorder+1).^repmat(0:modelorder,supportlength,1);
pinvA = pinv(A);
% we only need the linear term
coef = pinvA(2,:);
end % nested function end
