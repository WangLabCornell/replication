traces = dir("*.mat");

data = [];
for i = 1:length(traces)
    load(traces(i).name)
    item.name = trace_data.output_file_name;
    item.processivity = trace_data.velocity.processivity;

    item.velocity_30s = trace_data.velocity.overall_velocities(1);
    item.torque_30s = trace_data.velocity.overall_torques(1);

    item.velocity_60s = trace_data.velocity.overall_velocities(2);
    item.torque_60s = trace_data.velocity.overall_torques(2);

    item.velocity_120s = trace_data.velocity.overall_velocities(3);
    item.torque_120s = trace_data.velocity.overall_torques(3);

    item.velocity_all = trace_data.velocity.overall_velocities(4);
    item.torque_all = trace_data.velocity.overall_torques(4);

    item.pause_free_velocity = trace_data.velocity.pause_free_velocity;
    item.pause_free_torque = trace_data.velocity.pause_free_torque;
    
   
    data = [data; item];


end

writetable(struct2table(data),"velocity_summary.txt")