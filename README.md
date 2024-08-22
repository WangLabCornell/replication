# replication

This repository contains a collection of scripts that extract statistics from single molecule traces.

1. AOT_replisome_data_analysis.m
Script to load raw traces from the Angular Optical Trap (AOT) and calculate the fork position from the extension. This script plots the trace and can measure the stall torque of each trace.  Last, regions of the data are collected and export for further analysis, below.

2. replication_position_vs_time_analysis.m
Script to process replication traces to extract the mean velocity and pause free velocity of the replisome.

3. calculate_velocity_position.m
Function that is used to calculate the velocity and smoothed position from the position vs time data.  Used in (2)

4. MT_stall_and_restart_duration.m
Script to process single molecule traces from a Magnetic tweezers experiment.  Measures and exports timing information about when the replisome stalls and/or restarts.
