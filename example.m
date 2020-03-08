clc; clear; format;
%% test of least-square method
obj0 = powergrid_state_estimation('iSE30Bus.txt');
obj0.infer('least-square', 100, 1e-7, 'oStateEstimation.txt');

%% test of fast-decoupled method
obj1 = powergrid_state_estimation('iSE30Bus.txt');
obj1.infer('fast-decoupled', 100, 1e-5, 'oFDSE.txt');