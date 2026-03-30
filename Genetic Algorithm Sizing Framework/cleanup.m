% group 6
% clean up script for ga optimization
%--------------------------------------------------------------------------
clc; clear; close all

files1 = dir('Veh_WL*');
if ~isempty(files1)
delete(files1.name)
end
files2 = dir('Generation*.mat');
if ~isempty(files2)
delete(files2.name)
end
folders1 = dir("Generation*");

folders1 = folders1([folders1.isdir]);

if ~isempty(folders1)
for i = [string(folders1.name)]
    rmdir(i)
end
end