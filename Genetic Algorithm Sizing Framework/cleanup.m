% group 6
% clean up script for ga optimization
%--------------------------------------------------------------------------
clc; clear; close all

files1 = dir('Generation*.mat');
if ~isempty(files1)
    delete(files1.name)
end

files2 = dir('Aircraft_WL*.mat');
if ~isempty(files2)
    delete(files2.name)
end
