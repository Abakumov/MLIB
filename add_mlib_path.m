
%% ADD_MLIB_PATH.m
% adds MLIB library and its subfolders to the MATLAB search path 
%
%% Author: Ivan Abakumov
% Freie Universität Berlin
% E-mail: abakumov_ivan@mail.ru
% Publication date: 12th July 2019

addpath(genpath(mlibfolder))

% Figures always open with the tools
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))