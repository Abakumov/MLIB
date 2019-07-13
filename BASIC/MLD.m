function data = MLD(filename)

% My Load - load mat files
%% Author: Ivan Abakumov
% Freie Universität Berlin
% E-mail: abakumov_ivan@mail.ru
% Publication date: 12th July 2019

tempdata=load(filename);   
field=cell2mat(fieldnames(tempdata));
data=getfield(tempdata, field); %#ok<GFLD>