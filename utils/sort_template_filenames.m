function [filenames]=sort_template_filenames(filenames)
%% Syntax
%
% [filenames]=sort_template_filenames(filenames)
%
%% Inputs  
%
%
%
%
%% Computation/Processing     
% 
%
%
% 
%
%% Outputs  
% 
% 
%
%
%% Assumptions
%
%
%
%
% % % Triple percentage sign indicates that the code is part of the code
% template and may be activated if necessary in later versions. 
%% Version and Author Identity Notes  
% 
% Last modified by The Big Foot on 1/1/1400
% 
% previous version:
% next version: 
%% Related procedures and functions 
% 
%
%
%
%% Detailed notes
%
%
%
%
%% Processing inputs and beginning stuff

% putting in a stop for easier debugging
dbstop if error

% processing mandatory inputs
narg_min=1;

if nargin<narg_min
    error(['The number of inputs should at least be ' narg_min])
end

%% Body of the function
% the template # if the # between 6th and 7th undersocre
% 'template_syll_G_seq_G_chunks_10_outof_12.dat or template_syll_G_seq_G_chunks_10_outof_12_optimized.dat
template_nos=zeros(1,length(filenames));
for i=1:length(filenames)
    underscores=strfind(filenames{i},'_');
    if ~ismember(length(underscores),[8,9])
        error('The file name is not in a proper format')
    end
    temp_no=str2num(filenames{i}(underscores(6)+1:underscores(7)-1));
    if isempty(temp_no)|| ~isnumeric(temp_no)
        error('The file name is not in a proper format');
    end
    template_nos(1,i)=temp_no;   
end

[~,sorted_inds]=sort(template_nos);
filenames=filenames(sorted_inds);

