function []=write_optimized_templates_from_optimization_results(filename,filepath)
%% Syntax
%
% [outputs]=function_template(inp1,inp2,inp3,inp4,varargin)
%
%% Inputs  
% filename - name of the optimization results file
% filepath - location of the said file
%
%
%
%% Computation/Processing     
% It simply extracts the optimized template from the optimization results
% data structure and writes them as template files. 
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
narg_min=2;

prob_path=pwd;

in_message1='Please select the optimization grand results file';

if nargin<narg_min
     [filename,filepath]=uigetfile([prob_path filesep '*.mat'],in_message1);   % file input 
end

% putting file separators at the end of all input paths
if ~isempty(filepath)
    if ~strcmpi(filepath(end),filesep)
        filepath=[filepath,filesep];
    end
end

%% Body of the function
load([filepath filename])% laods a variable called results among others
results_opt=results;
no_templates=length(results_opt);

for i=1:no_templates
    template_file=results_opt(1,i).templatefullfile;
    [~,nm,~]=fileparts(template_file);
    template_metadata_file=[nm '_metadata.mat'];
    optimized_template=results_opt(1,i).optim_archive.template_vec(:,end);
    wrt_templ([filepath nm '_optimized'],optimized_template,0)
    copyfile([filepath template_metadata_file],[filepath nm '_optimized_metadata.mat'])
    load([filepath nm '_optimized_metadata.mat']) % loads varaible called template_metadata
    template_metadata.template=optimized_template;
    save([filepath nm '_optimized_metadata.mat'],'template_metadata','-append')
end


