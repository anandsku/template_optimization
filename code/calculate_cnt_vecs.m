function [cnt_vecs_aggregate]=calculate_cnt_vecs(dist_vecs_aggregate,curr_thr,amplitudes_aggregate,amp_thr)
%% Syntax
%
% cnt_vecs_aggregate=calculate_cnt_vecs(dist_vecs_aggregate,curr_thr)
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
narg_min=2;

% prob_path=pwd;
% prob_path2=pwd;
% options_available={'option1','option2','option3'};
% 
% % Assigning default values to supplementary inputs
% supp_inputs.sinp1=42;
% supp_inputs.sinp2='everything';
% supp_inputs.sinp3='';
% 
% list=[];
% for i=1:length(options_available)
%     list=[list,options_available{i} '\n'];
% end
% 
% in_message1='Please select a file';
% in_message2='Please select a directory';
% in_message3='Please enter a string input';
% in_message4='Please select a non-string input';
% in_message5='Please select among the following options:';
% in_message6='Please select the files';
if nargin<narg_min
    error(['The number of inputs should at least be ' narg_min])
%      [inp1,inp2]=uigetfile([prob_path filesep '*.abc'],in_message1);   % file input 
%      [inp7]=uipickfiles('filterspec',[prob_path filesep '*.abc'],'prompt',in_message6);  
%      inp3=uigetdir(prob_path2,in_message2); % directory input
%      inp4=input([in_message3 '\n-->  '],'s'); % string input
%      inp5=input([in_message4 '\n-->  ']); % non string input
%      inp6=input([in_message5 '\n' list '\n-->  '],'s');
%      
%      % sets defaults
%      supp_inputs.write_to_disk_q=1; % should the function write a file to disk containing its output  
%      supp_inputs.disk_write_dir=inp2;
%      % Processing supplementary inputs
%      [supp_inputs]=process_supplementary_inputs(supp_inputs);   
% else
%      % sets defaults
%      supp_inputs.write_to_disk_q=1; % should the function write a file to disk containing its output  
%      supp_inputs.disk_write_dir=inp2;
%     % processing supplementary inputs
%     supp_inputs=parse_pv_pairs(supp_inputs,varargin);
end

% packaging the inputs into the inputs structure. This can be useful in
% case you need to store the inputs as meta-data with the output. 
% inputs=struct('inp1',inp1,'inp2',inp2,'inp3',inp3,'inp4',inp4,'inp5',inp5,'inp6',inp6,'inp7',inp7,'spawning_func',mfilename('fullpath'));
% 
% % Checking if output directories need to specified and if they have been specified 
% 
% if supp_inputs.write_to_disk_q
%     if ~exist(supp_inputs.disk_write_dir,'dir')
%         supp_inputs.disk_write_dir=uigetdir(prob_path,'Please select the directory where to store the output mat file. Hit cancel if you don''t want the function to write a mat file');
%         if supp_inputs.disk_write_dir==0
%             supp_inputs.write_to_disk_q=0;
%         end
%     end
% end
% 
% % putting file separators at the end of all input paths
% if ~isempty(inp2)
%     if ~strcmpi(inp2(end),filesep)
%         inp2=[inp2,filesep];
%     end
% end
% if ~isempty(inp3)
%     if ~strcmpi(inp3(end),filesep)
%         inp3=[inp3,filesep];
%     end
% end
% 
% if supp_inputs.write_to_disk_q
%     if ~strcmpi(supp_inputs.disk_write_dir(end),filesep)
%         supp_inputs.disk_write_dir=[supp_inputs.disk_write_dir,filesep];
%     end
% end
% 
% 
% % creating a figure and axes
% fig1=figure;
% axes1=axes('parent',fig1);

%% Body of the function
no_sylls_in_aggregate=length(dist_vecs_aggregate);
cnt_vecs_aggregate(1,no_sylls_in_aggregate).cnts=cell(0);

for i=1:no_sylls_in_aggregate
    no_intervals=length(dist_vecs_aggregate(1,i).dists);
    cnt_vecs_aggregate(1,i).cnts=cell(1,no_intervals);
    for j=1:no_intervals
       no_lengths=length(dist_vecs_aggregate(1,i).dists{1,j});
       cnt_vecs_aggregate(1,i).cnts{1,j}=cell(1,no_lengths);
       for k=1:no_lengths
           no_instances=size(dist_vecs_aggregate(1,i).dists{1,j}{1,k},1);
           slice_length=size(dist_vecs_aggregate(1,i).dists{1,j}{1,k},2);
           cnt_vecs_aggregate(1,i).cnts{1,j}{1,k}=zeros(no_instances,slice_length);
           for l=1:no_instances
               temp_dists=dist_vecs_aggregate(1,i).dists{1,j}{1,k}(l,:);
               temp_amps=amplitudes_aggregate(1,i).amps{1,j}{1,k}(l,:);
               cnt_vec_temp=temp_dists<=curr_thr & temp_amps>amp_thr;
               cnt_vecs_aggregate(1,i).cnts{1,j}{1,k}(l,:)=cnt_vec_temp;               
           end          
       end      
    end   
end
