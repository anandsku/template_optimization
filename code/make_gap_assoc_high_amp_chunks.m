function []=make_gap_assoc_high_amp_chunks(batchfile,batchpath,amp_thr)
%% Syntax
%
% [outputs]=function_template(inp1,inp2,inp3,inp4,varargin)
%
%% Inputs  
% batchfile - name of bacth file
% batchpath -  location of the batch file
% amp_thr - value of amplitude threshold for this bird
%
%
%
%% Computation/Processing     
% It collated all gap assigned slices that have amplitide higher than the
% threshold and writes a file containing those slices. 
%
%
% 
%
%% Outputs  
% gap_assoc_high_amp_chunks - is just a matrix where columns are the
% various gap slices with high amplitude
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
narg_min=3;

prob_path=pwd;
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

in_message1='Please select the batch file';
% in_message2='Please select a directory';
% in_message3='Please enter a string input';
in_message4='Please select the value of the amplitude threshold';
% in_message5='Please select among the following options:';
% in_message6='Please select the files';
if nargin<narg_min
    % error(['The number of inputs should at least be ' narg_min])
     [batchfile,batchpath]=uigetfile([prob_path filesep '*.*'],in_message1);   % file input 
%      [inp7]=uipickfiles('filterspec',[prob_path filesep '*.abc'],'prompt',in_message6);  
%      inp3=uigetdir(prob_path2,in_message2); % directory input
%      inp4=input([in_message3 '\n-->  '],'s'); % string input
     amp_thr=input([in_message4 '\n-->  ']); % non string input
%      inp6=input([in_message5 '\n' list '\n-->  '],'s');
%      
%      % sets defaults
%      supp_inputs.write_to_disk_q=1; % should the function write a file to disk containing its output  
%      supp_inputs.disk_write_dir=batchfile;
%      supp_inputs.input_fields_clash_decision='current';
%      % Processing supplementary inputs
%      [supp_inputs]=process_supplementary_inputs(supp_inputs);   
% else
%      % sets defaults
%      supp_inputs.write_to_disk_q=1; % should the function write a file to disk containing its output  
%      supp_inputs.disk_write_dir=batchfile;
%      supp_inputs.input_fields_clash_decision='current';
%     % processing supplementary inputs
%     supp_inputs=parse_pv_pairs(supp_inputs,varargin);
end

% packaging the inputs into the inputs structure. This can be useful in
% case you need to store the inputs as meta-data with the output. 
% inputs=struct('inp1',batchpath,'inp2',batchfile,'inp3',inp3,'inp4',inp4,'inp5',inp5,'inp6',inp6,'inp7',inp7,'spawning_func',mfilename('fullpath'));

% Checking if output directories need to specified and if they have been specified 

% if supp_inputs.write_to_disk_q
%     if ~exist(supp_inputs.disk_write_dir,'dir')
%         supp_inputs.disk_write_dir=uigetdir(prob_path,'Please select the directory where to store the output mat file. Hit cancel if you don''t want the function to write a mat file');
%         if supp_inputs.disk_write_dir==0
%             supp_inputs.write_to_disk_q=0;
%         end
%     end
% end

% putting file separators at the end of all input paths
if ~isempty(batchpath)
    if ~strcmpi(batchpath(end),filesep)
        batchpath=[batchpath,filesep];
    end
end
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
filelist=make_filelist_from_batch(batchfile,batchpath);
no_files=size(filelist,1);
gap_assoc_high_amp_chunks=[];




for i=1:no_files
    file=filelist{i};
    [~,name,ext]=fileparts(file);    
    notmat_file=[batchpath name ext '.not.mat'];
    simdata_file=[batchpath name '_simdata.mat'];
    load(notmat_file)
    load(simdata_file)
    lens=cellfun(@length,simdata.chunk_assocs);
    gap_chunks=simdata.chunks(:,lens==2);
    gap_chunks_amp=simdata.chunk_amp(lens==2);
    high_amp_chunks=gap_chunks(:,gap_chunks_amp>amp_thr);
    gap_assoc_high_amp_chunks=[gap_assoc_high_amp_chunks,high_amp_chunks];    
end

save([batchpath 'gap_assoc_high_amp_chunks.mat'],'gap_assoc_high_amp_chunks')



% dummy plotting 
% plot(axes1,1,1,'or')
% 
% %% Processing outputs 
% arch_timestamp=datestr(now,'yyyy-mmm-dd HH:MM:SS');
% inputs=setfield('var1',var1); % any additional inputs collected during the function
% primitive_input_fields=fieldnames(arch_inputs);
% curr_input_fields=fieldnames(inputs);
% for i=1:length(primitive_input_fields)
%     fld=primitive_input_fields{i};
%     if ~isempty(find(strcmpi(fld,curr_input_fields),1))
%         switch supp_inputs.input_fields_clash_decision            
%             case 'current'               
%                 warning(['Archival input field ' fld ' is the same for inputs and arch_inputs. Going with the value in inputs'])
%             case 'primitive'                
%                 inputs=setfield(inputs,fld,arch_inputs.(fld)); 
%                 warning(['Archival input field ' fld ' is the same for inputs and arch_inputs. Going with the value in arch_inputs'])
%             otherwise                
%                 error('Incorrect input for var supp_inputs.input_fields_clash_decision')
%         end
%     else
%         inputs=setfield(inputs,fld,arch_inputs.(fld)); 
%     end
% end
% 
% 
% 
% arch_inputs=inputs;
% arch_supp_inputs=supp_inputs;
% 
% if supp_inputs.write_to_disk_q==1
%     matfile='function_output.mat';
%     matfullfile=[supp_inputs.disk_write_dir  matfile];
%     save(matfullfile,'outvar');
%     % alternative save: save(matfullfile,'outvar','arch_inputs','arch_supp_inputs','arch_timestamp');
%     figfile='function_plot.fig';
%     figfullfile={[supp_inputs.disk_write_dir  figfile]};
%     saveas(fig1,figfullfile);
% end

