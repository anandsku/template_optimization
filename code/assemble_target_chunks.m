function [target_chunks]=assemble_target_chunks(template_file,template_path,varargin)
%% Syntax
%
% [target_chunks]=assemble_target_chunks(template_file,template_path,input_method,target_file,target_path,varargin)
%
%% Inputs  
% 
% template_file - name of the template file for which you want to assmeble
% target chunks
%
% template_path - location of the template file
%
% input_method - where should the function find the target chunks. two
% options. 'template_metadata': load the syll associated chunks for the length and 
% the index listed in the template metadata. 'file_input' :load 
% the target chunks from a given mat file.
%
% if input_method = 'file_input'
% target_file - name of the file containing target chunks
% 
% target_path - location of the file containing target chunks
% 
% supp_inputs.target_verify - the file containing the target chunks 
% may contain some metadata (depends on if the user put it there while making that file). 
% this metadata allows the function to verify that the target chunks are 
% indeed meant for the template being referred to at the beginning of the function. 
% if this field is set to 1 (the default value), the function will try to
% perform that verification. if it is set to 0, the function will not
% verify. 
%
%
%
%
%% Computation/Processing     
% 
% This function simply loads the template metadata file and the syll assoc
% chunks file. It uses the information stored in the template metadata file
% to load the required chunks in the varible target_chunks. Alternatively
% it loads the file containing the target chunks and performs some basic
% checks on the compatibility between the loaded target chunks and the template. 
% If asked to verfiy (target_verify=1), it also comopares the metadata in
% the target chunks (if any) file to the template metadata to check further for compatibility.  
% 
%
%% Outputs  
% 
% target_chunks - is a matrix of chunks (frequecies are rows) from which
% the template was made. these are the target chunks for this template. 
%
%
%% Assumptions
% Assumes that the template metadata file and the syll_assoc_chunks file
% are both stored at the same location as the template itself. 
%
%
%
% % % Triple percentage sign indicates that the code is part of the code
% template and may be activated if necessary in later versions. 
%% Version and Author Identity Notes  
% 
% Last modified by Anand S. Kulkarni on 
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



if nargin<narg_min
   error(['The number of inputs should at least be ' narg_min])
end

% putting file separators at the end of all input paths
if ~strcmpi(template_path(end),filesep)
    template_path=[template_path,filesep];
end

%% Body of the function
% loading template metadata
[~,metadata_filename,~]=fileparts(template_file);
metadata_fullfile=[template_path metadata_filename '_metadata.mat'];
load(metadata_fullfile) % loads a variable called template_metadata
no_freqs_in_templ=size(template_metadata.template,1);

% loading syllable associated chunks
syll_assoc_chunks_fullfile=[template_path 'syll_assoc_chunks_syll_'...
                        upper(template_metadata.target_syll) '_seq_'...
                        lower(template_metadata.pre_syll) upper(template_metadata.target_syll)...
                        lower(template_metadata.post_syll) '.mat'];
load(syll_assoc_chunks_fullfile); % loads a bunch of variables. syll_assoc_chunks is the one we need

no_freqs_in_syll_assoc_chunks=size(syll_assoc_chunks{1,1}{1,1}{1,1},1);
if no_freqs_in_syll_assoc_chunks~=no_freqs_in_templ
     error('The number of freqs in the loaded chunks and in the template do not match')
end

% obtaining index_of_target_instances and no_target_instances

no_lengths=size(syll_assoc_chunks{1,1},2);
for i=1:no_lengths
    length_in_chunks=size(syll_assoc_chunks{1,1}{1,i}{1,1},2);
    if length_in_chunks==template_metadata.length_in_chunks_of_target_instances
        index_of_target_instances=i;
        no_target_instances=size(syll_assoc_chunks{1,1}{1,i},1);
        break
    end
end
if template_metadata.no_instances_used~=no_target_instances
   warning('There is a discrepancy between the number of instances mentioned in the template metadata and the number obtained from syll_assoc_chunks file'); 
end

% assmbling target chunks
no_target_chunks_per_instance=size(template_metadata.target_chunk_indices,2);
target_chunks=zeros(no_freqs_in_templ,no_target_instances*no_target_chunks_per_instance);
for i=1:no_target_instances
 target_chunks(:,((i-1)*no_target_chunks_per_instance)+1:i*no_target_chunks_per_instance)=...
             syll_assoc_chunks{1,1}{1,index_of_target_instances}{i,1}(:,template_metadata.target_chunk_indices);
end






