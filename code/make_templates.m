function []=make_templates(target_syll,distractor_sylls,data_dir,birdID,varargin)
%% Syntax
% 
% []=make_templates(target_syll,pre_syll,post_syll,distractor_sylls,data_dir,birdID,requested_lengths,chunk_indices_for_each_requested_length,varargin)
%
%% Inputs  
%
% target_syll - the symbol of the target syllable for which you want to
% make the template
%
% pre_syll - pre-syllable (if you want to make the template out of instances of the 
% target syllable from a particular sequence)
%
% post_syll - post-syllable (if you want to make the template out of instances of the 
% target syllable from a particular sequence)
% 
% distractor_sylls - symbols of other syllables from the bird's repertoire
% which will serve as distractors
% 
% data_dir - the directory where syll_assoc_chunks file is stored
%
% birdID - the symbol string denoting the bird's identity
%
% requested_lengths - the lengths in chunks of the instances of the target syllable
% from which you want to make the template**1
%
% chunk_indices_for_each_requested_length - the indices of the chunks from
% which you want to make the template**2
% 
%
%% Computation/Processing     
% This function looks for the syll_assoc_chunks file for the syllable and
% the sequence requested. It loads that and picks out all instances of the
% syllable with the requested # of chunks (**1). It then picks up the chunks 
% with the indices given in the input. All these chunks are assembled in
% target_chunks. The template is created by averaging all these chunks.
% 
%
%
% 
%
%% Outputs  
% This function write two files as its output. There is no output provided
% to the calling function. The first file is a .dat file which has the template 
% in the format required by EvTAF. The second file has the same name with
% _metadata.mat appended to the template name. This file contains the metadata for the
% template. Both files are stored in the same directory where the syll_assoc_chunks file
% is found. 
% 
%
%
%% Assumptions
% The template is written by default in the same direcotry as where the
% syll_assoc_chunks file is to be found. 
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
% One needs to run make_syll_assoc_chunks before running this function. 
%
%
%
%% Detailed notes
% **1 - due to variation in duration of syllables, different instances of the syllable
% contain different number of chunks in them. You can pick instances with
% specific lengths to make the template. If you enter 'modal',this function 
% picks the mode of this distribution (of the number of chunks) and works with the instances of
% the target syllable containing the modal number of chunks. this eliminates 
% (somewhat, atleast) the bother about alignment and allows you to average
% over the variation in a given portion of the syllable. 
%
% **2 - if one of your requested length is 12, you may want to make the
% template by using chunks 3 and 4. That is what you would enter here
%
%% Processing inputs and beginning stuff

% putting in a stop for easier debugging
dbstop if error

% processing mandatory inputs
narg_min=4;


if nargin<narg_min
     error(['The number of inputs should at least be ' narg_min])
end


% putting file separators at the end of all input paths
if ~strcmpi(data_dir(end),filesep)
    data_dir=[data_dir,filesep];
    supp_inputs.disk_write_dir=[supp_inputs.disk_write_dir,filesep];
end

%

%% Body of the function

% loading up the syll assoc chunks file
syll_assoc_chunks_fullfile=[data_dir 'syll_assoc_chunks_syll_' upper(target_syll) '_seq_' upper(target_syll) '.mat'];
if ~exist(syll_assoc_chunks_fullfile,'file')
    error(['A syllable associated chunks file for the ' upper(target_syll) ' does not exist at ' data_dir]);
else
    load(syll_assoc_chunks_fullfile) % loads a variable called syll_assoc_chunks 
end

if isempty(syll_assoc_chunks{1,1})
    error('The number of instances of is zero');    
else
    if  isempty(syll_assoc_chunks{1,1}{1,1})
        error('The number of instances of is zero');
    end
end

no_freqs_in_templ=size(syll_assoc_chunks{1,1}{1,1}{1,1},1);

% obtaining no_instances_for_each_length and lengths
no_lengths=size(syll_assoc_chunks{1,1},2);
no_instances_for_each_length=zeros(1,no_lengths);
lengths=zeros(1,no_lengths);
for i=1:no_lengths
    no_instances_for_each_length(1,i)=size(syll_assoc_chunks{1,1}{1,i},1);
    lengths(1,i)=size(syll_assoc_chunks{1,1}{1,i}{1,1},2);
end

% obtaining indices_for_requested_lengths, requested_lengths (re-obtaining), and no_lengths_requested 
[~,max_instances_ind]=max(no_instances_for_each_length);
modal_length=lengths(max_instances_ind(1));
requested_lengths=modal_length;
no_lengths_requested=length(requested_lengths); % this of course, has to be one, since its modal
indices_for_requested_lengths=max_instances_ind;



% temp variables for the following operation
temp_indices_for_requested_lengths=indices_for_requested_lengths;
temp_chunk_indices_for_each_requested_length=cell(1);

new_i=1;
for i=1:no_lengths_requested
   no_all_chunks=lengths(indices_for_requested_lengths(1,i));
   temp_indices_for_requested_lengths=[temp_indices_for_requested_lengths(1,1:new_i),repmat(temp_indices_for_requested_lengths(1,new_i),1,no_all_chunks-1),temp_indices_for_requested_lengths(1,new_i+1:end)];
   temp_chunk_indices=cell(1,no_all_chunks);
   for j=1:no_all_chunks
        temp_chunk_indices{1,j}=j;
   end
   temp_chunk_indices_for_each_requested_length=[temp_chunk_indices_for_each_requested_length{1,1:new_i-1},temp_chunk_indices,temp_chunk_indices_for_each_requested_length{1,new_i+1:end}];
   new_i=i+no_all_chunks;
  
end

indices_for_requested_lengths=temp_indices_for_requested_lengths;
chunk_indices_for_each_requested_length=temp_chunk_indices_for_each_requested_length;

% over the course of previous code we have obtained: indices_for_requested_lengths
% and  chunk_indices_for_each_requested_length. we will use these going ahead. 
no_lengths_requested=length(indices_for_requested_lengths);
for i=1:no_lengths_requested % equivalent to no_lengths_requested
    no_instances_this_length=no_instances_for_each_length(1,indices_for_requested_lengths(1,i));
    no_target_chunks_per_instance=length(chunk_indices_for_each_requested_length{1,i});       
    target_chunks_for_template=zeros(no_freqs_in_templ,no_instances_this_length*no_target_chunks_per_instance);
       for k=1:no_instances_this_length
         target_chunks_for_template(:,((k-1)*no_target_chunks_per_instance)+1:k*no_target_chunks_per_instance)=syll_assoc_chunks{1,1}{1,indices_for_requested_lengths(1,i)}{k,1}(:,chunk_indices_for_each_requested_length{1,i});
       end
   template=mean(target_chunks_for_template,2);
   template=template-min(template);
   template=template./max(template);

   template_name=['template_syll_' upper(target_syll) '_seq_' target_syll ...
                  '_chunks_'  num2str(chunk_indices_for_each_requested_length{1,i})...
                  '_outof_' num2str(lengths(indices_for_requested_lengths(1,i)))];
   template_name=strrep(template_name,'  ','_'); % replacing double space with an underscore 
   fid=fopen([data_dir template_name '.dat'],'w');
    for j=1:size(template,1)
        fprintf(fid,'%.5e\n',template(j,1));
    end
    fclose(fid);
    template_metadata.target_syll=target_syll;
    template_metadata.pre_syll='';
    template_metadata.post_syll='';
    template_metadata.distractor_sylls=distractor_sylls;
%         template_metadata.file_list=arch_inputs.file_list;
    template_metadata.template=template;
    template_metadata.length_in_chunks_of_target_instances=lengths(indices_for_requested_lengths(1,i));
    template_metadata.target_chunk_indices=chunk_indices_for_each_requested_length{1,i};
    template_metadata.no_instances_used=no_instances_for_each_length(indices_for_requested_lengths(1,i));
    template_metadata.birdID=birdID;
    template_metadata.is_template_synthetic=0;
    template_metadata.synthesis_details=[];

    template_metadata.name_of_struct='template_metadata';
%     check_struct_against_blank(template_metadata,@create_blank_datastruct_template_metadata);

    matfile=[template_name '_metadata.mat'];
    matfullfile=[data_dir matfile];
    save(matfullfile,'template_metadata');
       
end


