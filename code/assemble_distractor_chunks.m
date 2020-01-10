function [distractor_chunks]=assemble_distractor_chunks(template_file,template_path,varargin)
%% Syntax
%
% [distractor_chunks]=assemble_distractor_chunks(template_file,template_path,...
%                   input_method,distract_file,distract_path,varargin)
%% Inputs  
%
% template_file - name of the template file
% 
% template_path - location of the template file
%
% input_method - where should the function find the distractor chunks. two
% options. 'template_metadata': load the syll associated chunks for all the 
% listed distractor syllables in the template metadata. 'file_input' :load 
% the distractor chunks from a given mat file.
%
% if input_method = 'file_input'
% distract_file - name of the file containing distractor chunks
% 
% distract_path - location of the file containing distractor chunks
%
% supp_inputs.missing_distract_files -  what to do in case the syll assoc
% files for some of the distractor syllables give an error while loading.
% two options: 'throw_error' - this is the default value. the code will throw an
% error. 'ignore' - the code will ignore the error encountered while
% loading the file and will continue to the next distractor syllable in the
% list. 
%
% supp_inputs.distract_verify - the file containing the distractor chunks 
% may contain some metadata (depends on if the user put it there while making that file). 
% this metadata allows the function to verify that the distractor chunks are 
% indeed meant for the template being referred to at the beginning of the function. 
% if this field is set to 1 (the default value), the function will try to
% perform that verification. if it is set to 0, the function will not
% verify. 
%
%% Computation/Processing     
% assembles distractor chunks into variable called distractor_chunks
%
%
% 
%
%% Outputs  
%  distractor_chunks -  matrix containing the distractor chunks. the rows
%  represent frequencies and the columns are the different chunks. 
% 
%
%
%% Assumptions
%
%  The function assumes that the syll_assoc_chunks are located at the same
%  path as the template
%
%
% % % Triple percentage sign indicates that the code is part of the code
% template and may be activated if necessary in later versions. 
%% Version and Author Identity Notes  
% 
% Last modified by Anand Kulkarni on 
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
narg_min=5;

% Assigning default values to supplementary inputs
supp_inputs.missing_distract_files='throw_error';
supp_inputs.exclude_distractor_syll='';
supp_inputs.include_gaps=1;

if nargin<narg_min
     error(['The number of inputs should at least be ' narg_min])
else
    % processing supplementary inputs
    if ~isempty(varargin)
        if iscell(varargin{1})
            if strcmpi(varargin{1}{1},'unpackpvpairs')
                 varargin=varargin{1}(2:end);
            end
        end
    end
    supp_inputs=parse_pv_pairs(supp_inputs,varargin);    
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

distractor_sylls=template_metadata.distractor_sylls;
no_distract_sylls=size(distractor_sylls,2);
distractor_chunks=zeros(no_freqs_in_templ,1);
last_chunk_no=1;
% loading distractor_chunks
for i=1:no_distract_sylls
    distr_syll=distractor_sylls(1,i);
    if ~isempty(supp_inputs.exclude_distractor_syll)
       if isequal(supp_inputs.exclude_distractor_syll,distr_syll) 
          continue 
       end           
    end
    syll_assoc_chunks_fullfile=[template_path 'syll_assoc_chunks_syll_'...
                            upper(distr_syll) '_seq_' upper(distr_syll) '.mat'];
    try 
        load(syll_assoc_chunks_fullfile); % loads a bunch of variables. syll_assoc_chunks is the one we need
    catch exp
        if strcmpi(exp.identifier,'MATLAB:load:couldNotReadFile')
            if strcmpi(supp_inputs.missing_distract_files,'ignore')                    
                continue
            else
                rethrow(exp);
            end                
        else
           rethrow(exp);
        end
    end
    if isempty(syll_assoc_chunks{1,1})|| isempty(syll_assoc_chunks{1,1}{1,1})             
        warning(['Distractor syllable ' distr_syll ' syll_assoc_chunks is empty']);
        continue
    end

    no_freqs_in_syll_assoc_chunks=size(syll_assoc_chunks{1,1}{1,1}{1,1},1);




    if no_freqs_in_syll_assoc_chunks~=no_freqs_in_templ
         error('The number of freqs in the loaded chunks and in the template do not match')
    end

    straightened_chunks=straighten_syll_assoc_chunks(syll_assoc_chunks);
    no_st_chunks=size(straightened_chunks,2);
    distractor_chunks(:,last_chunk_no:last_chunk_no+no_st_chunks-1)=straightened_chunks;
    last_chunk_no=last_chunk_no+no_st_chunks;
end
if supp_inputs.include_gaps==1
     gap_assoc_chunks_fullfile=[template_path 'gap_assoc_high_amp_chunks.mat'];
     load(gap_assoc_chunks_fullfile); % loads variable gap_assoc_chunks
     distractor_chunks=[distractor_chunks,gap_assoc_high_amp_chunks];
end
    

if exist([template_path 'gap_assoc_high_amp_chunks.mat'],'file')
    load([template_path 'gap_assoc_high_amp_chunks.mat']); % loads variable called   gap_assoc_high_amp_chunks
    distractor_chunks=[distractor_chunks,gap_assoc_high_amp_chunks];
else
   warning('Gap assoc high amplitude chunks file not found. Not loading those chunks in distractors') 
   
end

