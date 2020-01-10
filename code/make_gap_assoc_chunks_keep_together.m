function make_gap_assoc_chunks_keep_together(batchfile,batchpath,varargin)
%% Inputs:
% batchfile - name of batch file
% batchpath - location of the bacth file

%% Outputs:
% gap_assoc_chunks - it conatins all gaps assigned slices. each gap is kept
% as a separate matrix and all such matrices are stored in a cell array

%% Processing inputs and beginning stuff

% putting in a stop for easier debugging
dbstop if error

% processing mandatory inputs
narg_min=2;

prob_path=pwd;

in_message1='Please select the batch file with the list of files you want to use';
if nargin<narg_min
     [batchfile,batchpath]=uigetfile([prob_path filesep 'batch*'],in_message1);    
end

% packaging the inputs into the inputs structure. This can be useful in
% case you need to store the inputs as meta-data with the output. 
inputs=struct('batchfile',batchfile,'batchpath',batchpath,'spawning_func',mfilename('fullpath'));


% processing supplementary inputs

% Assigning default values to supplementary inputs
supp_inputs.disk_write_dir=batchpath;
supp_inputs.write_to_disk_q=1; % should the function write a mat file to disk containing its output  

supp_inputs=parse_pv_pairs(supp_inputs,varargin);

% Checking if output directories need to specified and if they have been specified 

if supp_inputs.write_to_disk_q
    if ~exist(supp_inputs.disk_write_dir,'dir')
        supp_inputs.disk_write_dir=uigetdir('Please select the directory where to store the output mat file. Hit cancel if you don''t want the function to write a mat file');
        if supp_inputs.disk_write_dir==0
            supp_inputs.write_to_disk_q=0;
        end
    end
end


% putting file separators at the end of all input paths
if ~strcmpi(batchpath(end),filesep)
    batchpath=[batchpath,filesep];
    supp_inputs.disk_write_dir=[supp_inputs.disk_write_dir,filesep];
end
%% Body of the function

fid=fopen([batchpath batchfile],'r');
file=fgetl(fid);
% chunks=cell(0);
no_chunks_vec=[];
amplitudes_vec=cell(0);
% durations_vec=[];
file_list=cell(0);
chunks=cell(0);
% 
% no_intervals=2*length(target_syll)-1;
% srch_string=[pre_syll target_syll post_syll];
% no_missed_instances=0;
% tot_no_instances=0;

while ischar(file)
    [~,name,ext]=fileparts(file);    
    notmat_file=[batchpath name ext '.not.mat'];
    simdata_file=[batchpath name '_simdata.mat'];
    load(notmat_file)
    load(simdata_file)
    file_list=[file_list;file];
    
%    inds=strfind(lower(labels),lower(srch_string));
%    if exist('labels2','var')
%        broken_sylls=find(labels2~=1);
%        corrupt_sylls=[];
%        if supp_inputs.exclude_corrupt==1
%             corrupt_sylls=find(corrupted==1);
%        elseif supp_inputs.exclude_corrupt~=0
%            error('The value of supp_inputs.exclude_corrupt should be only 1 or 0')
%        end
%        inds=setdiff(inds,[broken_sylls,corrupt_sylls]);    
%        [parsed_labels,~,~]=parse_complex_label(labels,labels2,labels3,subclips);
%        all_inds=strfind(lower(parsed_labels),lower(srch_string));
%        no_missed_instances=no_missed_instances+length(all_inds)-length(inds);
%        tot_no_instances=tot_no_instances+length(all_inds);
%    else
%       no_missed_instances=no_missed_instances+0; 
%       tot_no_instances=tot_no_instances+length(inds);
%    end
%    no_isntances_in_file=length(inds);
%    temp_chunks=cell(length(inds),no_intervals);
%    temp_no_chunks_vec=zeros(length(inds),no_intervals);
%    temp_durations=zeros(length(inds),no_intervals);

%    for j=1:no_isntances_in_file
% 
%        for l=1:no_intervals
% 
%            if mod(l,2)==0
%                interval_ind=[inds(j)+l/2-1,inds(j)+l/2]+length(pre_syll);
%            else
%                interval_ind=inds(j)+((l+1)/2)-1+length(pre_syll);
%            end
           nosylls=length(labels);
           for idk=1:nosylls+1
               if idk==1
                  cup=[-1,1]; 
               elseif idk==nosylls+1
                   cup=[nosylls,-1];
               else
                   cup=[idk-1,idk];
               end

               assoc_chunks=[];

               for k=1:length(simdata.chunk_assocs)
                   if simdata.chunk_assocs{k}==cup 
                       assoc_chunks=[assoc_chunks,k];
                   end
               end 
               if ~isempty(assoc_chunks)
                    chunks=[chunks;simdata.chunks(:,assoc_chunks)];
                    no_chunks_vec=[no_chunks_vec;length(assoc_chunks)];
                    amplitudes_vec=[amplitudes_vec;simdata.chunk_amp(:,assoc_chunks)];
               end
               
                  
           end

          
%           temp_no_chunks_vec(j,l)=length(assoc_chunks);   
%           temp_durations(j,l)=offsets(interval_ind)-onsets(interval_ind);
%        end
%    end
   % looking for empy temp_chunks
%    empty_inds=find(temp_no_chunks_vec==0);
%    temp_chunks(empty_inds)=[];
%    temp_no_chunks_vec(empty_inds)=[];
%    temp_durations(empty_inds)=[];
%    
%    chunks=[chunks;temp_chunks];
%    no_chunks_vec=[no_chunks_vec;temp_no_chunks_vec];
%    durations_vec=[durations_vec;temp_durations];
   file=fgetl(fid) ;
end

fclose(fid);

no_instances=size(chunks,1);
unique_lengths=cell(1,1);
syll_assoc_chunks=cell(1,1);
syll_assoc_amplitudes=cell(1,1);


    [lengths_vec,~]=count_unique(no_chunks_vec);
    unique_lengths{1,1}=lengths_vec;    

    syll_assoc_chunks{1,1}=cell(1,length(unique_lengths{1,1}));
    syll_assoc_amplitudes{1,1}=cell(1,length(unique_lengths{1,1}));

    for k=1:no_instances 
        instance_chunks=chunks{k,1};
       instance_amplitudes=amplitudes_vec{k,1};

        no_chunks_in_instance=size(instance_chunks,2);

        len_index=find(unique_lengths{1,1}==no_chunks_in_instance);
        syll_assoc_chunks{1,1}{1,len_index}=[syll_assoc_chunks{1,1}{1,len_index};{instance_chunks}];
        syll_assoc_amplitudes{1,1}{1,len_index}=[syll_assoc_amplitudes{1,1}{1,len_index};instance_amplitudes];
    end

gap_assoc_chunks=syll_assoc_chunks;
gap_assoc_amplitudes=syll_assoc_amplitudes;

%
%% Processing outputs and ending stuff
arch_timestamp=datestr(now,'yyyy-mmm-dd HH:MM:SS');
inputs.file_list=file_list;
arch_inputs=inputs;
arch_supp_inputs=supp_inputs;

if supp_inputs.write_to_disk_q==1
    matfile=['gap_assoc_chunks_kept_together.mat'];
    matfullfile=[supp_inputs.disk_write_dir matfile];
    save(matfullfile,'gap_assoc_chunks','arch_inputs','arch_supp_inputs','arch_timestamp','gap_assoc_amplitudes');
end

% removing the stop that was put for easier debugging
dbclear if error



