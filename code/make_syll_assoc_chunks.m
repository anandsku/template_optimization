function [syll_assoc_chunks]=make_syll_assoc_chunks(batchfile,batchpath,target_syll)
%% Syntax
%
% [syll_assoc_chunks]=make_syll_assoc_chunks_file(batchpath,batchfile,target_syll,pre_syll,post_syll,varargin)
%
%% Inputs  
%
% batchpath -  location where the batch file is stored
% 
% batchfile -  name of the batch file
% 
% target_syll  - syllable for which you want to make a syllable associated
% chunks file
% 

% 
%% Computation/Processing     
%  This function will use the *.notmat and *.simdata files for the files listed
%  in the batchfile and find out all instances of the target_syll. For each
% instance, it will pickup all associated chunks. It will sort instances by
% the number of chunks in them and write them in syll_assoc_chunks. 
%
%
% 
%
%% Outputs  
% It writes a syll_assoc_chunks file containing the following cell array. 
% syll_assoc_chunks - is a cell array. it contains within itself a cell array,
% where each cell contains all instances with the same number
% of chunks.   
% 
%
%
%% Assumptions
% Assumes that nomat files and simdata files exist in the location as the batchfile.  
%
%
%
%
%% Version and Author Identity Notes  
% 
% Last modified by Anand S Kulkarni 
%  
%% Related procedures and functions 
% 
%
%
%
%% Detailed notes
% This function is capable of getting chunks associated with multiple intervals. 
% If the target_syll is entered as ab, it will give you all the asscoiated
% chunks for three intervals a, gap between a and b and b. 
%
%
%
%% Processing inputs and beginning stuff

% putting in a stop for easier debugging
dbstop if error

% processing mandatory inputs
narg_min=3;

if nargin<narg_min
     error(['The number of inputs should at least be ' narg_min])
end


% putting file separators at the end of all input paths
if ~strcmpi(batchpath(end),filesep)
    batchpath=[batchpath,filesep];
    supp_inputs.disk_write_dir=[supp_inputs.disk_write_dir,filesep];
end



%
%% Body of the function

fid=fopen([batchpath batchfile],'r');
file=fgetl(fid);
chunks=cell(0);
amplitudes=cell(0);
no_chunks_vec=[];
durations_vec=[];
onsets_vec=[];


no_missed_instances=0;
tot_no_instances=0;

while ischar(file)
    
    [~,name,ext]=fileparts(file);    
    notmat_file=[batchpath name ext '.not.mat'];
    simdata_file=[batchpath name '_simdata.mat'];
    load(notmat_file)
    load(simdata_file)
    
   inds=strfind(lower(labels),lower(target_syll));
   if exist('labels2','var')
       broken_sylls=find(labels2~=1);
       corrupt_sylls=find(corrupted==1);
       inds=setdiff(inds,[broken_sylls,corrupt_sylls]);    
       [parsed_labels,~,~]=parse_complex_label(labels,labels2,labels3,subclips);
       all_inds=strfind(lower(parsed_labels),lower(target_syll));
       no_missed_instances=no_missed_instances+length(all_inds)-length(inds);
       tot_no_instances=tot_no_instances+length(all_inds);
   else
      no_missed_instances=no_missed_instances+0; 
      tot_no_instances=tot_no_instances+length(inds);
   end
   no_isntances_in_file=length(inds);
   temp_chunks=cell(length(inds),1);
   temp_no_chunks_vec=zeros(length(inds),1);
   temp_durations=zeros(length(inds),1);
   temp_onsets=zeros(length(inds),1);
   temp_amplitudes=cell(length(inds),1);

   for j=1:no_isntances_in_file         
       interval_ind=inds(j);

       assoc_chunks=[];
       for k=1:length(simdata.chunk_assocs)
           if (simdata.chunk_assocs{k})==interval_ind 
               assoc_chunks=[assoc_chunks,k];
           end
       end 

      temp_chunks{j,1}=simdata.chunks(:,assoc_chunks);
      temp_amplitudes{j,1}=simdata.chunk_amp(assoc_chunks);
      temp_no_chunks_vec(j,1)=length(assoc_chunks);   
      temp_durations(j,1)=offsets(interval_ind)-onsets(interval_ind);
      temp_onsets(j,1)=(simdata.cfp_times(assoc_chunks(1))*10^3)+(simdata.chunk_duration*10^3)-onsets(interval_ind);  
      % the onset is relative to the first slice offset
       
   end
   % looking for empy temp_chunks
   empty_inds=find(temp_no_chunks_vec==0);
   temp_chunks(empty_inds)=[];
   temp_amplitudes(empty_inds)=[];
   temp_no_chunks_vec(empty_inds)=[];
   temp_durations(empty_inds)=[];
   temp_onsets(empty_inds)=[];
   
%    if ~isempty(find(temp_onsets<0 | temp_onsets>2*(simdata.chunk_duration*10^3), 1))
%       disp('this is wrong')        
%    end
   
   chunks=[chunks;temp_chunks];
   amplitudes=[amplitudes;temp_amplitudes];
   no_chunks_vec=[no_chunks_vec;temp_no_chunks_vec];
   durations_vec=[durations_vec;temp_durations];
   onsets_vec=[onsets_vec;temp_onsets];
   file=fgetl(fid) ;
end

fclose(fid);

no_instances=size(chunks,1);
unique_lengths=cell(1,1);
syll_assoc_chunks=cell(1,1);
syll_assoc_durations=cell(1,1);
syll_assoc_onsets=cell(1,1);
syll_assoc_amplitudes=cell(1,1);


[lengths_vec,~]=count_unique(no_chunks_vec(:,1));
unique_lengths{1,1}=lengths_vec;    

syll_assoc_chunks{1,1}=cell(1,length(unique_lengths{1,1}));
syll_assoc_durations{1,1}=cell(1,length(unique_lengths{1,1}));
syll_assoc_amplitudes{1,1}=cell(1,length(unique_lengths{1,1}));
syll_assoc_onsets{1,1}=cell(1,length(unique_lengths{1,1}));

for k=1:no_instances 
    instance_chunks=chunks{k,1};
    no_chunks_in_instance=size(instance_chunks,2);
    instance_duration=durations_vec(k,1);
    instance_onset=onsets_vec(k,1);
    instance_amplitudes=amplitudes{k,1};

    len_index=find(unique_lengths{1,1}==no_chunks_in_instance);
    syll_assoc_chunks{1,1}{1,len_index}=[syll_assoc_chunks{1,1}{1,len_index};{instance_chunks}];
    syll_assoc_durations{1,1}{1,len_index}=[syll_assoc_durations{1,1}{1,len_index};instance_duration];
    syll_assoc_onsets{1,1}{1,len_index}=[syll_assoc_onsets{1,1}{1,len_index};instance_onset];
    syll_assoc_amplitudes{1,1}{1,len_index}=[syll_assoc_amplitudes{1,1}{1,len_index};instance_amplitudes];
end



%
%% Processing outputs and ending stuff

matfile=['syll_assoc_chunks_syll_' upper(target_syll) '_seq_' upper(target_syll) '.mat'];
matfullfile=[batchpath matfile];
save(matfullfile,'syll_assoc_chunks','syll_assoc_durations','tot_no_instances',...
        'no_missed_instances','syll_assoc_amplitudes','syll_assoc_onsets');

% removing the stop that was put for easier debugging
dbclear if error



