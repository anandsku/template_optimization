function []=create_simdata_file(batchname,batchpath,varargin)
%% Inputs:
% batchname - Name of the batch file
% batchpath - Location  of the batch file
% 
% supp_inputs.chanspec - specify the channel to be read from the data file (useful for multichannel files)
% supp_inputs.nfft - size (in sample points) of the data slices 
% supp_inputs.overwrite - overwrite or not, an existing siumdata file

%% Outputs:
% It outputs a variable called simdata in a file ending in _simdata.mat. Simdata has the following fields:
% 'first_chunk_first_point' - The data point at which the code starts the
% first chunk/slice

% 'cbinfile' - name of the raw data file               
% 'chunks' - contains the filtered and normalized spectra of the chunks                  
% 'specs' - contains the raw absolute spectrum                  
% 'chunk_amp' - amplitude of the chunks               
% 'no_chunks'               
% 'chunk_first_points' - first data points in each chunk (samples)     
% 'cfp_times' - first data points in each chunk (times)          
% 'chunk_duration' - length of chunks (time)         
% 'chunk_len' - length of chunks (samplea)                     
% 'chunk_assocs' - denotes whether a chunk is associated with syllable or gap. if it 
% associated with a syllable, it denotes the syllable serial number

% 'chunk_assoc_type' - description of the chunk association       
% 'pc_in_syll' - % of the chunk that lies inside a syllable             
% 'templates' - ignore              
% 'simul' - ignore                  
% 'detect_settings_n_rates' - ignore

%%
% putting in a stop for easier debugging
dbstop if error

% processing mandatory inputs
narg_min=2;

supp_inputs.chanspec='0';
supp_inputs.nfft=256;
supp_inputs.overwrite='y';

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

% packaging the inputs into the inputs structure. This can be useful in
% case you need to store the inputs as meta-data with the output. 
inputs=struct('batchname',batchname,'batchpath',batchpath,'spawning_func',mfilename('fullpath'));

% putting file separators at the end of all input paths
if ~isempty(batchpath)
    if ~strcmpi(batchpath(end),filesep)
        batchpath=[batchpath,filesep];
    end
end
%% Body of the function
% loading the batch file
fid=fopen([batchpath batchname],'r');
fpath=batchpath; % assumption
fname=fgetl(fid);

fcfp=1; % fcfp= first chunk first point.
nfft=supp_inputs.nfft;
chanspec=supp_inputs.chanspec;
no_missed_instances=0;
tot_no_instances=0;

while ischar(fname)
    [~,fn,ext]=fileparts(fname);
    simdatafile=[fpath fn '_simdata.mat'];    
    if strcmpi(supp_inputs.overwrite,'n')
        if exist(simdatafile,'file')
            warning('Not writing the *_simdata.mat file for file %s as one already exists',fname);
            fname=fgetl(fid);
            continue
        end
    end
    
    simdata.first_chunk_first_point=fcfp;
    simdata.cbinfile=fname;
    
    
    [dat,fs,~,~]=ReadDataFile([fpath fname],chanspec);        
    
   fcfp= simdata.first_chunk_first_point;
    
    dat=dat(fcfp:end);
    [~,spec,real_spec]=evtafsim(dat,fs,nfft,'',0,0); % the last two inputs decide how the ...
    % spectrogram obtained is normalized and what frequencies are set to zero respectively
    spec=spec'; % each column is a chunk and the rows are different frequencies
    simdata.chunks=spec;
    simdata.specs=real_spec';
    filt_spec=simdata.specs;
    filt_spec(1:6,:)=0;
    simdata.chunk_amp=sum((filt_spec.^2),1);
    
    
    
    no_chunks=length(spec);
    chunk_first_points=fcfp+(0:nfft:nfft*(no_chunks-1));  % short for cfp
    cfp_times=chunk_first_points./fs; 
    chunk_duration=nfft/fs;
    
    simdata.no_chunks=no_chunks;
    simdata.chunk_first_points=chunk_first_points;
    simdata.cfp_times=cfp_times;
    simdata.chunk_duration=chunk_duration;
    simdata.chunk_len=nfft;
    
    % chunk association will come here
    load([fpath fn ext '.not.mat']);    
    
    unlabeled_sylls=strfind(lower(labels),'-');
    unlabeled_sylls2=strfind(lower(labels),'^');
    if length(unlabeled_sylls)==length(labels) || length(unlabeled_sylls2)==length(labels)
       error([[fpath fname] ' has not been labeled.']);        
    end
    
    % converting chunk duration to miliseconds to deal with onsets and
    % offsets
    chunk_duration=chunk_duration*1000;
    
    min_syll_dur=min(offsets-onsets);
    min_gap_dur=min(onsets(2:end)-offsets(1:end-1));
    syll_small=min_syll_dur<=chunk_duration;
    gap_small=min_gap_dur<=chunk_duration;
    
    if syll_small && gap_small 
        error(['Both syll and gaps can be smaller than slice duration for the file ' [fpath fname] ' . Aborting chunk assoc code.'])
   end
    
    
    
    
    simdata.chunk_assocs=cell(1,no_chunks);
    simdata.chunk_assoc_type=cell(1,no_chunks);
    simdata.pc_in_syll=zeros(1,no_chunks);
    simdata.templates=[];
    simdata.simul=[];
    simdata.detect_settings_n_rates=[];
    
    for i=1:no_chunks
        c_on=cfp_times(i)*1000; 
        c_off=c_on+chunk_duration-(chunk_duration/nfft);
        slice_within_syll=find((onsets<=c_on) & (offsets>=c_off));
        slice_within_gap=find((c_off<=onsets(2:end))&(c_on>=offsets(1:end-1)));
        slice_before_first_syll=c_off<=onsets(1);
        slice_after_last_syll=c_on>=offsets(end);
        onset_straddle=find((c_on<onsets)&(c_off>onsets));
        offset_straddle=find((c_on<offsets)&(c_off>offsets));
        syll_within_slice=find((onsets>=c_on)&(offsets<=c_off));
        gap_within_slice=find((offsets(1:end-1)>=c_on)&(onsets(2:end)<=c_off));
        scenarios_realized=0;
        
        if ~isempty(slice_within_syll)
            simdata.chunk_assocs{1,i}=slice_within_syll;
            simdata.chunk_assoc_type{1,i}='within';   
            simdata.pc_in_syll(1,i)=100;
            scenarios_realized=scenarios_realized+1;
        end
        if ~isempty(slice_within_gap)
            simdata.chunk_assocs{1,i}=[slice_within_gap,slice_within_gap+1];
            simdata.chunk_assoc_type{1,i}='between'; 
            simdata.pc_in_syll(1,i)=0;
            scenarios_realized=scenarios_realized+1;
        end
        if slice_before_first_syll
            simdata.chunk_assocs{1,i}=[-1,1];
            simdata.chunk_assoc_type{1,i}='between';  
            simdata.pc_in_syll(1,i)=0;
            scenarios_realized=scenarios_realized+1;
        end
        if slice_after_last_syll
            simdata.chunk_assocs{1,i}=[length(labels),-1];
            simdata.chunk_assoc_type{1,i}='between';   
            simdata.pc_in_syll(1,i)=0;
            scenarios_realized=scenarios_realized+1;
        end
        if ~isempty(onset_straddle) && isempty(offset_straddle) && isempty(syll_within_slice) && isempty(gap_within_slice)
            len_in_syll=c_off-onsets(onset_straddle);
            if len_in_syll/chunk_duration>0.5
                simdata.chunk_assocs{1,i}=onset_straddle;
                simdata.chunk_assoc_type{1,i}='near_onset';       
            else
                if onset_straddle==1
                    simdata.chunk_assocs{1,i}=[-1,onset_straddle];
                    simdata.chunk_assoc_type{1,i}='near_onset'; 
                else
                    simdata.chunk_assocs{1,i}=[onset_straddle-1,onset_straddle];
                    simdata.chunk_assoc_type{1,i}='near_onset'; 
                end
            end  
            simdata.pc_in_syll(1,i)=len_in_syll/chunk_duration*100;
            scenarios_realized=scenarios_realized+1;
        end
        if ~isempty(offset_straddle)&& isempty(onset_straddle)&& isempty(syll_within_slice) && isempty(gap_within_slice)
            len_in_syll=offsets(offset_straddle)-c_on;
            if len_in_syll/chunk_duration>0.5
                simdata.chunk_assocs{1,i}=offset_straddle;
                simdata.chunk_assoc_type{1,i}='near_offset';       
            else
                if offset_straddle==length(labels)
                    simdata.chunk_assocs{1,i}=[offset_straddle,-1];
                    simdata.chunk_assoc_type{1,i}='near_offset'; 
                else
                  simdata.chunk_assocs{1,i}=[offset_straddle,offset_straddle+1];
                  simdata.chunk_assoc_type{1,i}='near_offset'; 
                end
            end
            simdata.pc_in_syll(1,i)=len_in_syll/chunk_duration*100;
            scenarios_realized=scenarios_realized+1;
        end
        if ~isempty(syll_within_slice)
            len_in_syll=offsets(syll_within_slice)-onsets(syll_within_slice);
            simdata.chunk_assocs{1,i}=syll_within_slice;
            simdata.chunk_assoc_type{1,i}='within';
            simdata.pc_in_syll(1,i)=len_in_syll/chunk_duration*100;
            scenarios_realized=scenarios_realized+1;
        end
        if ~isempty(gap_within_slice)
            len_in_syll1=offsets(gap_within_slice)-c_on;
            len_in_syll2=c_off-onsets(gap_within_slice+1);
            len_in_syll=len_in_syll1+len_in_syll2;
            if len_in_syll1/chunk_duration<0.5 && len_in_syll2/chunk_duration<0.5 
                simdata.chunk_assocs{1,i}=[gap_within_slice,gap_within_slice+1];
                simdata.chunk_assoc_type{1,i}='between';       
            elseif len_in_syll1/chunk_duration>0.5
                simdata.chunk_assocs{1,i}=gap_within_slice;
                simdata.chunk_assoc_type{1,i}='near_offset';
            elseif len_in_syll2/chunk_duration>0.5
                simdata.chunk_assocs{1,i}=gap_within_slice+1;
                simdata.chunk_assoc_type{1,i}='near_onset';
            else
                error('Both pre and post syll are more than 50% of the slice. Something is wrong. Call Anand')
            end 
            simdata.pc_in_syll(1,i)=len_in_syll/chunk_duration*100;
            scenarios_realized=scenarios_realized+1;
        end    
        
        if scenarios_realized~=1
            if ~isempty(gap_within_slice) && ~isempty(slice_within_gap)
                if ~(c_on==offsets(gap_within_slice) && c_off==onsets(gap_within_slice+1))
                    error('Multiple/none chunk assignment scenarios have come out true. Something may be wrong. Call Anand')
                end               
            else
                error('Multiple/none chunk assignment scenarios have come out true. Something may be wrong. Call Anand')
            end
        end
        
    end
    
   

    %{
    notmat_file=[batchpath name ext '.not.mat'];
    load(notmat_file)
    inds=strfind(lower(labels),lower(target_syll));
    
    if exist('labels2','var') % this is the way broken or agglomerated syllables are excluded
       broken_sylls=find(labels2~=1);
       corrupt_sylls=[];
       if supp_inputs.exclude_corrupt==1
            corrupt_sylls=find(corrupted==1);
       elseif supp_inputs.exclude_corrupt~=0
           error('The value of supp_inputs.exclude_corrupt should be only 1 or 0')
       end
       inds=setdiff(inds,[broken_sylls,corrupt_sylls]);    
       [parsed_labels,~,~]=parse_complex_label(labels,labels2,labels3,subclips);
       % parse_complex_label uses the complex labeling scheme to arrive at
       % the exact number of different syllables. 
       all_inds=strfind(lower(parsed_labels),lower(srch_string));
       no_missed_instances=no_missed_instances+length(all_inds)-length(inds);
       tot_no_instances=tot_no_instances+length(all_inds);
    else
       % in the absence of complex labeling scheme, broken or agglomerated
       % syllables should be indicated by some special symbol so that they
       % can be easily excluded
      no_missed_instances=no_missed_instances+0; 
      tot_no_instances=tot_no_instances+length(inds);
    end
    
   no_isntances_in_file=length(inds);
   temp_chunks=cell(length(inds),1);
   temp_no_chunks_vec=zeros(length(inds),1);
   temp_durations=zeros(length(inds),1);
    
    for j=1:no_isntances_in_file      
          interval_ind=inds(j);
          assoc_chunks=[];
           for k=1:length(simdata.chunk_assocs)
               if (simdata.chunk_assocs{k})==interval_ind 
                   assoc_chunks=[assoc_chunks,k];
               end
           end 

          temp_chunks{j,l}=simdata.chunks(:,assoc_chunks);
          temp_no_chunks_vec(j,l)=length(assoc_chunks);   
          temp_durations(j,l)=offsets(interval_ind)-onsets(interval_ind);      
    end
   
    % looking for empy temp_chunks
   empty_inds=find(temp_no_chunks_vec==0);
   temp_chunks(empty_inds)=[];
   temp_no_chunks_vec(empty_inds)=[];
   temp_durations(empty_inds)=[];
   
   chunks=[chunks;temp_chunks];
   no_chunks_vec=[no_chunks_vec;temp_no_chunks_vec];
   durations_vec=[durations_vec;temp_durations];
    
    %}
    
    
    
    
    simdata_file=[fpath fn '_simdata.mat'];
    save(simdata_file,'simdata');   
    fname=fgetl(fid);
end








