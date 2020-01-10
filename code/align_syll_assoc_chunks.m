function [syll_assoc_chunks_aligned]=align_syll_assoc_chunks(filename,filepath,varargin)
%% Syntax
%
% []=align_syll_assoc_chunks(filename,filepath,file_fate,varargin)
%
%% Inputs  
% filename - name of the syll_assoc_chunks filethat needs to be aligned
% filepath - location of the said file
%
%
%
%% Computation/Processing     % 
% It makes all syllable instances (except those outside 2 standard deviations away from the mean duration)
% to have the same length in terms of number of slices. It does this by linearly stretching or compressing slices
% that are not of modal slice length. 
%
% 
%
%% Outputs  
% 'syll_assoc_chunks' - cell array containing the aligned syllable
% instances
% 'syll_assoc_amplitudes' - amplitudes for each slice within the syllable
% instances
% 'durations' - durations of all syllable instances
% 'no_missed_instances'- number of instances discarded due to their length being 
% outside the 2 standard deviations of mean
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



% Assigning default values to supplementary inputs
supp_inputs.exclusion_criterion=2; % in terms of # of SDs away from the mean
supp_inputs.exclusion_threshold=0.1; % in terns of the max proportion of data allowed to be thrown away under the above criterion
supp_inputs.ask_user=1;

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
if ~isempty(filepath)
    if ~strcmpi(filepath(end),filesep)
        filepath=[filepath,filesep];
    end
end

%% Body of the function
syll_assoc_chunks_aligned=-1;

load([filepath filename]);% should load various variables. importantly syll_assoc_chunks and syll_assoc_durations and no_missed_instances
    
if exist('durations','var') 
   error('It seems the chunks have already been aligned') 
end



syll_assoc_chunks_new=cell(1,1);
amps_new=cell(1,1);

no_lengths=length(syll_assoc_chunks{1,1}); 
% gathering the distribution of syll lengths
durations=[];
lens=zeros(1,no_lengths);
ins=zeros(1,no_lengths);
for j=1:no_lengths
   durations=[durations;syll_assoc_durations{1,1}{1,j}];      
   lens(1,j)=size(syll_assoc_chunks{1,1}{1,j}{1,1},2);
   ins(1,j)=size(syll_assoc_chunks{1,1}{1,j},1);
end
up_lim=mean(durations)+supp_inputs.exclusion_criterion*std(durations); 
lo_lim=mean(durations)-supp_inputs.exclusion_criterion*std(durations); 

exclude_inds=find((durations>=up_lim)|(durations<=lo_lim));
exclude_prop=length(exclude_inds)/length(durations);

if exclude_prop>supp_inputs.exclusion_threshold
   if supp_inputs.ask_user==1
        quest_msg=['You had chosen ' num2str(supp_inputs.exclusion_threshold)...
                  ' as the exclusion threshold. The proportion of instances '...
                  'that are being excluded for this current file are ' num2str(exclude_prop)...
                  '. Do you want continue anyway or abort the alignment for this file?'];
        quest_dialog_title='';
        quest_option1='Continue anyway';
        quest_option2='Abort';
        quest_default_option=quest_option2;

        quest_resp=questdlg(quest_msg,quest_dialog_title,quest_option1,quest_option2,quest_default_option);
        if strcmpi(quest_resp,quest_option2)||isempty(quest_resp)
           disp('Aborting the alignment since more than desired # of instacnes were being thrown away')
           return      
        end       
   elseif supp_inputs.ask_user==0
      warning(['You had chosen ' num2str(supp_inputs.exclusion_threshold)...
                  ' as the exclusion threshold. The proportion of instances '...
                  'that are being excluded for this current file are ' num2str(exclude_prop)]) 
   else
       error('Ask user can only have values 0 or 1')
   end
end
syll_assoc_chunks_new{1,1}{1,1}=cell(length(durations)-length(exclude_inds),1);
amps_new{1,1}{1,1}=[];

no_missed_instances=no_missed_instances+length(exclude_inds);

[~,modal_ind]=max(ins);
modal_len=lens(modal_ind);


cnt=0;   
for j=1:no_lengths
   curr_len=lens(1,j);
   if curr_len==0
       error(' Number of chunks should not be zero. Please recalculate simdata files with pretime=0')
   end


   chunk_contri=zeros(curr_len,modal_len);

   unit_size=curr_len/modal_len;
   chunk_times=(0:unit_size:curr_len);
   chunk_starts=chunk_times(1:end-1);
   chunk_ends=chunk_times(2:end);        
   chunk_starts_old=(0:1:curr_len-1);
   chunk_ends_old=(1:1:curr_len);

   % calculating chunk_contri
   for k=1:modal_len
        chunks_present_temp1=find((chunk_starts_old>=chunk_starts(1,k))&(chunk_starts_old<chunk_ends(1,k)));
        % is there a chunk that starts between these guys
        chunks_present_temp2=find((chunk_ends_old<=chunk_ends(1,k))&(chunk_ends_old>chunk_starts(1,k)));
        % is there a chunk that ends between these guys
        chunks_present=union(chunks_present_temp1,chunks_present_temp2);

        if isempty(chunks_present)% when the chunk old is entirely outsie the limits of the the new chunks
            ch1=find(chunk_starts_old<=chunk_starts(1,k));
            ch1=ch1(end);
            ch2=find(chunk_ends_old>=chunk_ends(1,k));
            ch2=ch2(1);
            if ch2==ch1
                chunks_present=ch1;
            else
                error('No chunks present in this one')
            end     

        end

        for l=1:length(chunks_present)
            chunk_ind=chunks_present(1,l);
           left_diff=chunk_starts_old(chunk_ind)-chunk_starts(1,k);
           right_diff=chunk_ends(1,k)-chunk_ends_old(chunk_ind);
           if left_diff<0
               left_diff=0;
           end

           if right_diff<0
               right_diff=0;
           end
           total_diff=left_diff+right_diff;
           chunk_contri(chunk_ind,k)=(unit_size-total_diff)/unit_size;              

        end           
    end


   no_instances_in_length=ins(1,j);

   for k=1:no_instances_in_length
       dur=syll_assoc_durations{1,1}{1,j}(k,1);
       if dur<=lo_lim || dur >=up_lim
            continue
       else
           cnt=cnt+1;
       end
       curr_spec=syll_assoc_chunks{1,1}{1,j}{k,1};
       curr_amp=syll_assoc_amplitudes{1,1}{1,j}(k,:);
       transformed_spec=curr_spec*chunk_contri;
       transformed_amp=curr_amp*chunk_contri;

       for o=1:size(transformed_spec,2)
          transformed_spec(:,o)=transformed_spec(:,o)-min(transformed_spec(:,o)); 
          transformed_spec(:,o)=transformed_spec(:,o)./max(transformed_spec(:,o));
       end  

%            if ~isempty(find(isnan(transformed_spec(1,:)),1))
%                disp('boo');
%            end
       syll_assoc_chunks_new{1,1}{1,1}{cnt,1}=transformed_spec;
       amps_new{1,1}{1,1}(cnt,:)=transformed_amp;
   end


end   
   
    


syll_assoc_chunks_aligned=syll_assoc_chunks_new;
syll_assoc_amplitudes=amps_new;
syll_assoc_chunks=syll_assoc_chunks_new;


%% Processing outputs 
    
matfullfile=[filepath filename];
save(matfullfile,'syll_assoc_chunks','syll_assoc_amplitudes','durations','no_missed_instances','-append'); % the durations store all durations and not just the ones aligned
  

