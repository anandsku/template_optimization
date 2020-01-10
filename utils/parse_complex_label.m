function [labels,orig_inds,new_inds_real,varargout]=parse_complex_label(labels1,labels2,labels3,subclips,varargin)
%% Syntax
% It is assumed that a file with labels2 and labels 3 also has subclips
% [outputs]=function_template(inp1,inp2,inp3,inp4,varargin)
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
narg_min=4;

supp_inputs.onsets=[];
supp_inputs.offsets=[];

if nargin<narg_min
    error(['The number of inputs should at least be ' narg_min])
else
   supp_inputs=parse_pv_pairs(supp_inputs,varargin); 
end


%% Body of the function
if ~(length(labels1)==length(labels2) && length(labels1)== length(labels3))
    error('The length of the three labels have to be equal');    
end 

if ~iscell(labels3)
   error('Labels 3 has to be a cell array')    
end
ascii_nums=(48:57);


orig_inds=(1:length(labels1));
new_inds=cell(0);

curr_ind=1;
while curr_ind<=length(labels1)
    if isequal(labels1(curr_ind),'=')
        no_subclips=length(subclips{curr_ind});
        temp_labels1=[];
        temp_labels2=[];  
        temp_subclips=cell(0);      
        for i=1:length(subclips{curr_ind})
            curr_subclip=subclips{1,curr_ind}{1,i};
            if ~ismember(double(curr_subclip(1)),ascii_nums)
                temp_labels1=[temp_labels1,curr_subclip(1)];
                temp_subclips=[temp_subclips,{curr_subclip(1)}];
            else
                error('The first element in the subclip string cannot be a number');
            end
            if length(curr_subclip)==1
                temp_labels2=[temp_labels2,1];
            elseif length(curr_subclip)==2
                if ismember(double(curr_subclip(2)),ascii_nums)
                    temp_labels2=[temp_labels2,str2num(curr_subclip(2))];
                else
                    error('The second element in the subclip string has to be a number');
                end
            else
               error('The size of a subclip string cannot exceed 2') 
            end
        end
        temp_labels3=repmat({''},1,length(temp_labels1));
        
        labels1=[labels1(1:curr_ind-1),temp_labels1,labels1(curr_ind+1:end)];
        labels2=[labels2(1:curr_ind-1),temp_labels2,labels2(curr_ind+1:end)];
        labels3=[labels3(1:curr_ind-1),temp_labels3,labels3(curr_ind+1:end)];
        subclips=[subclips(1:curr_ind-1),temp_subclips,subclips(curr_ind+1:end)];
        new_inds=[new_inds,{curr_ind:curr_ind+no_subclips-1}];
        curr_ind=curr_ind+no_subclips;
    else
        new_inds=[new_inds,{curr_ind}];
        curr_ind=curr_ind+1;
    end  
    
end

labels=[];
prev_broken=0;
newer_inds=zeros(1,length(labels1));
for ctr=1:length(labels1)       
    if isequal(labels1(ctr),'=')
        error('Compound labels should have been dealt with earlier in this code')  
    elseif labels2(ctr)~=1
        if ~prev_broken
            labels=[labels,labels1(ctr)];
            newer_inds(ctr)=length(labels);
            prev_broken=1;
        else
            if ~(labels1(ctr)==labels1(ctr-1) && labels2(ctr)-labels2(ctr-1)==1)
                labels=[labels,labels1(ctr)];
               newer_inds(ctr)=length(labels);
                prev_broken=1;
            else
                newer_inds(ctr)=length(labels);
            end
        end
    else
       labels=[labels,labels1(ctr)];
       newer_inds(ctr)=length(labels);
       prev_broken=0;  
    end    
end

new_inds_real=cell(1,length(orig_inds));
for i=1:length(orig_inds)
    for j=1:length(new_inds{i})
        new_inds_real{1,i}=[new_inds_real{1,i},newer_inds(new_inds{i}(j))];      
    end   
end

if ~isempty(supp_inputs.onsets) && ~isempty(supp_inputs.offsets)   

    % now for determining the new_onsets and the new_offsets
    onsets_new=zeros(1,length(labels))-1;
    offsets_new=zeros(1,length(labels))-1;

    for i=1:length(labels)
       containing_orig_inds=[];
       pos=[];
       for j=1:length(new_inds_real)
           if ismember(i,new_inds_real{j})
               containing_orig_inds=[containing_orig_inds,j];
               if new_inds_real{j}(1)==i && new_inds_real{j}(end)==i
                   pos=[pos,2];
               elseif new_inds_real{j}(1)==i 
                   pos=[pos,0];
               elseif new_inds_real{j}(end)==i
                   pos=[pos,1];
               else
                   pos=[pos,-1];
               end
           end

       end
       if length(containing_orig_inds)>2
           error('A syllable split into more than two parts cannot be handles by this system')
       end
       % determining onset
       if (pos(1)==2 || pos(1)==0) && labels2(new_inds{containing_orig_inds(1)}(1))~=3 % the first part should not have 3 as label2
           onsets_new(1,i)=supp_inputs.onsets(containing_orig_inds(1));
       end  

       % determining offset
       if (pos(end)==2 || pos(end)==1) && labels2(new_inds{containing_orig_inds(end)}(end))~=2 % the last part should not have 2 as label2
          offsets_new(1,i)=supp_inputs.offsets(containing_orig_inds(end)); 
       end

    end
    varargout{1}=onsets_new;
    varargout{2}=offsets_new;
else
   varargout{1}=[];
   varargout{2}=[];
    
end



