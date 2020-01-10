function [longest_streak_aggregate]=calculate_longest_streaks_cnt_vecs(cnt_vecs_aggregate,amplitudes_aggregate,onsets_aggregate,syll_assoc_aggregate,varargin)
%% Syntax
%
% [longest_streak_aggregate]=calculate_longest_streaks_cnt_vecs(cnt_vecs_aggregate,varargin)
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
narg_min=1;


% % Assigning default values to supplementary inputs
supp_inputs.collapse=0;

if nargin<narg_min
    error(['The number of inputs should at least be ' narg_min])  
else
    % processing supplementary inputs
    supp_inputs=parse_pv_pairs(supp_inputs,varargin);
end

%% Body of the function
no_sylls_in_aggregate=length(cnt_vecs_aggregate);
longest_streak_aggregate(1,no_sylls_in_aggregate).streaks=cell(0);
longest_streak_aggregate(1,no_sylls_in_aggregate).amps=cell(0); 
longest_streak_aggregate(1,no_sylls_in_aggregate).onsets=cell(0); 
longest_streak_aggregate(1,no_sylls_in_aggregate).slices=cell(0); 
longest_streak_aggregate(1,no_sylls_in_aggregate).durs=cell(0); 



for i=1:no_sylls_in_aggregate
    no_intervals=length(cnt_vecs_aggregate(1,i).cnts);
    longest_streak_aggregate(1,i).streaks=cell(1,no_intervals);
    longest_streak_aggregate(1,i).amps=cell(1,no_intervals);
    longest_streak_aggregate(1,i).onsets=cell(1,no_intervals);
    longest_streak_aggregate(1,i).slices=cell(1,no_intervals);
    longest_streak_aggregate(1,i).durs=cell(1,no_intervals);
     
    for j=1:no_intervals
       no_lengths=length(cnt_vecs_aggregate(1,i).cnts{1,j});
       longest_streak_aggregate(1,i).streaks{1,j}=cell(1,no_lengths);
       longest_streak_aggregate(1,i).amps{1,j}=cell(1,no_lengths);
       longest_streak_aggregate(1,i).onsets{1,j}=cell(1,no_lengths);
       longest_streak_aggregate(1,i).slices{1,j}=cell(1,no_lengths);
       longest_streak_aggregate(1,i).durs{1,j}=cell(1,no_lengths);
       
       for k=1:no_lengths
           no_instances=size(cnt_vecs_aggregate(1,i).cnts{1,j}{1,k},1);
           longest_streak_aggregate(1,i).streaks{1,j}{1,k}=zeros(no_instances,1);
           longest_streak_aggregate(1,i).amps{1,j}{1,k}=cell(no_instances,1);
           longest_streak_aggregate(1,i).onsets{1,j}{1,k}=zeros(no_instances,1);
           longest_streak_aggregate(1,i).slices{1,j}{1,k}=cell(no_instances,1);
           longest_streak_aggregate(1,i).durs{1,j}{1,k}=zeros(no_instances,1);
           for l=1:no_instances
               temp_cnts=cnt_vecs_aggregate(1,i).cnts{1,j}{1,k}(l,:);
               [temp_streak,streakslices]=calculate_longest_streak(temp_cnts);
               longest_streak_aggregate(1,i).streaks{1,j}{1,k}(l,1)=temp_streak; 
               longest_streak_aggregate(1,i).amps{1,j}{1,k}{l,:}=amplitudes_aggregate(1,i).amps{1,j}{1,k}(l,streakslices); 
               longest_streak_aggregate(1,i).onsets{1,j}{1,k}(l,1)=onsets_aggregate(1,i).onsets{1,j}{1,k}(l,1);
               longest_streak_aggregate(1,i).slices{1,j}{1,k}{l,:}=streakslices; 
               if ~isempty(syll_assoc_aggregate(1,i).durations)
                    longest_streak_aggregate(1,i).durs{1,j}{1,k}(l,:)=syll_assoc_aggregate(1,i).durations{1,j}{1,k}(l,:);               
               end
            end      
        end   
    end
end

% WE HAVE TAKEN CARE OF amplitudes in THE CALCULATION OF CNT_VEC. SO WHEN
% WE ARE COLLAPSING THIS, WE ARE NOT INCLUDING AMPLITUDES. 
if supp_inputs.collapse==1
    longest_streak_aggregate_temp(1,no_sylls_in_aggregate).streaks=[];
    longest_streak_aggregate_temp(1,no_sylls_in_aggregate).onsets=[];
    longest_streak_aggregate_temp(1,no_sylls_in_aggregate).slices=[];
    longest_streak_aggregate_temp(1,no_sylls_in_aggregate).durs=[];
    for i=1:no_sylls_in_aggregate
        no_intervals=length(cnt_vecs_aggregate(1,i).cnts);
        for j=1:no_intervals
           no_lengths=length(cnt_vecs_aggregate(1,i).cnts{1,j});
           for k=1:no_lengths
               longest_streak_aggregate_temp(1,i).streaks=...
               [longest_streak_aggregate_temp(1,i).streaks;longest_streak_aggregate(1,i).streaks{1,j}{1,k}];    
                longest_streak_aggregate_temp(1,i).onsets=...
               [longest_streak_aggregate_temp(1,i).onsets;longest_streak_aggregate(1,i).onsets{1,j}{1,k}];
               longest_streak_aggregate_temp(1,i).slices=...
               [longest_streak_aggregate_temp(1,i).slices;longest_streak_aggregate(1,i).slices{1,j}{1,k}]; 
           
           
               longest_streak_aggregate_temp(1,i).durs=...
               [longest_streak_aggregate_temp(1,i).durs;longest_streak_aggregate(1,i).durs{1,j}{1,k}];    
    
           end      
        end   
    end
    longest_streak_aggregate=longest_streak_aggregate_temp;
end
