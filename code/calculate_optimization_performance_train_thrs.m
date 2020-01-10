function []= cv_calculate_optimization_performance_train_thrs (targ_nt,trainpath,testpath,syll_assoc_aggregate,bird_id,chunk_duration)

%% Inputs 
% targ_nt - label fo the target syllable
% trainpath - training data location 
% testpath - test data location
% syll_assoc_aggregate - aggregate of all syllable associated chunks for the tareget data
% bird_id - 
% chunk_duration - duration of individual chunk in ms

%% Output
% It evaluates both average and optimized templates against test data. It
% writes their detection performance in a multichunk_results_* file. The
% data structure of this file is described in detail in  the document
% accompanying this code. 


supp_inputs.normalize=0;
supp_inputs.consec_chunks=[1:5];
supp_inputs.thr_increments=(-1:0.1:1); % corresponding to (0% to 200%)

% if ~isempty(varargin)
% if iscell(varargin{1})
%     if strcmpi(varargin{1}{1},'unpackpvpairs')
%          varargin=varargin{1}(2:end);
%     end
% end
% end
% supp_inputs=parse_pv_pairs(supp_inputs,varargin);

% getting the list of both pre and post templates
optimized_tag='_optimized';
template_ext='.dat';

load([trainpath 'optimization_grand_results_syll_' targ_nt '.mat']) % should load a variable called results
load(['T:\people\anand\optimization_paper_data\' bird_id '\optimized_unconstrained\amp_thr.mat'])% should load variable called amp_thr



% all_templates_struct=dir([trainpath 'template_syll_' targ_nt '_seq_' targ_nt '*' template_ext]);
% post_template_struct=dir([trainpath 'template_syll_' targ_nt '_seq_' targ_nt '*' optimized_tag template_ext]);
% no_post_templates=length(post_template_struct);
% pre_template_list=cell(1,no_post_templates);
% post_template_list=cell(1,no_post_templates);
% 
% for i=1:length(post_template_struct)
%     post_template_list{1,i}=post_template_struct(i).name;
%     [~,nm,~]=fileparts(post_template_struct(i).name);
%     pre_template_fullfile=[trainpath nm(1:end-length(optimized_tag)) template_ext];
%     if exist(pre_template_fullfile,'file')
%         pre_template_list{1,i}=[nm(1:end-length(optimized_tag)) template_ext];        
%     else
%         error('All of the pre and post templates for this syll cannot be found at this location')
%     end    
% end
% 
% all_template_list=cell(1,length(all_templates_struct));
% for i=1:length(all_templates_struct)
%     all_template_list{1,i}=all_templates_struct(i).name;     
% end
% 
% excluded_templates=setdiff(all_template_list,[pre_template_list,post_template_list]);
% optimization_skipped_templates=[trainpath 'optimization_skipped_templates_' targ_nt '.mat'];
% 
% if exist(optimization_skipped_templates,'file')
%     load(optimization_skipped_templates)% loads variable called skipped_templates and no_all_templates
%     if ~isequal(sort(excluded_templates),sort(skipped_templates))&& ~(isempty(excluded_templates)&& isempty(skipped_templates))
%        error('Excluded templates and skipped templates are not the same. Something may be wrong') 
%     end
%     if no_all_templates-length(excluded_templates)~=length(post_template_list)
%        error('There is a mismatch in the number of templates found after ignoring skipped templates') 
%     end        
% else
%    error('The optimization skipped templates for this target note does not exist') 
% end

no_templates=length(results);

% if length(results)~=no_templates    
%    error('The number of templates in results should be equalt to the number of templates') 
% end

% [pre_template_list]=sort_template_filenames(pre_template_list);
% [post_template_list]=sort_template_filenames(post_template_list);
% for ii=1:no_templates
%     prename=
%     template_list=[pre_template_list;post_template_list];
% 
% end
% Getting metadata for the templates
% We assume that the metadata for all the slices/templates are the same
% [~,nm,~]=fileparts(pre_template_list{1,1});
% metadata_file=[trainpath nm '_metadata.mat'];
% if ~exist(metadata_file,'file')
%     error('The metadata file does not exist at the specified location');
% end
% load(metadata_file); % loads a variable called template_metadata
% if template_metadata.length_in_chunks_of_target_instances-length(skipped_templates)~=length(pre_template_list)
%    error('These is a mismatch between the metadata and the actual number of template') 
% end

% Ontaining target_ind
% distractor_sylls=template_metadata.distractor_sylls;
sylls_in_aggregate=[syll_assoc_aggregate.syll];
% if ~isequal(sort([distractor_sylls,targ_nt]),sort(sylls_in_aggregate))
%     if  ~isequal(sort([distractor_sylls,targ_nt]),strrep(sort(sylls_in_aggregate),'*',''))
%         error('The set of syllables obtainesd from the template metadata and from the aggregate are different')
%     end
% end
target_ind=strfind(sylls_in_aggregate,targ_nt);



% now that we have syll_assoc_aggregate, lets obtain the optimal threshold
% for each template
% load('temp_vars_for_easier_debug.mat')
templates=cell(2,no_templates);
optimal_thresholds=zeros(2,no_templates);
template_list=cell(2,no_templates);

for i=1:no_templates    
    [pth,nm,xt]=fileparts(results(i).templatefullfile);
    pth=[pth filesep];
    for j=1:2
        if j==1
            curr_template=load([pth nm xt]);
             threshold=results(i).optim_archive.threshold_vec(1);   
             template_list{j,i}=[nm xt];
        else
            curr_template=load([pth nm optimized_tag xt]);
            threshold=results(i).optim_archive.threshold_vec(end);
            template_list{j,i}=[nm optimized_tag xt];
        end
        if supp_inputs.normalize==1
            curr_template=curr_template-min(curr_template);
            curr_template=curr_template./max(curr_template);
        elseif supp_inputs.normalize~=0
            error('The value of normalize can be only 1 or 0')            
        end
                
       templates{j,i}=curr_template;
       optimal_thresholds(j,i)=threshold^2; % squaring the threshold          
    end  
end

for i=1:no_templates    
    for j=1:2
        curr_template_file=template_list{j,i};
        curr_template=templates{j,i};
        [dist_vecs_aggregate,amplitudes_aggregate,onsets_aggregate]=calculate_dist_vecs(syll_assoc_aggregate,curr_template);
        for k=1:length(supp_inputs.thr_increments)
            curr_thr_increment=supp_inputs.thr_increments(k);
            curr_thr=optimal_thresholds(j,i)+curr_thr_increment*optimal_thresholds(j,i);
            cnt_vecs_aggregate=calculate_cnt_vecs(dist_vecs_aggregate,curr_thr,amplitudes_aggregate,amp_thr); 
            longest_streak_aggregate=calculate_longest_streaks_cnt_vecs(cnt_vecs_aggregate,amplitudes_aggregate,onsets_aggregate,syll_assoc_aggregate,'collapse',1); %% WRITE THIS FUNC
            target_streaks=[];
            distractor_streaks=[];
            trg_onsets=[];
            trg_slices=[];
            for m=1:length(sylls_in_aggregate)
                if m==target_ind
                   target_streaks=longest_streak_aggregate(1,m).streaks; 
                   trg_onsets=longest_streak_aggregate(1,m).onsets;
                   trg_slices=longest_streak_aggregate(1,m).slices;
                else
                   distractor_streaks=[distractor_streaks;longest_streak_aggregate(1,m).streaks];
                end                    
            end            
            for l=1:length(supp_inputs.consec_chunks)
                curr_consec_slice=supp_inputs.consec_chunks(l);
                hiinds=find(target_streaks>=curr_consec_slice);
                true_positive_onsets=trg_onsets(hiinds);
                if isempty(true_positive_onsets)
                    trg_latencies=[];
                else
                    trg_latencies=zeros(size(true_positive_onsets,1),1);
                    for po=1:size(true_positive_onsets,1)
                        trg_latencies(po)=true_positive_onsets(po)+(trg_slices{hiinds(po),1}(curr_consec_slice)-1)*chunk_duration*1000;
                    end
                end
                no_true_positive_trigs=length(find(target_streaks>=curr_consec_slice));
                no_false_negative_trigs=length(find(target_streaks<curr_consec_slice));
                no_true_negative_trigs=length(find(distractor_streaks<curr_consec_slice));
                no_false_positive_trigs=length(find(distractor_streaks>=curr_consec_slice));
                
                % archiving results
                if j==1
                    multichunk_results.pre(i,l,k).template_file=curr_template_file;
                    multichunk_results.pre(i,l,k).template=curr_template;
                    multichunk_results.pre(i,l,k).threshold=curr_thr;
                    multichunk_results.pre(i,l,k).threshold_incre=curr_thr_increment;
                    multichunk_results.pre(i,l,k).consec_chunk=curr_consec_slice;

                    multichunk_results.pre(i,l,k).no_true_positive_trigs=no_true_positive_trigs;
                    multichunk_results.pre(i,l,k).no_false_negative_trigs=no_false_negative_trigs;
                    multichunk_results.pre(i,l,k).no_true_negative_trigs=no_true_negative_trigs;
                    multichunk_results.pre(i,l,k).no_false_positive_trigs=no_false_positive_trigs; 
                    multichunk_results.pre(i,l,k).trg_jitter=std(trg_latencies); 
                else
                    multichunk_results.post(i,l,k).template_file=curr_template_file;
                    multichunk_results.post(i,l,k).template=curr_template;
                    multichunk_results.post(i,l,k).threshold=curr_thr;
                    multichunk_results.post(i,l,k).threshold_incre=curr_thr_increment;
                    multichunk_results.post(i,l,k).consec_chunk=curr_consec_slice;

                    multichunk_results.post(i,l,k).no_true_positive_trigs=no_true_positive_trigs;
                    multichunk_results.post(i,l,k).no_false_negative_trigs=no_false_negative_trigs;
                    multichunk_results.post(i,l,k).no_true_negative_trigs=no_true_negative_trigs;
                    multichunk_results.post(i,l,k).no_false_positive_trigs=no_false_positive_trigs;
                    multichunk_results.post(i,l,k).trg_jitter=std(trg_latencies); 
                    
                end      
                
            end
            
        end      
        
    end   
end
%%

for i=1:no_templates        
        for k=1:length(supp_inputs.consec_chunks)       
            for l=1:length(supp_inputs.thr_increments)
                multichunk_results.pre(i,k,l).pc_targets_missed=multichunk_results.pre(i,k,l).no_false_negative_trigs/...
                                           (multichunk_results.pre(i,k,l).no_false_negative_trigs+multichunk_results.pre(i,k,l).no_true_positive_trigs)*100;   

                multichunk_results.pre(i,k,l).pc_distractors_hit=multichunk_results.pre(i,k,l).no_false_positive_trigs/...
                                            (multichunk_results.pre(i,k,l).no_false_positive_trigs+multichunk_results.pre(i,k,l).no_true_negative_trigs)*100;    
                                                       
                multichunk_results.pre(i,k,l).pc_distractors_hit_wrt_targets=multichunk_results.pre(i,k,l).no_false_positive_trigs/...
                                           (multichunk_results.pre(i,k,l).no_false_negative_trigs+multichunk_results.pre(i,k,l).no_true_positive_trigs)*100;  
                                       
                                               
                multichunk_results.post(i,k,l).pc_targets_missed=multichunk_results.post(i,k,l).no_false_negative_trigs/...
                                           (multichunk_results.post(i,k,l).no_false_negative_trigs+multichunk_results.post(i,k,l).no_true_positive_trigs)*100;   
                multichunk_results.post(i,k,l).pc_distractors_hit=multichunk_results.post(i,k,l).no_false_positive_trigs/...
                    (multichunk_results.post(i,k,l).no_false_positive_trigs+multichunk_results.post(i,k,l).no_true_negative_trigs)*100;                 
                                        
                multichunk_results.post(i,k,l).pc_distractors_hit_wrt_targets=multichunk_results.post(i,k,l).no_false_positive_trigs/...
                                           (multichunk_results.post(i,k,l).no_false_negative_trigs+multichunk_results.post(i,k,l).no_true_positive_trigs)*100; 
                                   
                 
            end           
        end
end
    
   
     err_rates2_pre=([multichunk_results.pre(:).pc_targets_missed]+[multichunk_results.pre(:).pc_distractors_hit_wrt_targets])/2;
     allthrs=[multichunk_results.pre(:).threshold];
     pc_targets_missed_pre=[multichunk_results.pre(:).pc_targets_missed];
     pc_distractors_hit_wrt_targets=[multichunk_results.pre(:).pc_distractors_hit_wrt_targets];  
     all_jitter=[multichunk_results.pre(:).trg_jitter];
    [err_rates2_pre,ord]=sort(err_rates2_pre);
    pc_targets_missed_pre=pc_targets_missed_pre(ord);
    pc_distractors_hit_wrt_targets=pc_distractors_hit_wrt_targets(ord);   
    allthrs=allthrs(ord);
    all_jitter=all_jitter(ord);
    [ind1,ind2,~]=ind2sub(size(multichunk_results.pre),ord);

    descriptors={'Balanced error','Template slice #','# of consec. slices','threshold','% targets missed','% distr. hit wrt targets','targ jitter'};
    ttemp=[err_rates2_pre',ind1',supp_inputs.consec_chunks(ind2')',allthrs',pc_targets_missed_pre',pc_distractors_hit_wrt_targets',all_jitter'];
    ttemp=mat2cell(ttemp,ones(size(ttemp,1),1),ones(size(ttemp,2),1));
    ttemp=[descriptors;ttemp];
    multichunk_results.pre_sorted=ttemp;
    
     err_rates2_post=[multichunk_results.post(:).pc_targets_missed]+[multichunk_results.post(:).pc_distractors_hit_wrt_targets]/2;
     allthrs=[multichunk_results.post(:).threshold];
     pc_targets_missed_post=[multichunk_results.post(:).pc_targets_missed];
     pc_distractors_hit_wrt_targets=[multichunk_results.post(:).pc_distractors_hit_wrt_targets];    
     all_jitter=[multichunk_results.post(:).trg_jitter];
    [err_rates2_post,ord]=sort(err_rates2_post);
    pc_targets_missed_post=pc_targets_missed_post(ord);
    pc_distractors_hit_wrt_targets=pc_distractors_hit_wrt_targets(ord);   
    allthrs=allthrs(ord);
    all_jitter=all_jitter(ord);
    [ind1,ind2,~]=ind2sub(size(multichunk_results.post),ord);

    descriptors={'Balanced error','Template slice #','# of consec. slices','threshold','% targets missed','% distr. hit wrt targets','targ jitter'};
    ttemp=[err_rates2_post',ind1',supp_inputs.consec_chunks(ind2')',allthrs',pc_targets_missed_post',pc_distractors_hit_wrt_targets',all_jitter'];
    ttemp=mat2cell(ttemp,ones(size(ttemp,1),1),ones(size(ttemp,2),1));
    ttemp=[descriptors;ttemp];
    multichunk_results.post_sorted=ttemp;



%%


multichunk_results.tot_no_instances=syll_assoc_aggregate(1,target_ind).tot_no_instances;
multichunk_results.no_missed_instances_mis_seg=syll_assoc_aggregate(1,target_ind).no_missed_instances;
multichunk_results.bird_id=bird_id;
multichunk_results.consec_chunks=supp_inputs.consec_chunks;
multichunk_results.thr_increments=supp_inputs.thr_increments;

save([testpath 'op_multichunk_results_' targ_nt '.mat'],'multichunk_results')
 
 
