function []=calculate_optimization_performance_training_set(filename,path,bird_id) 
%% Inputs
% filename - name of the optimization results file
% path - location of the said file

%% Outputs
% It writes a structure called training_set_slice_performance that details the performance of each pre and post-optimization
%  template. 

[~,nm,~]=fileparts(filename);
syll_label=nm(end);

load([path filename])% loads a variable called results
no_templates=length(results);

for i=1:no_templates
    
    training_set_slice_performance.pre(i).pc_targ_slices_missed=results(i).optim_archive.fne_rate_vec(1);
    training_set_slice_performance.pre(i).pc_distractor_slices_hit=results(i).optim_archive.fpe_rate_vec(1);
    training_set_slice_performance.pre(i).threshold=results(i).optim_archive.threshold_vec(1);
    training_set_slice_performance.pre(i).template=results(i).optim_archive.template_vec(:,1);
    training_set_slice_performance.pre(i).sigma=results(i).optim_archive.fin_sigma_vec(1);

    training_set_slice_performance.post(i).pc_targ_slices_missed=results(i).optim_archive.fne_rate_vec(end);
    training_set_slice_performance.post(i).pc_distractor_slices_hit=results(i).optim_archive.fpe_rate_vec(end);
    training_set_slice_performance.post(i).threshold=results(i).optim_archive.threshold_vec(end);
    training_set_slice_performance.post(i).template=results(i).optim_archive.template_vec(:,end);
    training_set_slice_performance.post(i).sigma=results(i).optim_archive.fin_sigma_vec(end);
    
    training_set_slice_performance.bird_id=bird_id;
    training_set_slice_performance.syll_label=syll_label;
    
end

save([path 'training_set_slice_performance_' syll_label '.mat'],'training_set_slice_performance')
 
 
