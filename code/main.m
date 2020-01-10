%%
clear
close all

%% inputs start
birdID='340';
allsylls='ACDEFGHTWZ';
targsylls='ACDEFGHT';
root='T:\people\anand\optimization_paper_data\upload_data\resub_supp\';

stem='analysis_data\340\';
train_dir_name='train_dir';
test_dir_name='test_dir';

train_batch_name='batch_train';
test_batch_name='batch_test';

trainpath=[stem train_dir_name filesep];
testpath=[stem test_dir_name filesep];
%% inputs end

%create simdata file
create_simdata_file(train_batch_name,trainpath)
create_simdata_file(test_batch_name,testpath)

% calculate amplitude threshold
calculate_amplitude_threshold(train_batch_name,trainpath)


% create syll associated chunks files
for ii=1:length(allsylls)
    target_syll=allsylls(ii);
     make_syll_assoc_chunks(train_batch_name,trainpath,target_syll);
     make_syll_assoc_chunks(test_batch_name,testpath,target_syll);
end


% create gap assoc chunks
load([trainpath 'amp_thr.mat'])
make_gap_assoc_high_amp_chunks(train_batch_name,trainpath,amp_thr);
make_gap_assoc_chunks_keep_together(test_batch_name,testpath)


% align syll assoc chunks
for ii=1:length(targsylls)
    filename=['syll_assoc_chunks_syll_' targsylls(ii) '_seq_' targsylls(ii) '.mat'];
    align_syll_assoc_chunks(filename,trainpath);
end


% make templates
for ii=1:length(targsylls)
    target_syll=targsylls(ii);
    distractor_sylls=strrep(allsylls,target_syll,'');
    make_templates(target_syll,distractor_sylls,trainpath,birdID)
    
end


% optimize templates 

for ii=1:length(targsylls)
    currsyll=upper(targsylls(ii));
    syllassocfile=[trainpath 'syll_assoc_chunks_syll_' currsyll '_seq_' currsyll];
    load(syllassocfile)   
    no_templates=size(syll_assoc_chunks{1,1}{1,1}{1,1},2);
    grand_results=[];
    skipped_templates=cell(0);
    yyy=waitbar(0,['done with ' num2str(0) ' templates outof ' num2str(no_templates) ' for syll ' targsylls(ii)]);
    cum_time=0;
    for jj=1:no_templates
        template_name=['template_syll_' currsyll '_seq_' currsyll '_chunks_' num2str(jj) '_outof_' num2str(no_templates) '.dat'];      
        tic
        try
            [results]=optimize_template(template_name,trainpath);
                     if  results.optim_archive.force_stopped==1
                         disp('force stopped')
                     end
        catch err_exp
            if isequal(err_exp.message,'The sigma supplied is not sufficient for constructing smooth, monotonous density distributions')
                warning(['Skipping ' template_name ' due to low init sigma'])
                skipped_templates=[skipped_templates,template_name];
                continue
            elseif isequal(err_exp.message,'The shape of the distributions is not right') 
                warning(['Skipping ' template_name ' due to irregularly shaped distributions'])
                skipped_templates=[skipped_templates,template_name];
                continue
            else
                throw(err_exp)
            end
            
        end    
        

        grand_results=[grand_results,results];
        
        e_time=toc/60;
        cum_time=cum_time+e_time;
        mins_str=sprintf('%0.3f',cum_time);
        waitbar(jj/no_templates,yyy,['done with ' num2str(jj) ' templates outof ' num2str(no_templates) ' in ' mins_str ' mins for syll ' targsylls(ii)]);   
        
    end 
    
    results=grand_results;
    fullfilename=[trainpath 'optimization_grand_results_syll_' currsyll '.mat'];
    save(fullfilename,'results')
    save([trainpath 'optimization_skipped_templates_' currsyll '.mat'],'skipped_templates','no_templates');
end


% write optimized templates

for ii=1:length(targsylls)    
    res_filename=['optimization_grand_results_syll_' targsylls(ii) '.mat'];
    write_optimized_templates_from_optimization_results(res_filename,trainpath)
end


% calculate training set performance

for ii=1:length(targsylls)    
    res_filename=['optimization_grand_results_syll_' targsylls(ii) '.mat'];
    calculate_optimization_performance_training_set(res_filename,trainpath,birdID) 
end






% evaluating 

    syll_assoc_aggregate=[];
    examplesimdatafile=dir([testpath '*_simdata.mat']);
    simstruct=load([examplesimdatafile(1).folder filesep examplesimdatafile(1).name]);
    chunk_duration=simstruct.simdata.chunk_duration;


    for ii=1:length(allsylls)

        load([testpath 'syll_assoc_chunks_syll_' allsylls(ii) '_seq_' allsylls(ii) '.mat']);% loads a variable called syll_assoc_chunks

        syll_assoc_aggregate(ii).chunks=syll_assoc_chunks;
        syll_assoc_aggregate(ii).syll=allsylls(ii);
        syll_assoc_aggregate(ii).durations=syll_assoc_durations;
        syll_assoc_aggregate(ii).amps=syll_assoc_amplitudes;
        syll_assoc_aggregate(ii).onsets=syll_assoc_onsets;
        syll_assoc_aggregate(ii).no_missed_instances=no_missed_instances;
        syll_assoc_aggregate(ii).tot_no_instances=tot_no_instances; 


    end
    % loading gaps assoc chunks
    load([testpath 'gap_assoc_chunks_kept_together.mat']);
    syll_assoc_aggregate(ii+1).chunks=gap_assoc_chunks;
    syll_assoc_aggregate(ii+1).syll='*';
    syll_assoc_aggregate(ii+1).durations=[];
    syll_assoc_aggregate(ii+1).no_missed_instances=0;
    syll_assoc_aggregate(ii+1).tot_no_instances=Inf; 
    syll_assoc_aggregate(ii+1).amps=gap_assoc_amplitudes;

    
    for jj=1:length(targsylls)
        targ_nt=targsylls(jj);
        calculate_optimization_performance_train_thrs(targ_nt,trainpath,testpath,syll_assoc_aggregate,birdID,chunk_duration)
    end



