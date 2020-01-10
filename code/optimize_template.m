function [results]=optimize_template(template_file,template_path,varargin)
%% Syntax
% [results]=optimize_template(template_file,template_path,target_input_method,...
%                  target_file,target_path,distractor_input_method,distractor_file,...
%                 distractor_path,distractor_factor,varargin)
%              
%% Inputs and Supplemental Inputs  
%
% template_file - name of the template file you want to optimize
%
% template_path - location of the template file
%
% target_input_method - where should the function find the target chunks. two
% options. 'template_metadata': load the syll associated chunks for the length and 
% the index listed in the template metadata. 'file_input' :load 
% the target chunks from a given mat file.
%
% if input_method = 'file_input'
% target_file - name of the file containing target chunks
% 
% target_path - location of the file containing target chunks 
%
% distractor_input_method -where should the function find the distractor chunks. two
% options. 'template_metadata': load the syll associated chunks for all the 
% listed distractor syllables in the template metadata. 'file_input' :load 
% the distractor chunks from a given mat file.
%
% if input_method = 'file_input'
% distractor_file - name of the file containing distractor chunks
% 
% distractor_path - location of the file containing distractor chunks
%
% distractor_factor - the factor by which the probability density function
% of the distractor distances should be multiplied **1
%
% supp_inputs.sigma_init - the initial sigma of the Gaussians used to
% construct the continuous distributions.default = 0.2;
%
% supp_inputs.verify_only - vefify if the initial sigma yields properly shaped
% Gaussians. if it does, don't bother finding the lowest allowable sigma.
% default = 0;
%
% supp_inputs.collation_function - function used to collate individual
% gaussians and construct the continuous density distributions. default = @mean.
% An alternative can be @sum. 
%
% supp_inputs.target_factor - the factor by which the probability density function
% of the target distances is multipled. similar to distractor_factor.**2. default=1.  
%
% supp_inputs.optimize_in_0_1_box - should the values of the optimized
% template be restricted to between 0 and 1. default=0 **3
%
% supp_inputs.plotting_on - should the dynamic plotting of the progress of
% optimization be displayed. default=1 
%
% supp_inputs.c and supp_inputs.gamma - multiplicative factors that help
% determine the step size in the back-tracking line search. defaults:c=0.1
% and gamma = 0.65
%
% supp_inputs.history_length - the number of previous iterations to
% consider while comparing  termination parameters to the termination condition. default: 10 
%
% supp_inputs.history_change_thr_grad_mag - the lower-bound threshold below which 
% the termination parameter (mean absolute deviation of gradient magnitude as a 
% percentage of the mean gradient magnitude over the last supp_inputs.history_length iterations )
% has to go to satisfy part of the termination condition. default = 0.5
% 
% supp_inputs.history_change_thr_error - the lower-bound threshold below which 
% the termination parameter (mean absolute deviation of error as a 
% percentage of the mean error over the last supp_inputs.history_length iterations )
% has to go to satisfy part of the termination condition. default = 0.5
% 
% supp_inputs.init_stepsize - initial stepsize to be used for determining
% the actual step size using back-tracking line search. this is just for the purpose of 
% making the use of that value (1) explicit. it is not to be changed. default = 1
% 
% supp_inputs.show_decision_params - should the termination parameters and
% how they compare to the termination condition be displayed at each
% iteration. default = 1;
%
% supp_inputs.iter_limit - the maximum number of iterations allowed. after
% this number is reached, the code simply retuns the values that exist at
% that iteration and switch the flag force_stopped to 1. default =  =1000.  
% 
% supp_inputs.grad_mag_lim - the lower-bound threshold below which 
% the termination parameter (gradient magnitude)
% has to go to satisfy part of the termination condition. default = 10^-5;
%
%% Computation   
% 
% Modifies the template in a manner dictated by the gradient descent algorithm
% to achieve minimum overlap (error) between the distance distributions of
% target and distractor chunks. It uses back-tracking line search to calculate 
% the step size. 
%
%
% 
%
%% Outputs  
%
% optim_archive.error_vec= vector containing magnitude of error through all
% the steps of the optimization
%
% optim_archive.threshold_vec=vector containing value of the threshold through all
% the steps of the optimization
%
% optim_archive.fne_instances_vec=vector containing # of false negative chunks through all
% the steps of the optimization
%
% optim_archive.fpe_instances_vec=vector containing # of false positive chunks through all
% the steps of the optimization
%
% optim_archive.fne_rate_vec=vector containing % of false negative chunks through all
% the steps of the optimization
%
% optim_archive.fpe_rate_vec=vector containing % of false positive chunks through all
% the steps of the optimization
%
% optim_archive.fin_sigma_vec=vector containing value of the sigma used through all
% the steps of the optimization
%
% optim_archive.template_grad_vec=matrix containing the gradient vector through all
% the steps of the optimization
%
% optim_archive.template_vec=matrix containing the template vector through all
% the steps of the optimization
% 
% optim_archive.template_grad_mag_vec=vector containing magnitude of the gradient of the template
% through all the steps of the optimization
%
%% Assumptions
%
%
%
%
%
%% Version and Author Identity Notes  
% Last modified by Anand S Kulkarni on 4/25/17
%
% This is the most current version. It was quite heavily modified compared
% to the previous version. Spcifically the step size chosing part was added
% and the function for calculating the gradient was written (by taking code distributed thru the earliert version).
% 
% 
% previous version: optimize_template_v1
% next version: 
%% Related procedures and functions 
% 
%
%
%
%% Detailed Notes
% 
%**1 - this effectively determines where the threshold gets placed by
% shifting the crossing point of the two distributions (in the valley). lower 
% values shift the crossing point to the right and higher values shift it
% to the left. It also determines the relative contribution (relative to target chunks)
% of the distractors to the gradient and the error calculation. 
% 
%**2 - the target factor is by default kept at 1. Moving around of threhsold is
% achieved by controlling the distractor_factor itself. Other aspects involving relative 
% contributions of target and distractor chunks are also similarly
% controlled. 
%
%**3 - this is useful when the version of evtaf that you are using
%normalizes the template by default
%
%% Processing inputs and beginning stuff

% putting in a stop for easier debugging
dbstop if error

% processing mandatory inputs
narg_min=2;



% Assigning default values to supplementary inputs
supp_inputs.sigma_init=0.2;
supp_inputs.verify_only=0;
supp_inputs.collation_function=@mean;
supp_inputs.target_factor=1;
supp_inputs.distractor_factor=1;
supp_inputs.optimize_in_0_1_box=0;
supp_inputs.plotting_on=0;
supp_inputs.c=0.1;
supp_inputs.gamma=0.65;
supp_inputs.history_length=10;
supp_inputs.history_change_thr_grad_mag=0.5;
supp_inputs.history_change_thr_error=0.5;
supp_inputs.init_stepsize=1;
supp_inputs.show_decision_params=0;
supp_inputs.iter_limit=1000;
supp_inputs.grad_mag_lim=10^-5;
supp_inputs.skip_optimization=0;
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

if ~isempty(template_path)
   if ~strcmpi(template_path(end),filesep)
        template_path=[template_path,filesep];
   end
end

%% Body of the function

distractor_factor=supp_inputs.distractor_factor;
% loading template, target_chunks, and distractor_chunks

template=load([template_path template_file]);

[target_chunks]=assemble_target_chunks(template_file,template_path);
[distractor_chunks]=assemble_distractor_chunks(template_file,template_path,...
                'exclude_distractor_syll',supp_inputs.exclude_distractor_syll,...
                'include_gaps',supp_inputs.include_gaps);


% initiating figure and axes
if supp_inputs.plotting_on
    grand_fig=figure;
    densities_axes=subplot(1,3,1,'parent',grand_fig);
    error_axes=subplot(1,3,2,'parent',grand_fig);
    t_grad_mag_axes=subplot(1,3,3,'parent',grand_fig);
end
% initaitng archive variables
optim_archive.error_vec=[];
optim_archive.threshold_vec=[];
optim_archive.fne_instances_vec=[];
optim_archive.fpe_instances_vec=[];
optim_archive.fne_rate_vec=[];
optim_archive.fpe_rate_vec=[];
optim_archive.fin_sigma_vec=[];
optim_archive.template_grad_vec=[];
optim_archive.template_vec=[];
optim_archive.template_grad_mag_vec=[];

% other initializations
first_entry=1;
iteration_no=0;

% PV pairs for the function to be called inside the loop
param1='verify_only';
val1=supp_inputs.verify_only;
param2='sigma_init';
val2=supp_inputs.sigma_init;
param3='target_factor';
val3=supp_inputs.target_factor;
param4='collation_function';
val4=supp_inputs.collation_function;


while true        
    iteration_no=iteration_no+1;
    if iteration_no>supp_inputs.iter_limit
        break
    end
    
    iter_vec=(1:iteration_no);
   
    
    if supp_inputs.verify_only==1   
        if first_entry==1
             param5='calc_gradient';
             val5=1;
             param6='curr_thr';
             val6=-144;% just a specific value it can be tested against
             param7='curr_sigma';
             val7=-144;
             param8='curr_iter';
             val8=iteration_no;
            [sigma,density_target,density_distractors,dists_vec,...
            threshold,fne_instances,fpe_instances,fne_rate,fpe_rate,...
            total_error,template_grad,template_grad_mag]=calculate_error_function_and_gradient(template,target_chunks,distractor_chunks,distractor_factor,...
                                                     param1,val1,param2,val2,param3,val3,param4,val4,param5,val5,param6,val6,param7,val7,param8,val8);
                                                     
        elseif first_entry==0
            param5='calc_gradient';
             val5=1;
              param6='curr_thr';
            val6=threshold;% just a specific value it can be tested against
             param7='curr_sigma';
             val7=sigma_init;
             param8='curr_iter';
             val8=iteration_no;
            [sigma,density_target,density_distractors,dists_vec,...
            threshold,fne_instances,fpe_instances,fne_rate,fpe_rate,...
            total_error,template_grad,template_grad_mag]=calculate_error_function_and_gradient(template,target_chunks,distractor_chunks,distractor_factor,...
                                                    param1,val1,param2,val2,param3,val3,param4,val4,param5,val5,param6,val6,param7,val7,param8,val8);
                                                     
        else
                   error('Incorrect value for supp_inputs.verify_only')           
        end
    elseif supp_inputs.verify_only==0
        if first_entry==1
            param5='calc_gradient';
            val5=1;
            param6='curr_thr';
            val6=-144;% just a specific value it can be tested against
             param7='curr_sigma';
             val7=-144;
             param8='curr_iter';
             val8=iteration_no;
            [sigma,density_target,density_distractors,dists_vec,...
        threshold,fne_instances,fpe_instances,fne_rate,fpe_rate,...
        total_error,template_grad,template_grad_mag]=calculate_error_function_and_gradient(template,target_chunks,distractor_chunks,distractor_factor,...
                                                     param1,val1,param2,val2,param3,val3,param4,val4,param5,val5,param6,val6,param7,val7,param8,val8);
                                                     
        elseif first_entry==0
            try
                param1='verify_only';
                val1=1;
                param2='sigma_init';
                val2=sigma_init;
                param5='calc_gradient';
                val5=1;
                param6='curr_thr';
                val6=threshold;% just a specific value it can be tested against
                 param7='curr_sigma';
                 val7=sigma_init;
                 param8='curr_iter';
                 val8=iteration_no;
                [sigma,density_target,density_distractors,dists_vec,...
                threshold,fne_instances,fpe_instances,fne_rate,fpe_rate,...
                total_error,template_grad,template_grad_mag]=calculate_error_function_and_gradient(template,target_chunks,distractor_chunks,distractor_factor,...
                                                        param1,val1,param2,val2,param3,val3,param4,val4,param5,val5,param6,val6,param7,val7,param8,val8);
                                                     
            catch err_exp
                if isequal(err_exp.message,'The sigma supplied is not sufficient for constructing smooth, monotonous density distributions')
                    param1='verify_only';
                    val1=0;
                    param2='sigma_init';
                    val2=sigma_init;
                    param5='calc_gradient';
                    val5=0;
                     param6='curr_thr';
                    val6=threshold;% just a specific value it can be tested against
                     param7='curr_sigma';
                     val7=sigma_init;
                     param8='curr_iter';
                     val8=iteration_no;
                    [sigma,~,~,~,...
                    ~,~,~,~,~,...
                    ~]=calculate_error_function_and_gradient(template,target_chunks,distractor_chunks,distractor_factor,...
                     param1,val1,param2,val2,param3,val3,param4,val4,param5,val5,param6,val6,param7,val7,param8,val8);
                                                     
                    
                    sigma_init=sigma;
                    template=optim_archive.template_vec(:,1);
                    iteration_no=1;
                    iter_vec=(1:iteration_no);
                   % reinitiating
                    optim_archive.error_vec=[];
                    optim_archive.threshold_vec=[];
                    optim_archive.fne_instances_vec=[];
                    optim_archive.fpe_instances_vec=[];
                    optim_archive.fne_rate_vec=[];
                    optim_archive.fpe_rate_vec=[];
                    optim_archive.fin_sigma_vec=[];
                    optim_archive.template_grad_vec=[];
                    optim_archive.template_vec=[];
                    optim_archive.template_grad_mag_vec=[];
                    
                    param1='verify_only';
                    val1=1;
                    param2='sigma_init';
                    val2=sigma_init;
                    param5='calc_gradient';
                    val5=1;
                     param6='curr_thr';
                    val6=-144;% since iteration_no is set to one above this does not matter
                     param7='curr_sigma';
                     val7=-144;
                     param8='curr_iter';
                     val8=iteration_no;
                    [sigma,density_target,density_distractors,dists_vec,...
                    threshold,fne_instances,fpe_instances,fne_rate,fpe_rate,...
                    total_error,template_grad,template_grad_mag]=calculate_error_function_and_gradient(template,target_chunks,distractor_chunks,distractor_factor,...
                    param1,val1,param2,val2,param3,val3,param4,val4,param5,val5,param6,val6,param7,val7,param8,val8);
                                                     

                
                else
                    throw(err_exp)
                end
            end
        else
            error('Incorrect value for first_entry')
        end
    else
       error('Incorrect value for supp_inputs.verify_only')  
    end
    % flagging the improper shape of disttributions                                         
    if threshold==-1
       error('The shape of the distributions is not right')
    end
    % setting the sigma_init of the next iteration as the current optimal
    % sigma
    sigma_init=sigma;
       
    
    % archiving 
    optim_archive.template_vec=[optim_archive.template_vec,template];
    optim_archive.threshold_vec=[optim_archive.threshold_vec,threshold];
    optim_archive.fin_sigma_vec=[optim_archive.fin_sigma_vec,sigma];
    optim_archive.fne_instances_vec=[optim_archive.fne_instances_vec,fne_instances];
    optim_archive.fpe_instances_vec=[ optim_archive.fpe_instances_vec,fpe_instances];
    optim_archive.fne_rate_vec=[optim_archive.fne_rate_vec,fne_rate];
    optim_archive.fpe_rate_vec=[optim_archive.fpe_rate_vec,fpe_rate];
    optim_archive.error_vec=[optim_archive.error_vec,total_error];
    optim_archive.template_grad_vec=[optim_archive.template_grad_vec,template_grad];  
    optim_archive.template_grad_mag_vec=[optim_archive.template_grad_mag_vec,template_grad_mag];
     
    
    % plotting densities and threshold   
    if supp_inputs.plotting_on
        if first_entry==1 
            densities_in_h=plot(densities_axes,dists_vec,density_target,'g');
            hold(densities_axes,'on')
            densities_out_h=plot(densities_axes,dists_vec,density_distractors,'r');
           
            th_ylim=[0,max([density_target;density_distractors])];
            th_xlim=[threshold,threshold];
            th_line_h=plot(densities_axes,th_xlim,th_ylim,'color','k','linestyle','--');
            hold(densities_axes,'off')
            title(densities_axes,'Density of Distance Distributions and the Threshold')
            set(densities_in_h,'XDataSource','dists_vec','YdataSource','density_target');
            set(densities_out_h,'XDataSource','dists_vec','YdataSource','density_distractors');
            set(th_line_h,'Xdatasource','th_xlim','Ydatasource','th_ylim');
            
            err_line_h=plot(error_axes,optim_archive.error_vec,'*-'); 
            title(error_axes,'Total Error')
            set(err_line_h,'ydatasource','optim_archive.error_vec','xdatasource','iter_vec');
            
            template_grad_mag_h=plot(t_grad_mag_axes,optim_archive.error_vec,'*-'); 
            title(t_grad_mag_axes,'Template Gradient Magnitude')
            set(template_grad_mag_h,'ydatasource','optim_archive.template_grad_mag_vec','xdatasource','iter_vec');
        elseif first_entry==0 
            
            % refreshing
            th_ylim=[0,max([density_target;density_distractors])];
            th_xlim=[threshold,threshold];
            refreshdata(grand_fig,'caller')
            drawnow  
        else
           error('Incorrect value for first_entry') 
        end
    end       
   
    % calculating termination decision variables
    
    if iteration_no<=supp_inputs.history_length
        start_ind=1;
    else
        start_ind=iteration_no-supp_inputs.history_length+1;
    end
    
    hist_change_grad_mag=mad(optim_archive.template_grad_mag_vec(start_ind:end))/mean(optim_archive.template_grad_mag_vec(start_ind:end))*100;
    hist_change_error=mad(optim_archive.error_vec(start_ind:end))/mean(optim_archive.error_vec(start_ind:end))*100;
    % resetting the flag
    if first_entry==1
        first_entry=0;
    end
    % termination 	
    condition_1=hist_change_grad_mag<supp_inputs.history_change_thr_grad_mag && hist_change_grad_mag>0 &&...
       hist_change_error<supp_inputs.history_change_thr_error && hist_change_error>0; % calculating mad makes this condition (>0) redundant
    condition_2=optim_archive.template_grad_mag_vec(end)<supp_inputs.grad_mag_lim &&...
     hist_change_error<supp_inputs.history_change_thr_error && hist_change_error>0; % calculating mad makes this condition (>0) redundant
 
    if condition_1||condition_2
       if supp_inputs.show_decision_params==1     
            dd={'iteration_no','grad_mag','hist_change_grad_mag','hist_change_error','template_step_mag';...
                      iteration_no,template_grad_mag,hist_change_grad_mag,hist_change_error,'N.A., conditions met'}        
        elseif supp_inputs.show_decision_params~=0
            error('The value of supp_inputs.show_decision_params should be 0 or 1')        
        end
        break    
    elseif supp_inputs.skip_optimization==1
        break
    else        
        % determining the template_step_mag using line search
        fc=total_error;
        delta_f=supp_inputs.c*template_grad_mag;
        ttt=supp_inputs.init_stepsize;
        while true
%                disp('rr') 
%             end
            param5='calc_gradient';
            val5=0;
            param2='sigma_init';
            val2=sigma_init;
             param6='curr_thr';
            val6=threshold;% just a specific value it can be tested against
             param7='curr_sigma';
             val7=sigma_init;
             param8='curr_iter';
             val8=iteration_no;
            
            try
                param1='verify_only';
                val1=1;   
                %[sigma_init,~,~,~,threshold,~,~,~,~,total_error]
                [~,~,~,~,~,~,~,~,~,total_error]=calculate_error_function_and_gradient(template-(ttt*template_grad),...
                                                                target_chunks,distractor_chunks,distractor_factor,...
                                                            param1,val1,param2,val2,param3,val3,param4,val4,param5,val5,param6,val6,param7,val7,param8,val8);
                                                     
            catch err_exp
                if isequal(err_exp.message,'The sigma supplied is not sufficient for constructing smooth, monotonous density distributions')
                    param1='verify_only';
                    val1=0;
                    [~,~,~,~,~,~,~,~,~,total_error]=calculate_error_function_and_gradient(template-(ttt*template_grad),...
                                                                target_chunks,distractor_chunks,distractor_factor,...
                                                             param1,val1,param2,val2,param3,val3,param4,val4,param5,val5,param6,val6,param7,val7,param8,val8);
                                                     
                else
                    throw(err_exp)
                end
                
                
            end
                                                        
           new_f=total_error; 
           if new_f>fc+ttt*delta_f
               ttt=ttt*supp_inputs.gamma;             
           else
              template_step_mag=ttt;
              break
           end           
        end
        % recalculating the template
        template=template-template_step_mag*template_grad;        
        if supp_inputs.optimize_in_0_1_box
                template(template<0)=0;
                template(template>1)=1;
        end
        
    end
    if supp_inputs.show_decision_params==1     
        dd={'iteration_no','grad_mag','hist_change_grad_mag','hist_change_error','template_step_mag';...
                  iteration_no,template_grad_mag,hist_change_grad_mag,hist_change_error,template_step_mag}        
    elseif supp_inputs.show_decision_params~=0
        error('The value of supp_inputs.show_decision_params should be 0 or 1')        
    end
        
    
end

%% Processing outputs 

results.optim_archive=optim_archive;
if iteration_no>supp_inputs.iter_limit
    results.optim_archive.force_stopped=1;
else
    results.optim_archive.force_stopped=0;
end    
results.templatefullfile=[template_path template_file];

