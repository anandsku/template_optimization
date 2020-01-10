function [sigma,density_target,density_distractors,dists_vec,...
       threshold,fne_instances,fpe_instances,fne_rate,fpe_rate,...
       total_error,varargout]=calculate_error_function_and_gradient(template,target_chunks,distractor_chunks,distractor_factor,varargin)
%% Syntax
%
% 
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
narg_min=3;


% % Assigning default values to supplementary inputs
supp_inputs.verify_only=1;
supp_inputs.sigma_init=0.3;
supp_inputs.target_factor=1;
supp_inputs.collation_function=@mean;
supp_inputs.calc_gradient=1;
supp_inputs.curr_thr=[];
supp_inputs.curr_sigma=[];
supp_inputs.curr_iter=[];
supp_inputs.exfac=9;


if nargin<narg_min
    error(['The number of inputs should at least be ' narg_min])
else
    % processing supplementary inputs
    supp_inputs=parse_pv_pairs(supp_inputs,varargin);
end

%% Body of the function
    
% calculating distances
dists_target=sqrt(sum((target_chunks(:,:)-repmat(template(:,:),1,size(target_chunks(:,:),2))).^2,1));
dists_distractor=sqrt(sum((distractor_chunks(:,:)-repmat(template(:,:),1,size(distractor_chunks(:,:),2))).^2,1));

% calculating optimal sigma/ verifying optimal sigma
[sigma,density_target,density_distractors,dists_vec]...
=calculate_optimal_sigma(dists_target,dists_distractor,'verify_only',supp_inputs.verify_only,...
'sigma_init',supp_inputs.sigma_init);

sigma_target=sigma;
sigma_distractor=sigma;


% calculating optimal threshold
[threshold]=calculate_optimal_threshold(density_target,density_distractors,dists_vec,distractor_factor); 

% multiplying the density distributions with appropriate factors 
density_target=density_target*supp_inputs.target_factor;
density_distractors=density_distractors*distractor_factor;

%  error rates and instances
fne_instances=length(find(dists_target>=threshold));
fpe_instances=length(find(dists_distractor<=threshold));
fne_rate=length(find(dists_target>=threshold))/length(dists_target);
fpe_rate=length(find(dists_distractor<=threshold))/length(dists_distractor);



% calculating total error 
fn_error=supp_inputs.collation_function(1-cdf('norm',threshold.*ones(size(dists_target)),dists_target,sigma_target*ones(size(dists_target))));
fp_error=supp_inputs.collation_function(cdf('norm',threshold.*ones(size(dists_distractor)),dists_distractor,sigma_distractor*ones(size(dists_distractor))));
total_error=supp_inputs.target_factor*fn_error+distractor_factor*fp_error;

varargout={};
if supp_inputs.calc_gradient==1
 % template gradient calculation
    template_grad_fne=supp_inputs.collation_function((repmat(template,1,size(target_chunks,2))-target_chunks) .*repmat(pdf('norm',threshold*ones(size(dists_target)),...
        dists_target,sigma_target*ones(size(dists_target)))./dists_target,size(target_chunks,1),1),2);% fne=false negative error

    template_grad_fpe=supp_inputs.collation_function((repmat(template,1,size(distractor_chunks,2))-distractor_chunks) .*repmat(pdf('norm',threshold*ones(size(dists_distractor)),...
         dists_distractor,sigma_distractor*ones(size(dists_distractor)))./dists_distractor,size(distractor_chunks,1),1),2);% fpe=false positive error

    template_grad=supp_inputs.target_factor*template_grad_fne-distractor_factor*template_grad_fpe;
    template_grad_mag=sqrt(sum(template_grad.^2));
    varargout{1}=template_grad;
    varargout{2}=template_grad_mag;
elseif supp_inputs.calc_gradient~=0
   error('The value of supp_inputs.calc_gradient should be either 1 or 0') 
end
