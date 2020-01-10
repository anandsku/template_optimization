function [supp_inputs]=process_supplementary_inputs(supp_inputs)

% this function processes supplementary inputs in case of manual inputs
% instances of a function. 

if nargin~=1
    error('This function needs exactly one input argument')
end


fprintf('The following default values are being used for supplementary inputs\n')
supp_inputs
 change_required=input('\nDo you want to change any of them? y or n \n-->  ','s');
 if strcmpi(change_required,'y')
     supp_inp_fields=fieldnames(supp_inputs);
     no_fields=length(supp_inp_fields);
     
     fields_to_change=input(['\nEnter the names of the supplementary input fields that you want to change?'...
                      '\nNote that this input is a cell array  e.g. {''' supp_inp_fields{1} ''',''' supp_inp_fields{no_fields} '''} \n-->  ']);     
     no_fields_to_change=length(fields_to_change);             
     for i=1:no_fields
         for j=1:no_fields_to_change
             if strcmpi(supp_inp_fields{i},fields_to_change{j}) % if this is a field to change         
                 response3=input(['Enter the new value for the supplementary input ' supp_inp_fields{i} '. Enter ''leave_blank'' to leave it blank or hit Enter to leave it unchanged \n-->  '],'s');
                 if ~isempty(response3) % if the response if not empty change the default
                     if strcmpi(response3,'leave_blank')                
                        supp_inputs.(supp_inp_fields{i})='';    
                     else             
                         if ~isempty(supp_inputs.(supp_inp_fields{i})) % if the default value is not empty
                             if ischar(supp_inputs.(supp_inp_fields{i})) % if the default value is a string
                                 string_resp=1;
                             else % if the default value is a not a string
                                  string_resp=0;
                             end
                         else % if the default value is empty, ask user how he wants his entry to be interpreted
                            response4=input('Do you want your entry to be interpreted as a string? 1 or 0 \n-->');
                            if response4==1
                                string_resp=1;
                            elseif response4==0
                                string_resp=0;
                            else
                               error('Your response should be either 1 or 0') 
                            end

                         end

                         if string_resp==1
                                supp_inputs.(supp_inp_fields{i})=response3;
                         else
                                supp_inputs.(supp_inp_fields{i})=eval(response3);
                         end
                     end

                 end
             end
         end    
         
     end  
     fprintf('This is how the values for supplementary inputs look like after the change\n')
     supp_inputs

 elseif strcmpi(change_required,'n')
     fprintf('You chose not to change the default values of supplementary inputs\n\n')
 else
    error('Your response should be either y or n') 
 end

