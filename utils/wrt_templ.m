function wrt_templ(TEMPLNAME,templ,ask_info)
% wrt_templ(TEMPLNAME,templ);
%ask_info tells the code if it should ask for template info. if it is == 1,
%the code will ask for the template info and it will also write a
%TEMPLNAME_temp_info.mat file. 
if ask_info==1
    template.targetnt=input('Please enter the target note for this template\n','s');
    template.prents=input('Please enter the pre notes for this template\n','s');
    template.postnts=input('Please enter the post notes for this template\n','s');
    template.tag=input('Please enter the tag for this template\n','s');    
    template.notes=input('Please enter the notes for this template\n','s');
    template.batch=input('Please enter the name of the batch file from which this template was derived\n','s');
    template.birdname=input('Please enter the name of the bird for whom this template was made\n','s');
end

% setting other fields of template to []
template.slices=[];
template.avn_file=[];
template.synshift=[];

    
fid=fopen([TEMPLNAME '.dat'],'w');
for ii=1:size(templ,1)
	fprintf(fid,'%.5e',templ(ii,1));
	for jj=2:size(templ,2)
		fprintf(fid,' %.5e',templ(ii,jj));
	end
	fprintf(fid,'\n');
end
fclose(fid);

if ask_info==1
    save([TEMPLNAME '_temp_info.mat'],'template')
end

return;
