function []=calculate_amplitude_threshold(batchname,batchpath)
%% Inputs:
% batchname - name of the bacth file
% batchpath - location of the batch file
% 
%% Outputs
% It wites a file called amp_thr.mat. The key variable in that file is 
% amp_thr - it is the amplitude value that make sure that the numbers of
% the following kinds of slices are minimized:
% 1. syllable assigned slices with amplitude lower that amplitude threshold
% 2. gap assigned slices with amplitude higher than amplitrude threshold 

filelist=make_filelist_from_batch(batchname,batchpath);
am=[];
pc=[];
extr_inds_hi=true(0);
extr_inds_low=true(0);
equal_2=@(x)isequal(x,true);
for i=1:length(filelist)
    [~,nm,ext]=fileparts(filelist{i});
    simdata_name=[batchpath nm '_simdata.mat'];
    notmat_name=[batchpath nm ext '.not.mat'];
    load(simdata_name)
    load(notmat_name)    
    
    extr_inds_temp_low=simdata.pc_in_syll==0;
    extr_inds_low=[extr_inds_low,extr_inds_temp_low];
    extr_inds_temp_hi=simdata.pc_in_syll==100;
    extr_inds_hi=[extr_inds_hi,extr_inds_temp_hi];
    am=[am,simdata.chunk_amp];
    pc=[pc,simdata.pc_in_syll];
end



fiftyplus=pc>=50;
fiftyminus=pc<50;

minusprc=prctile(am(fiftyminus),(0:1:100));
plusprc=prctile(am(fiftyplus),(0:1:100));



minusall=sort(am(fiftyminus));
plusall=sort(am(fiftyplus));
olapst=minusall(find(minusall>plusall(1),1)); % overlap start
olapen=plusall(find(minusall(end)<plusall,1)); % overlap end
nbins=1000;
testthrs=olapst:(olapen-olapst)/nbins:olapen;
fpr=zeros(1,length(testthrs));
fnr=zeros(1,length(testthrs));
fpc=zeros(1,length(testthrs));
fnc=zeros(1,length(testthrs));

for ki=1:length(fpr)
    currthr=testthrs(ki);
   fpr(ki)=length(find(minusall>currthr))/length(minusall)*100; 
   fnr(ki)=length(find(plusall<currthr))/length(plusall)*100;   
   fpc(ki)=length(find(minusall>currthr)); 
   fnc(ki)=length(find(plusall<currthr));  
end

[~,minind]=min(fpr+fnr);
amp_thr=testthrs(minind);
falsepos=fpr(minind);
falseneg=fnr(minind);
falseposcount=fpc(minind);
falsenegcount=fnc(minind);
save(fullfile(batchpath,'amp_thr.mat'),'amp_thr')

% optional plotting
%{
ax1=axes;
semilogx(minusprc,(0:1:100),'r')
hold on
semilogx(plusprc,(0:1:100),'b')

grid minor
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';



figure
plot(fpr,fnr,'linewidth',2)
hold on

plot(fpr(minind),fnr(minind),'*k')

xlabel('false positive rate')
ylabel('false negative rate')
plot(ax1,[amp_thr,amp_thr],[0,100],'g-')
%}


 