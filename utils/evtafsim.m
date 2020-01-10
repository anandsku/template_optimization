function [vals,spec,real_spec] = evtafsim(rsong,fs,nfft,templates,OLDWAY,USEERROR,normalize_temps,take_sqrt);
%[vals,spec] = evtafsim(rsong,fs,templates,OLDWAY,USEERROR);
% returns the threshold in case it was set in the program
% simtaf with the clocks and multiple templates put in
% templates is a matrix with each column is one template
% nfft size of chunk whose fft needs to be taken
% normalize_temps - tells evtafsim  


if (~exist('USEERROR','var'))
    USEERROR=0;
end
if (~exist('OLDWAY','var'))
    OLDWAY=0;
end

if (~exist('normalize_temps','var'))
    normalize_temps=1;
end

if (~exist('take_sqrt','var'))
    take_sqrt=0;
end

blen=nfft/2;
hammy  = hamming(2*blen);

ntempl = size(templates,2);
%nrep = floor(length(rsong)/nfft)-1; % so here you exclude the last partial chunk and a full chunk before that
nrep = floor(length(rsong)/nfft);

vals = zeros([nrep,ntempl]);
spec = zeros([nrep,blen]);
real_spec=zeros(size(spec));

% zeroing the first six frequencies in the template... just in case
for jj = 1:ntempl        
    templates(1:6,jj) = 0;           
end


%templates = templates;
% templates are normalized just in case
if normalize_temps
    for jj = 1:ntempl
        if (OLDWAY==0)
            templates(:,jj) = templates(:,jj)-min(templates(:,jj));
            templates(:,jj) = templates(:,jj)./max(templates(:,jj));
        else
            normtmp = templates(:,jj).'*templates(:,jj);
            templates(:,jj) = templates(:,jj)./sqrt(normtmp);
        end
    end
end

for ii = 1:nrep
	ind1 = (ii-1)*nfft+ 1;
	ind2 = ind1 + nfft - 1;
	datchunk = rsong(ind1:ind2) - mean(rsong(ind1:ind2));
	fdatchunk = abs(fft(hammy.*datchunk));
    %fdatchunk = abs(fft(datchunk));

	sp = abs(fdatchunk(1:blen));
    real_spec(ii,:)=sp.';
    if (USEERROR==0)
        sp(1:6)=0.0;
    else
        %sp(2:end)=sp(1:end-1);
        %sp(1:6)=0.0;
        sp = [zeros([6,1]);sp(6:end-1)];
    end
    if (OLDWAY==1)
        normtmp = sqrt(sp.'*sp);
        sp = sp./normtmp;
    else
        sp = sp-min(sp);
        sp = sp./max(sp);
    end
    for jj = 1:ntempl
        if (OLDWAY==1)
            vals(ii,jj) = acos(sp.'*templates(:,jj));
        else
            if take_sqrt
             vals(ii,jj) =sqrt((sp-templates(:,jj)).'*(sp-templates(:,jj)));
            else
             vals(ii,jj) = (sp-templates(:,jj)).'*(sp-templates(:,jj));
            end
        end
    end
    spec(ii,:) = sp.';
end
return;
