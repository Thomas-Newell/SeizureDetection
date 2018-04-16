%THIS SCRIPT WILL RUN MANUALLY TO DO THE AUTOMATED SEIZURE DETECTION!!!

clear,clc
direc = uigetdir('Select the directory to analyze');        %Selects the drectory for testing.
disp('Directory Started')                                   %Displays that analysis has been initiated.
dirlist = dir(direc);                                       %Displays that analysis has been initiated.
for k = 3:size(dirlist,1)
    name = dirlist(k,1).name;
    if strcmp(name(end-2:end),'acq') == 1
        n = name;
        break
    end
end
sss = dir([direc '\' n]);
if sss.bytes < 7.5986e+06 % File must be > 5 mins in length (this many bytes).
    disp('ACQ file is too small to analyze!')
    break
end
if exist('n','var') == 0
    error('No ACQ file in directory')
end
keep('direc','n');
disp(['Seizure Detection initiated on ' date ' for File: ' n]);
%% ANALYSIS
started = 1;
info =  acqreader07092013([direc '\' n]);
ACe = 3.5; FFTe = 1; LLe = 2; SSe = 2.5;
ACw = 0.4; FFTw = 2.5; LLw = 2; SSw = 1.9;
tic;
% AC
AC = AC_V1([direc '\' n],25,1,60);
%FFT
FFT = FFT_V1([direc '\' n],25,1,1);
%LL
LL = LL_V1([direc '\' n],25,1);
%SS
SS = SS_V1([direc '\' n],25,1,7.6e10);
duration = toc;
disp(['Analysis took ' num2str(duration/60) ' minutes']);
%Interpolate...
smallest = 1:min([size(AC,2),size(FFT,2),size(LL,2),size(SS,2)]);
biggest = 1:max([size(AC,2),size(FFT,2),size(LL,2),size(SS,2)]);
PLAY = zeros(size(AC,1),smallest(end));
if size(AC,2) > smallest(end)
    for int = 1:size(LL,1)
        PLAY(int,:) = interp1(biggest,AC(int,:),smallest);
    end
    AC = PLAY;
end
PLAY = zeros(size(AC,1),smallest(end));
if size(FFT,2) > smallest(end)
    for int = 1:size(LL,1)
        PLAY(int,:) = interp1(biggest,FFT(int,:),smallest);
    end
    FFT = PLAY;
end
PLAY = zeros(size(AC,1),smallest(end));
if size(LL,2) > smallest(end)
    for int = 1:size(LL,1)
        PLAY(int,:) = interp1(biggest,LL(int,:),smallest);
    end
    LL = PLAY;
end
PLAY = zeros(size(AC,1),smallest(end));
if size(SS,2) > smallest(end)
    for int = 1:size(LL,1)
        PLAY(int,:) = interp1(biggest,SS(int,:),smallest);
    end
    SS = PLAY;
end
%Normalizing
empchan = find(info.ChannelNames == 'e','E');
for ec = 1:size(empchan,1)
    AC(empchan(ec),:) = 0;       %PREALLOCATE THIS SOMEHOW...Save 11ish seconds...
    FFT(empchan(ec),:) = 0;
    LL(empchan(ec),:) = 0;
    SS(empchan(ec),:) = 0;
end
ACNEWnorm1 = AC;        %Preallocating...
FFTNEWnorm1 = FFT;      %Preallocating...
LLNEWnorm1 = LL;        %Preallocating...
SSNEWnorm1 = SS;        %Preallocating...
for C = 1:size(AC,1)
    ACNEWnorm1(C,:) = (AC(C,:)-min(AC(C,:)))...
        /(max(AC(C,:))-min(AC(C,:)));     %The max value in that channel is 1.
    FFTNEWnorm1(C,:) = (FFT(C,:)-min(FFT(C,:)))...
        /(max(FFT(C,:))-min(FFT(C,:)));     %The max value in that channel is 1.
    LLNEWnorm1(C,:) = (LL(C,:)-min(LL(C,:)))...
        /(max(LL(C,:))-min(LL(C,:)));     %The max value in that channel is 1.
    SSNEWnorm1(C,:) = (SS(C,:)-min(SS(C,:)))...
        /(max(SS(C,:))-min(SS(C,:)));     %The max value in that channel is 1.
end
for C = 1:size(AC,1)
    if isnan(ACNEWnorm1(C,1)) == 1
        ACNEWnorm1(C,:) = 0;
    end
    if isnan(FFTNEWnorm1(C,1)) == 1
        FFTNEWnorm1(C,:) = 0;
    end
    if isnan(LLNEWnorm1(C,1)) == 1
        LLNEWnorm1(C,:) = 0;
    end
    if isnan(SSNEWnorm1(C,1)) == 1
        SSNEWnorm1(C,:) = 0;
    end
end
%COMB
COMB = ( (ACNEWnorm1*ACw).^ACe + (FFTNEWnorm1*FFTw).^FFTe+(LLNEWnorm1*LLw).^LLe+(SSNEWnorm1*SSw).^SSe )...
    / (ACw.^ACe + FFTw.^FFTe + LLw.^LLe + SSw.^SSe);
%% SPIKEDETECTION
THRESH = 2.5;
cutoff = 5;
[combDETECTIONS, info, combszseconds] = SzDetectAlt(COMB,info,THRESH,cutoff);
%% FILE WRITING
fclose('all');
for tv = 1:size(COMB,1)
    Scopy = COMB(tv,:);
    Scopy(find(Scopy > cutoff*mean(Scopy))) = cutoff*mean(Scopy);   
    thrval(tv,1) = THRESH*mean(Scopy);
end
DetFile = [direc '\' n(1:end-8) '_Version3_' date '.det'];
DetFID = fopen(DetFile,'w');
if size(combszseconds,2) > 0
    combszseconds = sortrows(combszseconds,2);
    for ccc = 1:size(combszseconds,1)
        fprintf(DetFID,['%d, ' num2str(floor(combszseconds(ccc,2))) '\n'],combszseconds(ccc,1));
    end;
end
fclose(DetFID);
% RETURN RESULTS FIGURE
clf
fl = (info.EndOfFileInHours);
timevar = linspace(0,fl,size(COMB,2));
for ccc=1:size(AC,1)
    plot(timevar, COMB(ccc,:));
    line([timevar(1) timevar(end)],[thrval(ccc,1) thrval(ccc,1)],'color','g','linewidth',3);
    xlabel('time'), ylabel('COMB');title([n ' channel ' num2str(ccc)]);
    text(1,.8, [info.ChannelNames(ccc,:)]);
    if mean(COMB(ccc,:)) > 0;
        axis([0 info.EndOfFileInHours -.1*min(COMB(ccc,:)) 1.1*(max(COMB(ccc,:)))])
    else
        axis([0 info.EndOfFileInHours -.1 1.1])
    end
    saveas(gcf,[direc '\' n(1:end-8) '_SzDetectionV3_Chan_' num2str(ccc) '_results'],'png')
    close('all')
end;