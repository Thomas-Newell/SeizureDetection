%Automated Sz Detection V3
%Written by Thomas Newell
%Last Update: 10/27/2015

clear,clc
NOW = clock;

%Every directory to be analyzed MUST be in the cell structure below!!!!!!
% paths = {'\\3222C-RECORDING\Data','\\3222C-recording\Data2','\\3222C-RECORDING\Data3','\\B0335-EEG1\Data','\\B0335-EEG1\Data2','\\B0335-EEG1\Data3','\\B0335-EEG2\Data','\\B0335-EEG2\Data2','\\B0335-EEG2\Data3'};

paths = {'\\3222C-RECORDING\Data','\\3222C-RECORDING\Data2','\\3222C-RECORDING\Data3','\\B0335-EEG1\Data','\\B0335-EEG1\Data2','\\B0335-EEG1\Data3','\\B0335-EEG2\Data','\\B0335-EEG2\Data2','\\B0335-EEG2\Data3','\\EEG_MOUSE_BLA1\Data','\\EEG_MOUSE_BLA1\Data2','\\EEG_MOUSE_BLA1\Data3','\\MOUSE-EEG\Data','\\MOUSE-EEG\Data2','\\MOUSE-EEG\Data3'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for PATH = 1:size(paths,2)      %Goes through each directory in 'paths'
    keep('NOW','paths','PATH')
    direc = paths{1,PATH};
    disp(['Starting analysis of ' direc]);
    foldlist = dir(direc);
    fold = 3;
    while fold < size(foldlist,1) + 1    %Make some fix here so that it only leaves the last folder up on the recording directory
        folddate = foldlist(fold,1).name;
        direc = paths{1,PATH};
        folder = [direc '\' folddate];
        if size(folddate,2) > 8         %Folder size must be longer than 8... As it should be
            if folddate(9) == '-';      %Ninth index should be a '-' for the ones we're using.
                if (365.25*NOW(1) + 30.4167*NOW(2) + NOW(3) + NOW(4)/24 + NOW(5)/1440 + NOW(6)/86400)...
                        - (365.25*str2double(folddate(1:4)) + 30.4167*str2double(folddate(5:6)) + str2double(folddate(7:8))...
                        + str2double(folddate(10:11))/24 + str2double(folddate(12:13))/1440 + str2double(folddate(14:15))/86400) < 14 %Must be less than 14 days old
                    if NOW(3) ~= str2double(folddate(7:8)) %As long as the date isn't the same...
                        %THERE MUST BE NO .DET FILE
                        dirlist = dir(folder);
                        isDET = -1;
                        for D = 3:size(dirlist,1)
                            if dirlist(D,1).name(end-2:end) == 'det'
                                isDET = D;
                            end
                        end
                        if isDET > 0
                            disp(['There is already a .det file for folder ' folddate]);
                            fold = fold + 1;
                            continue
                        else
                            for k = 3:size(dirlist,1)
                                name = dirlist(k,1).name;
                                if strcmp(name(end-2:end),'acq') == 1
                                    n = name;
                                    break
                                end
                            end
                            if exist('n','var') == 0
                                disp(['No ACQ file in directory ' folder])
                                break
                            end
                            disp(['Seizure Detection initiated on ' date ' for File: ' n]);
                            %% ANALYSIS
                            direc = [direc '\' n(1:15)];
                            sss = dir([direc '\' n]);
                            if sss.bytes < 7.5986e+06 % File must be > 5 mins in length (this many bytes).
                                disp('ACQ file is too small to analyze!')
                                fold = fold + 1;
                                continue
                            end
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
%                             THRESHexp = 2.5;
%                             THRESHexp2 = 2.0;
                            cutoff = 5;
                            [combDETECTIONS, info, combszseconds] = SzDetectAlt(COMB,info,THRESH,cutoff);
                            %% FILE WRITING
                            fclose('all');
                            for tv = 1:size(COMB,1)
                                Scopy = COMB(tv,:);
                                Scopy(find(Scopy > cutoff*mean(Scopy))) = cutoff*mean(Scopy);
                                thrval(tv,1) = THRESH*mean(Scopy);
%                                 thrvalexp(tv,1) = THRESHexp*mean(Scopy);
%                                 thrvalexp2(tv,1) = THRESHexp2*mean(Scopy);
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
%                                 line([timevar(1) timevar(end)],[thrvalexp(ccc,1) thrvalexp(ccc,1)],'color','r','linewidth',3);
%                                 line([timevar(1) timevar(end)],[thrvalexp2(ccc,1) thrvalexp2(ccc,1)],'color','b','linewidth',3);
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
                            nowza = clock;
                            disp(['Finished at ' num2str(nowza(4)) ':' num2str(nowza(5))]);
                        end
%                         fold = fold + 1;
                    else
                        disp([folddate ' is or was being recorded today... This folder will be analyzed tomorrow']);
%                         fold = fold + 1;
                    end
                else
                    disp([folddate ' is too old to review']);
%                     fold = fold + 1; %Keep this commented out. Testing
%                     7/27/2017
                end
            end
%             fold = fold + 1; %Bug fix 4/6/2016 (Comment out if up-to-date on Sz detection)
        end
        fold = fold + 1; %Bug fix 7/21/2017
    end
end
disp('Done with Seizure Detection!!!');