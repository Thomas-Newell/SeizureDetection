function [DETECTIONS, info, szseconds] = SzDetectAlt(SPIKES,info,absminthresh,cutoff)
%This function is what will actually detect the 'spikes' in the results
%from A SINGLE algorithm. This will be inserted into each algorithm 
%function and will returned in addition to the raw data calculated by each 
%individual algorithm... Written by Thomas Newell.


%%


chandiff = info.nChannels-size(SPIKES,1);           %Finds the difference between the number of channels and number of channels analyzed by algorithms.
if chandiff > 0                                     %If there is a difference.
    zerotack = zeros(chandiff,size(SPIKES,2));      %Makes a vector of zeros the length of one file.
    SPIKES = [SPIKES;zerotack];                     %Tacks the zeros to the bottom of SPIKES.
end
for r = 1:size(SPIKES,1)                            %For each line of SPIKES...  
    if isnan(SPIKES(r,:))                           %If a line is full of NANs...
        SPIKES(r,:) = 0;                            %Replace the line with zeros...        
    end
end
secondsperindex = info.EndOfFileInSeconds/size(SPIKES,2);       %Conversion factor to turn spike times to seconds.
szseconds = [];                                                 %Makes space for where seizure times (in seconds) will go.
sznumber = 1;                                                   %Starts a seizure counter.
for SpikeChan = 1:size(SPIKES,1)                                %For each channel
    chanmean = mean(SPIKES(SpikeChan,:));
    if chanmean == 0 
        continue; 
    end;            %If the channel is full of zeros, stop and move on.
    Scopy = SPIKES(SpikeChan,:);
    Scopy(find(Scopy > cutoff*mean(Scopy))) = cutoff*mean(Scopy);
    minthresh = absminthresh*mean(Scopy);         %THIS IS WHAT YOU WILL END UP EDITING AND TRYING TO OPTIMIZE OMG OMG ZOMG AHHHGHGHGH!!!                                            %This is the absolute max. Might as well start there...
    [maxval,maxindex] = max(SPIKES(SpikeChan,:)); %Finds max peak in results.
    while maxval > minthresh %while the max peak is above the algorithm threshold:
        szseconds(sznumber,:) = [SpikeChan,maxindex*secondsperindex, maxval];   %Stores seizure [Channel,Time(sec), val], (can't be preallocated).
        lowindex = (maxindex - floor(90/secondsperindex));      %Finds point one minute before peak.
        highindex = (maxindex + floor(90/secondsperindex));     %Finds point one minute after peak.
        if lowindex <= 0                                        %This makes it so you can't zero before the start.
            lowindex = 1;                                       %If the index is negative, make it the start of the results.
        end
        if highindex >= size(SPIKES,2)                          %This makes it so you can't zero after the end.
           highindex = size(SPIKES,2);                          %If the index is too large, make it the end of the file.
        end
        SPIKES(SpikeChan,lowindex:highindex) = 0;               %Actually zeros spike.
        sznumber = sznumber + 1;                                %Counts seizure.
        [maxval,maxindex] = max(SPIKES(SpikeChan,:));           %Finds a next highest peak.
    end
end
DETECTIONS = sznumber - 1;                                      %Displays the number of detections.