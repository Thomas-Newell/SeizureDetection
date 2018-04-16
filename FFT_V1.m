function [FFTTOT,info] = FFT_V1(n,window,step,scale) 
%This is a fully commented function that runs the latest version of the
%FFT Power Algorithm as of 6/22/2015... n (string) is the title of
%the acq file, window(#) is the size of the window, and step(#) is the step
%size. This function requires an assciated text file for each ACQ, but 
%future versions will not... Written by Thomas Newell



close('all'); fclose('all');                        %Closes all figures if any are open.
load('C:\Users\SeizureDetection\Desktop\WeightFunc')                  %Gets appropriate text file to test against.
info = acqreader07092013(n);                        %Loads information from .acq file (requires function acqreader).
emptychans = find(info.ChannelNames == 'e','E');    %Gets empty channel names from file (any channel with an 'e' in an animal name).
%% This portion of the code runs the actual FFT Power Algorithm.
disp('Starting Power Analysis');            %Displays the start of the autocorr analysis.
windowsize = window*500;                    %Sampling rate is 500Hz, so window size is window*500 samples.
stepsize = step*500;                        %Step size is 500 samples/sec times user-defined step.
FFTchan = [];                               %Preallocation might be tricky for partial hours.
emp = 0;                                    %Creates empty channel counter. 
endvalue = 0;                               %Endvalue set to 0... Since we havent started yet. Makes sence at end of code.
FFTTOT = 0;                                 %Initializes a space to return the algorithm results.
for H = 0:floor(info.EndOfFileInHours)      %Stepping through the file in hours...
    selected_data = acqdatareader(info,H*3600,3600);    %Loads in 1 hour of data.
    for chan = 1:info.nChannels                         %Stepping through channels...
        if any(chan == emptychans)                      %If the channel is one of the empty ones.
            emp = emp + 1;                              %Adds one to the empty channel counter.
            continue;                                   %Skips to the next channel if the current one is empty. 
        end;
        if size(selected_data.data,1) == 0              %If there is no data in the first channel...
            break                                       %Break out of this loop (some bug fix).            
        end
        x = selected_data.data(chan,:);                 %Assigns hour of EEG from current channel to x (Can be filtered if wanted).
        itterations = floor((length(x)-windowsize+stepsize)/stepsize);  %Finds the number of windows that will need to be analyzed.
        if itterations < 1                                              %Continues with code as long as there is a window to be analysed.
            break
        end
        blockstarts = 0:step:length(x)/500-window+1;                    %Finds the indices that each window will start at.
        blockends = blockstarts + window;                               %Finds the indices that each window will end at.
        if length(x)/500 == 3600                                        %For COMPLETE hours...
            data = zeros(itterations,windowsize+1);                     %Allocates space for each window of data for the current hour.
            fftdata = zeros(itterations,windowsize+1);                  %Allocates space for the fft of each window.
        else
            data = []; fftdaa = [];                                     %Makes space for INCOMPLETE hour... preallocating may be tricky...
        end
        data(1,1:windowsize+1) = x(1:blockends(1)*500+1);               %Puts first EEG data window into data matrix.
        for itts = 2:itterations                                        %Puts the rest of the EEG data windows into the data matrix.
            data(itts,1:windowsize+1) = x(blockstarts(itts)*500:blockends(itts)*500);   %Actually places data.
        end
        fftdata(1,1:windowsize+1) = x(1:blockends(1)*500+1);            %Puts first EEG data window into fftdata matrix.
        for itts = 1:itterations                                        %For each window... 
            fftdata(itts,1:size(data,2)) = fft(data(itts,:));           %Performs a fast fourier transform.
        end
        fftdata = fliplr(abs(fftdata));                                 %Flips the data left/right.
        fftdata = fftdata(:,1:floor((size(data,2)-1)/2));               %Only uses half the data since the FFT spectrum is symmetrical.
        newWeightFunc = interp1(linspace(0,250,length(WeightFunc)),WeightFunc,linspace(0,250,size(fftdata,2)));
        for row = 1:size(fftdata,1)
            fftdata(row,:) = fftdata(row,:).*(scale.*newWeightFunc);
        end
        fftvals = sum(fftdata,2); fftvals = fftvals';                   %Sums up the spectrum (0-250 Hz) and transposes it to make a list of sums.
        FFTchan(chan,1:size(fftvals,2)) = fftvals';                     %Puts the results into a matrix titled FFTchan.
    end
    if chan == info.nChannels                                           %If on the last channel...
        FFTTOT(1:size(FFTchan,1),H*endvalue+1:H*endvalue+size(FFTchan,2)) = FFTchan; %Stores hour of data in FFTTOT
        endvalue = size(FFTchan,2);                                     %Creates an endvalue so the next hour of data can be tacked on correctly.
        FFTchan = [];                                                   %Clears FFTchan for next hour. Preallocation might be tricky for partial hours...
    end
end
goods = find(~mean(FFTTOT,2) == 0);                                     %Finds channels where the mean is not 0 (good channels).
[a,b] = find(FFTTOT(goods(1),:) == 0);                                  %Finds zeros at the end of the file.
if size(b,2) == 0                                                       %If all channels are good... Do nothing.
else                                                                    %Otherwise...
    FFTTOT = FFTTOT(:,1:min(b)-1);                                      %Eliminates zeros.
end

%All results are stored in FFTTOT!
disp('Finished Power Analysis');                                        %Complete!