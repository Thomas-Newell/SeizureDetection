function [SPIKETOT,info] = SS_V1(n,window,step,spikethresh)
%This is a fully commented function that runs the latest version of the
%Line Length Algorithm as of 6/22/2015... n (string) is the title of
%the acq file, window(#) is the size of the window, and step(#) is the step
%size. This function requires an assciated text file for each ACQ, but
%future versions will not... Written by Thomas Newell



close('all'); fclose('all');                        %Closes all figures if any are open.\
load('C:\Users\SeizureDetection\Desktop\SpikeFunc.mat');
info = acqreader07092013(n);                        %Loads information from .acq file (requires function acqreader).
emptychans = find(info.ChannelNames == 'e','E');    %Gets empty channel names from file (any channel with an 'e' in an animal name).
%% This portion of the code runs the Slope Analysis Algorithm
disp('Starting Slope Analysis');                                  %Displays the start of the autocorr analysis.
winsize = window*500;                                   %Sampling rate is 500Hz, so window size is window*500 samples.
stepsize = step*500;                                    %Step size is 500 samples/sec times user-defined step.
hoursize = floor((3600-window+step)/step);              %Finds out how many windows are needed for one hour.
emp = 0;                                                %Initiates empty channel counter.
SPIKETOT = 0;                                              %Initializes a space to return the algorithm results.
% spikethresh = 1.0e11;
for H = 0:floor(info.EndOfFileInHours)                  %Stepping through the file in hours...
    selected_data = acqdatareader(info,H*3600,3600);    %Loads in 1 hour of data.
    for channel = 1:info.nChannels                      %Stepping through channels...
        if any(channel == emptychans)                      %If the channel is one of the empty ones.
            emp = emp + 1;                              %Adds one to the empty channel counter.
            continue;                                   %Skips to the next channel if the current one is empty.
        end;
        x = selected_data.data(channel,:);              %Assigns hour of EEG from current channel to x (Can be filtered if wanted).
        s = abs(diff(x)/.002);
        s(end+1) = s(end);
        S = minmaxfilt(s,5,'max','same');
        itterations = floor((length(S)-winsize+stepsize)/stepsize); %Finds the number of windows that will need to be analyzed.
        SPIKECOUNTchan = zeros(1,itterations);
        if itterations < 1                              %Continues with code as long as there is a window to be analysed.
            break
        end
        result = 1;
        for k = 1:itterations                           %Takes chunks that overlap by (step) seconds...
            loc=(k-1)*stepsize;                         %Location of the beginning of the data to be analyzed.
            A = S(loc+1:loc+winsize) - spikethresh;
            crossings = find(A(1:end-1).*A(2:end)<0);
            numspikes = (size(crossings,2)/2)/window;
            [val,ind] = min(abs(f_spike-numspikes));
            if numspikes*(w_spike(ind)) < 0.1
                SPIKECOUNTchan(result) = 0;
            else
                SPIKECOUNTchan(result) = numspikes*(w_spike(ind));
            end
            result = result+1;
        end
        SPIKETOT(channel,(H*hoursize)+1:(H*hoursize)+length(SPIKECOUNTchan)) = SPIKECOUNTchan; %Stores hour of data in SPIKETOT.
        clear('SPIKECOUNTchan');                           %Clears linelength vector for next hour.
    end
end

%ALL RESULTS ARE STORED IN THE MATRIX SPIKETOT!!!
disp('Finished Slope Analysis');                   %Complete!
