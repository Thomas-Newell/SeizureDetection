function [M3tot,info] = AC_V1(n,window,step,subwinlen)
%This is a fully commented function that runs the latest version of the
%Autocorrelation Algorithm as of 6/22/2015... n (string) is the title of
%the acq file, window(#) is the size of the window, and step(#) is the step
%size. This function requires an assciated text file for each ACQ, but
%future versions will not... Written by Thomas Newell



close('all'); fclose('all');                        %Closes all figures if any are open.
info = acqreader07092013(n);                        %Loads information from .acq file (requires function acqreader).
emptychans = find(info.ChannelNames == 'e','E');    %Gets empty channel names from file (any channel with an 'e' in an animal name).
%% This portion of the code runs the actual Autocorrelation Algorithm.
disp('Starting *NEW* AutoCorrelation Analysis');          %Displays the start of the autocorr analysis.
winsize = window*500;                               %Sampling rate is 500Hz, so window size is window*500 samples.
stepsize = step*500;                                %Step size is 500 samples/sec times user-defined step.
M3chan = [];                                        %Preallocation might be tricky for partial hours.
endvalue = 0;                                       %Endvalue set to 0... Since we havent started yet. Makes sence at end of code.
emp = 0;                                            %Starts counter for empty channels in the file.
M3tot = 0;                                          %Initializes a space to return the algorithm results.
for hour = 0:floor(info.EndOfFileInHours)           %Stepping through the file in hours...
    selected_data = acqdatareader(info,hour*3600,3600.002);         %Loads in 1 hour of data (plus one data point).
    for chan = 1:info.nChannels                                     %Stepping through channels...
        if any(chan == emptychans)                                  %If the channel is one of the empty ones.
            emp = emp + 1;                                          %Adds one to the empty channel counter.
            continue;                                               %Skips to the next channel if the current one is empty.
        end;
        if size(selected_data.data,1) == 0                          %If there is no data in the first channel...
            break                                                   %Break out of this loop (some bug fix).
        end
        x = selected_data.data(chan,:);                             %Assigns hour of EEG from current channel to x (Can be filtered if wanted).
        winvec = 1:stepsize:size(x,2)-winsize;                      %Creates a starting point for each window in that hour according to step size.
        wincount = 1;                                               %Initializes window counter.
        windata = zeros(size(winvec,2),winsize);                    %Allocates space for each window.
        for k = 1:size(winvec,2)                                    %Stepping through each window start...
            windata(wincount,1:winsize+1) = x(1,winvec(wincount):winvec(wincount)+winsize);  %Creates matrix of windows.
            wincount = wincount + 1;                                %Counts to next window...
        end
%         subwinlen = 60;                                             %This 60 corresponds to the number of samples in each subwindow... (60/500 = 120 millisecond window to catch a spike).
        subwinvec = 1:subwinlen:size(windata,2);          %Creates vector of start points for subwindows
        %%
        ACval = zeros(1,size(windata,1));
        for BigWin = 1:size(windata,1) %For each big window
            EEG = windata(BigWin,:);    %Saves the big window EEG to a new variable
            
            subwincount = 1;
            start = 1;
            Smax = zeros(1,size(subwinvec,2)); Smin = zeros(1,size(subwinvec,2));
            while start < size(EEG,2)
                if size(EEG,2)-start > 60
                subwindata = EEG(1,start:start+subwinlen);
                else
                    subwindata = EEG(1,start:size(EEG,2));
                end
                Smax(1,subwincount) = max(subwindata);
                Smin(1,subwincount) = min(subwindata);
                start = start + subwinlen;
                subwincount = subwincount + 1;
            end
            %Now we have the max and min of each subwindow
            count = 1;
            while count <= size(Smax,2)-2
                HV(count) = min([Smax(count),max([Smax(count+1),Smax(count+2)])]);
                LV(count) = max([Smin(count),min([Smin(count+1),Smin(count+2)])]);
                count = count + 1;
            end
            ACval(1,BigWin) = sum(HV-LV);
            
        end
        
        M3chan(chan,1:size(ACval,2)) = ACval;                              %Results for each channel saved in M3chan.
    end
    if chan == info.nChannels                                                   %If on the last channel...
        M3tot(1:size(M3chan,1),hour*endvalue+1:hour*endvalue+size(M3chan,2)) = M3chan; %Stores hour of data in M3tot
        endvalue = size(M3chan,2);                                              %Creates an endvalue so the next hour of data can be tacked on correctly.
        M3chan = [];                                                            %Clears M3chan for next hour. Preallocation might be tricky for partial hours...
    end
end

%All results stored in M3tot!!!
disp('Finished *NEW* AutoCorrelation Analysis Test');                                       %Complete!
