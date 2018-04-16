function structure = acqreader07092013(n)


structure.filestarttime = n;
structure.FID = fopen(n);
if structure.FID == -1 
    error('File Open Failed') 
    
end;
structure.nItemHeaderLen = fread(structure.FID,1,'short');
structure.IVersion = fread(structure.FID,1,'long');
structure.IExtItemHeaderLen = fread(structure.FID,1,'long');
structure.nChannels = fread(structure.FID,1,'short');
fseek(structure.FID,16,'bof');structure.dSampleTime = fread(structure.FID,1,'double'); structure.SampleRate = 1000/structure.dSampleTime;
structure.ChannelNames = zeros(structure.nChannels,1);
fseek(structure.FID,structure.IExtItemHeaderLen,'bof');
structure.IChanHeaderLen = fread(structure.FID,1,'long');
fseek(structure.FID,(structure.IExtItemHeaderLen + (structure.IChanHeaderLen * structure.nChannels)),'bof');
structure.nlength = fread(structure.FID,1,'int32');

structure.DataStart = structure.nlength + (4*structure.nChannels) + (structure.IChanHeaderLen * structure.nChannels) + structure.IExtItemHeaderLen;
fseek(structure.FID,0,'eof');
structure.EndOfFileInBits = ftell(structure.FID);
structure.EndOfFileInSeconds = (structure.EndOfFileInBits - structure.DataStart)/structure.nChannels/4/500;
structure.EndOfFileInHours = structure.EndOfFileInSeconds/3600;


fseek(structure.FID,structure.IExtItemHeaderLen,'bof');
for channels = 1:structure.nChannels
    structure.IChanHeaderLen = fread(structure.FID,1,'long');
    structure.nNum = fread(structure.FID,1,'short');
    structure.szCommentText = fread(structure.FID,40,'char'); 
    structure.szCommentText = structure.szCommentText';
    structure.szCommentText = char(structure.szCommentText);
    structure.ChannelNames(channels,1:40) = structure.szCommentText;
    fseek(structure.FID,structure.IExtItemHeaderLen + (channels*structure.IChanHeaderLen),'bof');
end

clear('structure.szCommentText');
structure.ChannelNames = char(structure.ChannelNames);
% fclose(structure.FID);
end