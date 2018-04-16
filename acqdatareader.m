function selected_data = acqdatareader(structuretoload,t_start,length)

structure = structuretoload;
fseek(structure.FID,structure.DataStart + (4*t_start*structure.SampleRate*structure.nChannels),'bof');
selected_data.data = fread(structure.FID, [structure.nChannels (length*structure.SampleRate)], 'int32');
selected_data.time = t_start:1/structure.SampleRate:length + t_start - (1/structure.SampleRate);

end

