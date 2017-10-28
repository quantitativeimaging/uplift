function [time,value] = ReadDataLog(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   File name:    ReadDatalog.m
%
%   Purpose  :    Convert datalogger output to time vector and value
%                 vector.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_id = fopen(filename);
% read file line by line
line= fgetl(file_id);
lines = {};
while ischar(line)
    lines{end+1} = line;
    line = fgetl(file_id);
end
fclose(file_id);

% Reject data if
% 1. the line contains alphabets (beginning, end, interruption messages)
% 2. the length of line is wrong (odd datum point, sometimes the 1st)
lines_lengths = [];
for row = 1:length(lines)
    % Rejection 1
    % Checks whether every character in a line is alphabetic
    if sum(isstrprop(lines{row},'alpha')) == 0
        lines_lengths(row) = length(lines{row});
    else %reject line if yes
        lines_lengths(row) = -1;
    end
end

% find the mode of length of line, this is assummed to be the main data
% format
line_length = mode(lines_lengths);
raw_data = {};
for row = 1:length(lines_lengths)
    % Rejection 2
    % Only copy lines consistent with length of the bulk data
    if lines_lengths(row) == line_length
        raw_data{end+1} = lines{row};
    end
end

% read time and values for each row

time = [];
value = [];

for row = 1:length(raw_data)
   % Time and Date is stored in the first 12 characters of a line
   time_cell = textscan(raw_data{row}(1:12),'%f:%f:%f:%f');
   % convert time to seconds; should be fine unless data is collected
   % during mid-night
   time(end+1) = time_cell{1}*3600 + time_cell{2}*60 + time_cell{3} + time_cell{4}/1000;
   % 14 char in front of the logged value
   value(end+1) = str2num(raw_data{row}(14:end));    
   
end

end