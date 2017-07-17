%% Read Scan data
% ScanData = [EweV, ImA]
clear

% Import data from text file.
filename = '/Users/xiaotingzhong/Desktop/SEMs/Polish/Cu_071217_potentiostat/071217_thickPitted_scan4_C01.csv';
delimiter = '\t';
startRow = 2;

% Read columns of data as strings:
formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

% Convert the contents of columns containing numeric strings to numbers.
    % Replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
for col=[1,2,3,4,5,6,7,8,9,10,11]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end
    % Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); 
raw(R) = {NaN}; 

% Allocate imported array to column variable names
mode1 = cell2mat(raw(:, 1));
oxred = cell2mat(raw(:, 2));
error1 = cell2mat(raw(:, 3));
controlchanges = cell2mat(raw(:, 4));
times1 = cell2mat(raw(:, 5));
controlV = cell2mat(raw(:, 6));
EweV = cell2mat(raw(:, 7));
ImA = cell2mat(raw(:, 8));
dqmAh = cell2mat(raw(:, 9));
QQoC = cell2mat(raw(:, 10));
PW = cell2mat(raw(:, 11));

clearvars delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me R;

plot = [EweV,ImA];
[pathstr,name,ext] = fileparts(filename);
v = matlab.lang.makeValidName(name,'ReplacementStyle','delete');
eval([v '= plot;']);
    % Create a file to save the files fullfile(pathstr,[name ext versn])
saveFile = fullfile(pathstr,[name '.mat' '']);
    % When saving variable to file, save the string of variable name
save(saveFile,v);


%% Read polish data
% PolishData=[time, ImA, EweV]
clear

% Import data from text file.
filename = '/Users/xiaotingzhong/Desktop/SEMs/Polish/Cu_071117_potentiostat/071117_2thick_polish2_C01.csv';
delimiter = '\t';
startRow = 2;

% Read columns of data as strings:
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

% Convert the contents of columns containing numeric strings to numbers.
    % Replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
for col=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end

% Allocate imported array to column variable names
mode1 = cell2mat(raw(:, 1));
oxred = cell2mat(raw(:, 2));
error1 = cell2mat(raw(:, 3));
controlchanges = cell2mat(raw(:, 4));
Nschanges = cell2mat(raw(:, 5));
counterinc = cell2mat(raw(:, 6));
Ns = cell2mat(raw(:, 7));
times1 = cell2mat(raw(:, 8));
controlV = cell2mat(raw(:, 9));
EweV = cell2mat(raw(:, 10));
ImA = cell2mat(raw(:, 11));
dqmAh = cell2mat(raw(:, 12));
QQomAh = cell2mat(raw(:, 13));
dQC = cell2mat(raw(:, 14));
QQoC = cell2mat(raw(:, 15));
QchargedischargemAh = cell2mat(raw(:, 16));
halfcycle = cell2mat(raw(:, 17));
EnergychargeWh = cell2mat(raw(:, 18));
EnergydischargeWh = cell2mat(raw(:, 19));
CapacitancechargeF = cell2mat(raw(:, 20));
CapacitancedischargeF = cell2mat(raw(:, 21));
QdischargemAh = cell2mat(raw(:, 22));
QchargemAh = cell2mat(raw(:, 23));
CapacitymAh = cell2mat(raw(:, 24));
Efficiency = cell2mat(raw(:, 25));
cyclenumber = cell2mat(raw(:, 26));
PW = cell2mat(raw(:, 27));

clearvars delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me;

plot = [times1,ImA,EweV];
[pathstr,name,ext] = fileparts(filename);
v = matlab.lang.makeValidName(name,'ReplacementStyle','delete');
eval([v '= plot;']);
    % Create a file to save the files fullfile(pathstr,[name ext versn])
saveFile = fullfile(pathstr,[name '.mat' '']);
    % When saving variable to file, save the string of variable name
save(saveFile,v);

