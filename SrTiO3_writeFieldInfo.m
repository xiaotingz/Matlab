orientation = textread('r_redone.txt');
FieldID = [1:length(orientation)].';
PhaseID = ones(length(orientation),1);

% the way fprintf works is writing the first row into the first column
data = [FieldID, PhaseID, orientation(:,1:3)].';

% NOTICE! the first line of orientation file should be #fields. Do it directly on the txt file
fileID = fopen('SrTiO3_orientation.txt','w');
fprintf(fileID,'%-4i %-4i %-15.8f %-15.8f %-15.8f\n', data);
fclose(fileID);
