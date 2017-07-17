% [times1,ImA,EweV];
plot_tmp = x071117_2thick_polish2_C01;
endTime = ceil(plot_tmp(length(plot_tmp),1));
plot_data = zeros(endTime,3);
times = 1:endTime;
plot_data(:,1) = times.';


timeCnt = 1;
I_tmp = 0;
E_tmp = 0;
cnt_tmp = 0;
for i = 1:length(plot_tmp)
    if plot_tmp(i,1) <= plot_data(timeCnt, 1)
        I_tmp = I_tmp + plot_tmp(i,2);
        E_tmp = E_tmp + plot_tmp(i,2);
        cnt_tmp = cnt_tmp + 1;
    else
        plot_data(timeCnt, 2) = I_tmp/cnt_tmp;
        plot_data(timeCnt, 3) = E_tmp/cnt_tmp;
        timeCnt = timeCnt + 1;
        I_tmp = plot_tmp(i,2);
        E_tmp = plot_tmp(i,3);
        cnt_tmp = 1;
    end
end
plot_data(timeCnt, 2) = I_tmp/cnt_tmp;
plot_data(timeCnt, 3) = E_tmp/cnt_tmp;


[pathstr,name,ext] = fileparts(filename);
v = matlab.lang.makeValidName(name,'ReplacementStyle','delete');
eval([v '= plot_data;']);
    % Create a file to save the files fullfile(pathstr,[name ext versn])
saveFile = fullfile(pathstr,[name '.mat' '']);
    % When saving variable to file, save the string of variable name
save(saveFile,v);