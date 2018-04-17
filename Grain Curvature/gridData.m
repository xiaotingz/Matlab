function data_grid = gridData(xAxis, data_grain)

grainCurv = data_grain(:,5);
if strcmp(xAxis, 'D')
    x_data = data_grain(:,2);
    % - Austenite & Ferrite Sample start with 0.9; STO1470 start with 0.8
    start = 0; 
    width = 0.2;
    step = ceil((max(x_data) - start)/width)+1;
elseif strcmp(xAxis, 'numFaces')
    x_data = data_grain(:,3);
    start = 0; 
    width = 1;
    step = ceil((max(x_data) - start)/width)+1;
elseif strcmp(xAxis, 'numEdges')
    x_data = data_grain(:,4);
    start = 0; 
    width = 3;
    step = ceil((max(x_data) - start)/width)+1; 
end

    data_grid = zeros(step,6);
    s_cur = 0;
    s_x = 0;
    cnt = 0;
    s_sqdiff = 0;
    for i = 1:step
    % mean Curv for each Bin
        for j = 1:length(data_grain)
            if x_data(j) >= start && x_data(j) < start+width
                s_cur = s_cur + grainCurv(j);
                s_x = s_x + x_data(j);
                cnt = cnt + 1;
            end
        end
        data_grid(i,1) = start;
        data_grid(i,2) = start + width;
        data_grid(i,3) = s_x/cnt;
        data_grid(i,4) = s_cur/cnt;
        data_grid(i,5) = cnt;
    % standard mean deviation
        for j = 1:length(data_grain)
            if x_data(j)>= start && x_data(j) < start+width
                s_sqdiff = s_sqdiff + (grainCurv(j) - data_grid(i,4))^2;
            end
        end
        data_grid(i,6) = sqrt(s_sqdiff/cnt);

        start = start + width;
        cnt = 0;
        s_cur = 0;
        s_x = 0;
        s_sqdiff = 0;
    end
end



