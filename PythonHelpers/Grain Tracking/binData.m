function data_grid = binData(x, stepsize)
% ------------------ load data for debug --------------------
%     x = log(D_An4_TG/(sum(D_An4_TG)/length(D_An4_TG)));
%     x = numNeigh_An4_TCG;
%     stepsize = 2;
% -----------------------------------------------------------
    tmp = min(x):stepsize:max(x);
    data_grid = [tmp', zeros(length(tmp), 1)];
    for i = 1:length(data_grid)-1
        lowLim = data_grid(i);
        upLim = data_grid(i + 1);
        if i < length(data_grid)-1
            mask2 = (x >= lowLim & x < upLim);
        elseif i == length(data_grid)-1
            mask2 = (x >= lowLim & x <= upLim);
        end
        data_grid(i,2) = sum(mask2)/length(x);
    end
    data_grid = data_grid(data_grid(:,2)~=0,:);
end
    