function data_grain = calcGrainCurvature(data_face, grain_ForCal)

s_curv = 0;
s_area = 0;
tmp = zeros(length(grain_ForCal),4);

for i = 1:length(grain_ForCal)
    for j = 1:length(data_face)
        if data_face(j,1) == grain_ForCal(i,1)
            s_area = s_area + data_face(j,3);
            s_curv = s_curv - data_face(j,3)*data_face(j,4);
        elseif data_face(j,2) == grain_ForCal(i,1)
            s_area = s_area + data_face(j,3);
            s_curv = s_curv + data_face(j,3)*data_face(j,4);
        end
    end
    
    tmp(i,1) = grain_ForCal(i,1);
    tmp(i,2) = s_area;
    tmp(i,3) = s_curv;
    tmp(i,4) = s_curv/s_area;
    
    s_curv = 0;
    s_area = 0;
end
data_grain = [grain_ForCal,tmp(:,3)];

end
