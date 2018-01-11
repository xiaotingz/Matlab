% in this script, "Grain Curvature is calculated as curvature*area, and then normalized by /R"


% Ausentie ---  ('/Volumes/RESEARCH/Grain Curvature/MinSize=100/setTo0/Jan31_A_setTo0.dream3d');
% Ferrite --- ('/Volumes/RESEARCH/Grain Curvature/MinSize=100/setTo0/Jan31_Fca0.dream3d');
% file = ('/Users/xiaotingzhong/Desktop/Datas/setTo0/Jan31_A_setTo0.dream3d');


%get delete_ID(index for grains to remove) --- remove a outer grain if its centroid is falls within 2<R> from outer frame  
% V4
centroids(1,:) = [];
grain_diameter_raw(1) = [];
num_of_neigh(1) = [];
aveD_raw = sum(grain_diameter_raw) / length(grain_diameter_raw);
%V6
% centroids(1,:) = [];
% grain_diameter_raw(1) = [];
% grain_diameter_raw = grain_diameter_raw.';
% num_of_neigh(1) = [];
% num_of_neigh = num_of_neigh.';
% aveD_raw = sum(grain_diameter_raw) / length(grain_diameter_raw);


total = 0;
delete_id = ones(length(centroids),1);
delete_nlist = ones(length(neighborList),1);

for i = 1:length(centroids)
% check centroid of the grain i 
    if centroids(i,1) < aveD_raw || X - centroids(i,1) < aveD_raw
        delete_id(i) = 0;
    elseif centroids(i,2) < aveD_raw || Y - centroids(i,2) < aveD_raw
        delete_id(i) = 0;
    elseif centroids(i,3) < aveD_raw || Z - centroids(i,3) < aveD_raw
        delete_id(i) = 0;
    end
    
% check centroid of neighbors of the ith grain
    % get neighbor list of the grain i
    nlist_start = sum(num_of_neigh(1:(i-1))) + 1;
    nlist_end = sum(num_of_neigh(1:i));
    % for the jth neighbor grain k, check its centroid
    for j = nlist_start : nlist_end
        k = neighborList(j);
        total = total + k;
        if centroids(k,1) < aveD_raw || X - centroids(k,1) < aveD_raw
            delete_id(i) = 0;
        elseif centroids(k,2) < aveD_raw || Y - centroids(k,2) < aveD_raw
            delete_id(i) = 0;
        elseif centroids(k,3) < aveD_raw || Z - centroids(k,3) < aveD_raw
            delete_id(i) = 0;
        end
    end
        
end
delete_bool = logical(delete_id);

% % get the neighbor list for grains within the final range
% for i = 1 : length(delete_id)
%     if delete_id(i) == 0;
%         nlist_start = sum(num_of_neigh(1:(i-1))) + 1;
%         nlist_end = sum(num_of_neigh(1:i));
%         
%         for j = nlist_start : nlist_end
%             delete_nlist(j) = 0;
%         end
%     end
% end
% nlist_bool = logical(delete_nlist);


ID_list = (1:length(grain_diameter_raw)).';
ID_ForCal = ID_list(delete_bool);
D_ForCal = grain_diameter_raw(delete_bool);

F_mF_diff = zeros(size(ID_ForCal));
for i = 1 : length(ID_ForCal)
    grainID = ID_ForCal(i);
    total = 0;
    nlist_start = sum(num_of_neigh(1:(grainID-1))) + 1;
    nlist_end = sum(num_of_neigh(1:grainID));
        for j = nlist_start : nlist_end
            k = neighborList(j);
            total = total + num_of_neigh(k);
        end
    F_mF_diff(i) = num_of_neigh(grainID) - total/num_of_neigh(grainID);
end

%%
grain_diameter_raw = grain_diameter_raw.';
num_of_neigh = num_of_neigh.';
D_ForCal = grain_diameter_raw(delete_bool);
NFaces_ForCal = num_of_neigh(delete_bool);
grain_ForCal = [ID_ForCal,D_ForCal,NFaces_ForCal];

% calculate grain curvature
s_curv = 0;
tmp4 = zeros(length(grain_ForCal),2);

for i = 1:length(tmp4)
    for j = 1:length(data_final)
        if data_final(j,1) == grain_ForCal(i,1)
            s_curv = s_curv - data_final(j,3)*data_final(j,4);
        elseif data_final(j,2) == grain_ForCal(i,1)
            s_curv = s_curv + data_final(j,3)*data_final(j,4);
        end
    end
    
    tmp4(i,1) = grain_ForCal(i,1);
    tmp4(i,2) = s_curv;
    
    s_curv = 0;
end

% normalized grain curvature
grainCurv = tmp4(:,2) ./ (D_ForCal.*((4/3*pi)^(1/3)/2));


curv_FmF = [F_mF_diff, grainCurv];
scatter(curv_FmF(:,1),curv_FmF(:,2),10,[0.5,0.5,0.5],'filled');
hold on




start = -40;
width = 2;
step = 40;
data_grid = zeros(step,2);
cnt = 0;
total = 0;

for i = 1:step
% mean Curv for each Bin
    for j = 1:length(curv_FmF)
        if curv_FmF(j,1) >= start && curv_FmF(j,1) < start+width
            total = total + curv_FmF(j,2);
            cnt = cnt + 1;
        end
    end
    data_grid(i,1) = start + width/2;
    data_grid(i,2) = total/cnt;
    
    start = start + width;
    total = 0;
    cnt = 0;
end

scatter(data_grid(:,1),data_grid(:,2),50,[1,0,0],'s','filled')
% axis([-40,30,-1,2]);
set(gca,'fontsize',13)
xlabel('Number of Faces Difference (F - mF)','FontSize',15);
ylabel('Normalized Grain Curvature','FontSize',15);
box on
% disable tick on top and right of the box
a = gca;set(a,'box','off','color','none');
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
axes(a)
linkaxes([a b])

line([-40,50],[0,0],'LineStyle','--', 'Color','k')
line([0,0],[-30,30],'LineStyle','--', 'Color','k')
