clear all

% generate random points in homochoric space
points = rand(10000,3)*(sqrt(2)-1);

% get the random points in FZ
cond = zeros(length(points),1);
for i = 1:length(points)
    cond_1 = ((points(i,1) >= points(i,2)) && (points(i,2) >= points(i,3)));
    cond_2 = (points(i,1) <= sqrt(2)-1);
    cond_3 = (points(i,1) + points(i,2) + points(i,3) <= 1);
    cond(i) = (cond_1 && cond_2 && cond_3);
end
cond = logical(cond);
tmp1 = points(:,1); 
tmp2 = points(:,2); 
tmp3 = points(:,3);
rodri_v(:,1) = tmp1(cond);
rodri_v(:,2) = tmp2(cond);
rodri_v(:,3) = tmp3(cond);

% convert the points to homochoric space
angle = zeros(length(rodri_v),1);
axis = zeros(length(rodri_v),3);
homo_v = zeros(length(rodri_v),3);
homo_num = zeros(length(rodri_v),1);
for i = 1:length(rodri_v)
    angle(i) = 2*atan(norm(rodri_v(i,:)));
    axis(i,:) = rodri_v(i,:)/norm(rodri_v(i,:));
    homo_num(i) =  (3/4 * (angle(i) - sin(angle(i))))^(1/3);
    homo_v(i,:) = axis(i,:) * homo_num(i);
end

%% sort the vectors in homochoric space
homo = [[1:length(homo_num)].', homo_num, zeros(length(homo_num),1), homo_v];
homo = sortrows(homo, 2);
distance = 0.1;
groups = zeros(1000, 100);
cnt_r = 1;
cnt_c = 1;
start_v = [0,0,0];
for i = 1:length(homo)
    diff = homo(i, 4:6) - start_v;
    if norm(diff) < distance
        groups(cnt_r,cnt_c) = homo(i,1);
        cnt_c = cnt_c + 1;
    else
        start_v = homo(i, 4:6);
        cnt_r = cnt_r + 1;
        cnt_c = 1;
        groups(cnt_r,cnt_c) = homo(i,1);
    end
end

cnt_r = 1;
while groups(cnt_r,1) ~= 0
    cnt_c = 1;
    group_homo(cnt_r,:) = [0,0,0];
    while (groups(cnt_r, cnt_c) ~= 0 && cnt_c < size(groups,2))
        index = groups(cnt_r, cnt_c);
        group_homo(cnt_r,:) = group_homo(cnt_r,:) + homo_v(index,:);
        cnt_c = cnt_c + 1;
    end
    group_homo(cnt_r,:) = group_homo(cnt_r,:)/cnt_c;
    cnt_r = cnt_r + 1;
end

    


%%
% figure(1)
% CS1 = crystalSymmetry('m-3m');
% CS2 = crystalSymmetry('m-3m');
% oR = fundamentalRegion(CS1,CS2,'antipodal');
% plot(oR)


figure(1)
scatter3(homo_v(:,1),homo_v(:,2),homo_v(:,3),'.')
title('Homochoric Space')
xlabel('X')
ylabel('Y')
zlabel('Z')
rotate3d on
% 
% 
% figure(2)
% scatter3(rodri_v(:,1),rodri_v(:,2),rodri_v(:,3),'.')
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% rotate3d on

% scatter3(test(:,1),test(:,2),test(:,3),'.')

figure(2)
scatter3(group_homo(:,1),group_homo(:,2),group_homo(:,3),'.')
title('Group in Homochoric Space')
xlabel('X')
ylabel('Y')
zlabel('Z')
rotate3d on


figure(3)
scatter3(rodri_v(:,1),rodri_v(:,2),rodri_v(:,3),'.')
title('Rodrigues Space')
xlabel('X')
ylabel('Y')
zlabel('Z')
rotate3d on
