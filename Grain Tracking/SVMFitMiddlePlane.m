load('/Users/xiaotingzhong/Dropbox/grainTracking_forCluster/180822_FaceCorresp');
load('/Users/xiaotingzhong/Dropbox/grainTracking_forCluster/180822');
idx = 1329;

x_to_y = X_to_Y{idx};
obj_facelabel_an4 = tracked_uniqueface_an4(idx, :);
obj_facelabel_an5 = tracked_uniqueface_an5(idx, :);
face_node_info = getSingleFaceNodes(file_an4, obj_facelabel_an4, file_an5, obj_facelabel_an5);

visualizeFace(face_node_info, x_to_y)

X = [face_node_info{4,1}; face_node_info{4,2}];
Y = [ones(length(face_node_info{4,1}), 1); 2*ones(length(face_node_info{4,2}), 1)];
SVMModel = fitcsvm(X, Y,'KernelFunction','linear',...
    'Standardize',false,'ClassNames',{'1','2'});

x_range = [min(X(:,1)), max(X(:,1))];
y_range = [min(X(:,2)), max(X(:,2))];
normal = SVMModel.Beta;
bias = SVMModel.Bias;
plotPlane(normal, bias, x_range, y_range)

Y_pre = str2double(cell2mat(predict(SVMModel,X)));
error = sum(Y ~= Y_pre)
%%
rng(1); % For reproducibility
r = sqrt(rand(100,1)); % Radius
t = 2*pi*rand(100,1);  % Angle
data1 = [r.*cos(t), r.*sin(t)]; % Points

r2 = sqrt(3*rand(100,1)+1); % Radius
t2 = 2*pi*rand(100,1);      % Angle
data2 = [r2.*cos(t2), r2.*sin(t2)]; % points

data3 = [data1;data2];
theclass = ones(200,1);
theclass(1:100) = -1;

%Train the SVM Classifier
cl = fitcsvm(data3,theclass,'KernelFunction','rbf',...
    'BoxConstraint',Inf,'ClassNames',[-1,1]);

% Predict scores over the grid
d = 0.02;
[x1Grid,x2Grid] = meshgrid(min(data3(:,1)):d:max(data3(:,1)),...
    min(data3(:,2)):d:max(data3(:,2)));
xGrid = [x1Grid(:),x2Grid(:)];
[~,scores] = predict(cl,xGrid);

% Plot the data and the decision boundary
figure;
h(1:2) = gscatter(data3(:,1),data3(:,2),theclass,'rb','.');
hold on
ezpolar(@(x)1);
h(3) = plot(data3(cl.IsSupportVector,1),data3(cl.IsSupportVector,2),'ko');
contour(x1Grid,x2Grid,reshape(scores(:,2),size(x1Grid)),[0 0],'k');
legend(h,{'-1','+1','Support Vectors'});
axis equal
hold off


cl2 = fitcsvm(data3,theclass,'KernelFunction','rbf');
[~,scores2] = predict(cl2,xGrid);

figure;
h(1:2) = gscatter(data3(:,1),data3(:,2),theclass,'rb','.');
hold on
ezpolar(@(x)1);
h(3) = plot(data3(cl2.IsSupportVector,1),data3(cl2.IsSupportVector,2),'ko');
contour(x1Grid,x2Grid,reshape(scores2(:,2),size(x1Grid)),[0 0],'k');
legend(h,{'-1','+1','Support Vectors'});
axis equal
hold off



%%
class_1 = [0,0,0; 0,1,0; 1,0,0; 1,1,0];
class_2 = class_1;
class_2(:,3) = class_2(:,3) + 1;
% class_1 = [0,0,0; 1,0,0; 0,0,1; 1,0,1];
% class_2 = class_1;
% class_2(:,2) = class_2(:,2) + 1;
% class_1 = [0,0,0; 0,1,0; 0,0,1; 0,1,1];
% class_2 = class_1;
% class_2(:,1) = class_2(:,1) + 1;

X = [class_1; class_2];
Y = [ones(4,1); 2*ones(4,1)];

SVMModel = fitcsvm(X,Y,'KernelFunction','linear', 'Standardize',false);

figure
h = scatter3(X(:,1), X(:,2), X(:,3));
hold on
h(3) = plot3(X(cl.IsSupportVector,1),X(cl.IsSupportVector,2), X(cl.IsSupportVector,3),'ko');




