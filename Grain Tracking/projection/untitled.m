X = randi(10,10, 3);
Y = ones(10, 1);


SVMModel = fitcsvm(X, Y,'KernelFunction','linear', 'Standardize',false);
normal = SVMModel.Beta;
bias = SVMModel.Bias;


scatter3(X(:,1),X(:,2), X(:,3), 80, 'filled')
hold on
x_range = [min(X(:,1)), max(X(:,1))];
y_range = [min(X(:,2)), max(X(:,2))];
[xx,yy]=ndgrid(x_range(1) : (x_range(2) - x_range(1))/2 : x_range(2), y_range(1) : (y_range(2) - y_range(1))/2 : y_range(2));
    % --- calculate corresponding z: Ax + By + Cz + D = 0  --- 
z = (-normal(1)*xx - normal(2)*yy - bias)/normal(3);
scale = ceil(max(sum(features(:,1)==1), sum(features(:,1)==2))/100);
scale = min(scale, 6);

% --- plot the SVM seperation plane and plane normal--- 
colors = get(gca,'colororder');
color = colors(1,:);
surf(xx,yy,z, 'FaceColor', color, 'EdgeColor', color, 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);




mdl = fitlm(X, Y);
coeffs = mdl.Coefficients.Estimate;
norm_normal = norm([coeffs(2:3); -1]);
normal = [coeffs(2:3); -1]./norm_normal;
bias = coeffs(1)/norm_normal;
z = (-normal(1)*xx - normal(2)*yy - bias)/normal(3);
scale = ceil(max(sum(features(:,1)==1), sum(features(:,1)==2))/100);
scale = min(scale, 6);

% --- plot the SVM seperation plane and plane normal--- 
colors = get(gca,'colororder');
color = colors(2,:);
surf(xx,yy,z, 'FaceColor', color, 'EdgeColor', color, 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);


daspect([1 1 1])

