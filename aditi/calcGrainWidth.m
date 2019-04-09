%calculate grain width
function data_face = calcGrainWidth(triangle_data,Edge_data)
% -- calc mean width for each grain
% -- get rid of bad datas
    % 1.facelabel <= 0; 2.extreme curvature value; 3.NaN(sometimes)
%also filtering data (deleting bad grains and abnormal curvature values)
index = (triangle_data(:,1)<=0 | triangle_data(:,2)<=0 | triangle_data(:,10) > 100 | triangle_data(:,10) < -100);
triangle_data(index,:)=[];

%now we calculate the mean width using the formula=(1/2pi) * summation(edgelength*exterior turning angle)
sort_tri = sortrows(triangle_data,1); %sorting the matrix according to the first face label
M = unique(sort_tri(:,1));%listing total num. of grains
R = histc(sort_tri(:,1),unique(sort_tri(:,1)));%count of each face label
R2= histc(sort_tri(:,2),unique(sort_tri(:,2)));%count of each face label in col 2
k=1;
sum=0;
ctr=0;
for i = 1:length(M)
    i
while(sort_tri(k,1)==i | sort_tri(k,2)==i)   % to find the same face label from both column1 and column2
      for j=1:R(i) %running it as many times as face label 'i' exists
          for p=j+1:R(i)
              com = intersect(sort_tri(j,3:5),sort_tri(p,3:5));%vertices of the common edge between 2triangles of a grain
              if length(com) == 2 % if there are 2 common vertices
                  len=Edge_data(find(Edge_data(:,1)==com(:,1) & Edge_data(:,2)==com(:,2)),3);%edge length of the 2 common triangles
                  %now we find the external turning angle:
                  if (sort_tri(j,1)<sort_tri(j,2) & sort_tri(p,1)<sort_tri(1,2)) | (sort_tri(j,1)>sort_tri(j,2) & sort_tri(p,1)>sort_tri(1,2))
                      %check if face label is on the same side for both
                      %triangles
                     angle = 6.28319 - atan2(norm(cross(sort_tri(j,6:8),sort_tri(p,6:8))),dot(sort_tri(j,6:8),sort_tri(p,6:8)));% (360 deg - theta) =exterior turning angle(radians)
                     sum=sum+(len*angle);
                     ctr=ctr+1;
                  else
                      angle = atan2(norm(cross(sort_tri(j,6:8),sort_tri(p,6:8))),dot(sort_tri(j,6:8),sort_tri(p,6:8)));%if face labels are on opposite sides, then no need to subtract form 360deg
                      sum=sum+(len*angle);
                      ctr=ctr+1;
                  end                           
                 
              end
          end
      end
      width(k,1)=sum/(2*pi); %storing the mean grain width for each grain here
      k=k+1;
end
end

               
               
        
end    