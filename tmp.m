property = zeros(3^3,3);

size = 3;
i = 1;
cnt3 = 1;
for a = 1 : size
    cnt2 = 1;
    for b = 1 : size
        cnt1 = 1;
        for c = 1 : size
            property(i,1) = cnt1;
            property(i,2) = cnt2;
            property(i,3) = cnt3;
            
            cnt1 = cnt1 + 1;
            i = i + 1;
        end
        cnt2 = cnt2 + 1;
    end
    cnt3 = cnt3 + 1;
end

index1 = zeros(length(property),1);
for i = 1 : length(property)
    index1(i) = property(i,1) + property(i,2)*size + property(i,3)*size^2;
end
index1 = index1 - 13;

index2 = zeros(length(property),1);
for i = 1 : length(property) 
    index2(i) = property(i,1)*size^2 + property(i,2)*size + property(i,3);
end
index2 = index2 - 13;

full = [property,index1,index2];

%point is possibly giving different column the right weight