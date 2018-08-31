function intersection = Intersect(list1, list2)

length1 = length(list1) + 1;
length2 = length(list2) + 1;
i = 1;
j = 1;
intersection = [];
while i < length1 && j < length2
    if list1(i) < list2(j)
        i = i + 1;
    elseif list1(i) > list2(j)
        j = j + 1;
    elseif list1(i) == list2(j)
        intersection = vertcat(intersection, list1(i));
        i = i + 1;
        j = j + 1;
    end
end


end
