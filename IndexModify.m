NumStep = 9;
index = zeros(NumStep^5,5);

% list all the bins as increasing Phi(e),theta(d),phi_2(c),PHI(b),phi_1(a)
% in the mean time, 
% (inner loop is excuted first.)
n1 = NumStep;
n1n2 = NumStep^2;
n1n2n3 = NumStep^3;
n1n2n3n4 = NumStep^4;
i = 1;
% NOTICE, Phi is from 0 to 2Pi,so the size of e is 4*NumStep 
for a = 1 : NumStep
    for b = 1 : NumStep
        for c = 1 : NumStep
            for d = 1 : NumStep
                for e  = 1 : NumStep*4
                    index(i,1) = a;
                    index(i,2) = b;
                    index(i,3) = c;
                    index(i,4) = d;
                    index(i,5) = e;
                    
                    gbcd_index(i) = a + b*n1 + c*n1n2 + d*n1n2n3 + e*n1n2n3n4;
                    i = i + 1;
                end
            end
        end
    end
end

gbcd_index = (gbcd_index - 7380).';

% index1 = zeros(length(property),1);
% for i = 1 : length(property)
%     index1(i) = property(i,1) + property(i,2)*NumStep + property(i,3)*NumStep^2 + property(i,4)*NumStep^3 + property(i,5)*NumStep^4;
% end
% index1 = index1 - 7380;

% % can't just give different digit different weight, because of Phi
% index2 = zeros(length(property),1);
% for i = 1 : length(property) 
%     index2(i) = property(i,1)*NumStep^4 + property(i,2)*NumStep^3 + property(i,3)*NumStep^2 + property(i,4)*NumStep + property(i,5);
% end
% index2 = index2 - 7380;

% full = [property,index1,index2];
