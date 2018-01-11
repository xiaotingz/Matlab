clear
O = CrysSym; 

angle_twin = deg2rad(60);
axis_twin = [1,1,1]; axis_twin = axis_twin/norm(axis_twin);
g_twin = AAToG(angle_twin, axis_twin);

angle1 = deg2rad(22.9);
axis1 = [-4,16,-5]; axis1 = axis1/norm(axis1);
angle2 = deg2rad(51.6);			
axis2 = [-2,6,9]; axis2 = axis2/norm(axis2);
g1 = AAToG(angle1, axis1);
g2 = AAToG(angle2, axis2);

% apply sigma3 rotation onto one orientation(arbitarily chose g1)
g1_simga3s = zeros(3,3,24*2);
for i = 1:24
    g1_simga3s(:,:,i) = (O(:,:,i)*g_twin)*g1;
end
for i = 1:24
    g1_simga3s(:,:,(i+24)) = g1*(O(:,:,i)*g_twin).';
end

        

for i = 1:48
    disA = pi/2;
    for j = 1:24
        for k = 1:24
            gg1 = O(:,:,j) * g1_simga3s(:,:,i);
            gg2 = O(:,:,k);

            dg = gg1 * gg2.';
            misA = acos(0.5*(trace(dg)-1));
            RFvec = [(dg(2,3)-dg(3,2)),(dg(3,1)-dg(1,3)),(dg(1,2)-dg(2,1))] /(2*sin(misA)) *tan(misA/2);
    %         if ( misA < disA && sum(misAxis > 0) == 3 )        
            axisInFZ = ( (sum(RFvec >= 0) == 3) && (RFvec(1) <= (sqrt(2) - 1)) && ...
                        (RFvec(1) >= RFvec(2)) && (RFvec(2) >= RFvec(3)) && (sum(RFvec) <= 1));
            if ( misA < disA && axisInFZ )
                g1_simga3s(:,:,i) = dg;
                disA = misA;
%                 message = ['misA_1=', num2str(rad2deg(misA)), ' disA_1=', num2str(rad2deg(disA)), ' axis_1=',num2str(RFvec)];
%                 disp(message);
            end
            dg = gg2 * gg1.';
            misA = acos(0.5*(trace(dg)-1));
            RFvec = [(dg(2,3)-dg(3,2)),(dg(3,1)-dg(1,3)),(dg(1,2)-dg(2,1))] /(2*sin(misA)) *tan(misA/2);
    %         if ( misA < disA && sum(misAxis > 0) == 3 )
            axisInFZ = ( (sum(RFvec >= 0) == 3) && (RFvec(1) <= (sqrt(2) - 1)) && ...
                        (RFvec(1) >= RFvec(2)) && (RFvec(2) >= RFvec(3)) && (sum(RFvec) <= 1));
            if ( misA < disA && axisInFZ )
                g1_simga3s(:,:,i) = dg;
                disA = misA;
%                 message = ['misA_1=', num2str(rad2deg(misA)), ' disA_1=', num2str(rad2deg(disA)), ' axis_2=',num2str(RFvec)];
%                 disp(message);
            end       
        end
    end
end
disA = pi/2;
for j = 1:24
    for k = 1:24
        gg1 = O(:,:,j) * g2;
        gg2 = O(:,:,k);

        dg = gg1 * gg2.';
        misA = acos(0.5*(trace(dg)-1));
        RFvec = [(dg(2,3)-dg(3,2)),(dg(3,1)-dg(1,3)),(dg(1,2)-dg(2,1))] /(2*sin(misA)) *tan(misA/2);
        axisInFZ = ( (sum(RFvec >= 0) == 3) && (RFvec(1) <= (sqrt(2) - 1)) && ...
                    (RFvec(1) >= RFvec(2)) && (RFvec(2) >= RFvec(3)) && (sum(RFvec) <= 1));
        if ( misA < disA && axisInFZ )
            g2 = dg;
            disA = misA;
%             message = ['misA_2=', num2str(rad2deg(misA)), ' disA_2=', num2str(rad2deg(disA)), ' axis=',num2str(RFvec)];
%             disp(message);
        end
        dg = gg2 * gg1.';
        misA = acos(0.5*(trace(dg)-1));
        RFvec = [(dg(2,3)-dg(3,2)),(dg(3,1)-dg(1,3)),(dg(1,2)-dg(2,1))] /(2*sin(misA)) *tan(misA/2);
        axisInFZ = ( (sum(RFvec >= 0) == 3) && (RFvec(1) <= (sqrt(2) - 1)) && ...
                    (RFvec(1) >= RFvec(2)) && (RFvec(2) >= RFvec(3)) && (sum(RFvec) <= 1));
        if ( misA < disA && axisInFZ )
            g2 = dg;
            disA = misA;
%             message = ['misA_2=', num2str(rad2deg(misA)), ' disA_2=', num2str(rad2deg(disA)), ' axis=',num2str(RFvec)];
%             disp(message);
        end       
    end
end


minA = 10000;
for i = 1:48
    for k = 1:24
        for l = 1:24
            gg1 = O(:,:,k)*g1_simga3s(:,:,i);
            gg2 = O(:,:,l)*g2;
            delg_sym = gg1*gg2.';
            angle = acosd(0.5*(trace(delg_sym)-1));
            if angle < minA
                minA = angle;
            end
            delg_sym = gg2*gg1.';
            angle = acosd(0.5*(trace(delg_sym)-1));
            if angle < minA
                minA = angle;
            end
        end
    end
end
minA 