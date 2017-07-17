function euler = GToE(delg)
    euler = zeros(1,3);
    eps = 1.e-4;
    
    if abs(delg(3,3) - 1) < eps
        euler(2) = 0;
        euler(1) = atan2(delg(1,2),delg(1,1))/2;
        euler(3) = euler(1);
    else
        euler(2) = acos(delg(3,3));
        euler(1) = atan2(delg(3,1)/sin(euler(2)), -delg(3,2)/sin(euler(2)));
        euler(3) = atan2(delg(1,3)/sin(euler(2)), delg(2,3)/sin(euler(2)));
    end
    
    if euler(1) < 0
        euler(1) = 2*pi + euler(1);
    elseif euler(3) < 0
        euler(3) = 2*pi + euler(3);
    end 
end

% from Dr. Rohrer's sub.f
%     eps = 1.e-6;
%     euler = zeros(1,3);
%     delg = eye(3);
%     
% 	xz=delg(3,3);
% 	
% 	if (xz*xz < 1.0)
%         sf=sqrt(1.0-xz*xz);
%     else
%         sf=0.0;
%     end
% 
% 	if sf > eps
%         euler(1) =-delg(3,2)/sf;
%         euler(1) = acos(euler(1));
%     end
%     
% 	if (delg(3,1) < 0.0)
%         euler(1)=2.0*pi-euler(1);
%         euler(2)=acos(xz);
%         euler(3)=delg(2,3)/sf;
%         euler(3)=acos(euler(3));
%     end
%     
% 	if (delg(1,3) < 0.0)
%         euler(3)=2.0*pi-euler(3);
%     else
%         euler(1)=delg(2,2);
%         euler(1)=acos(euler(1));
%     end
%     
% 	if (delg(1,2) < 0.0)
%         euler(1)=2.0*pi-euler(1);
%         euler(2)=0.0;
%         euler(3)=0.0;
%     end

    
