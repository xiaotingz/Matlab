function g = EAtoG(EA)
%     """
%     Input: a set of Euler Angle
%                 size=(3,)
%     Output: the corresponding orientation matrix g
%                 size = (3, 3)
%     """
    g = zeros(3,3);
    
    g(1,1)=cosd(EA(1))*cosd(EA(3))-sind(EA(1))*sind(EA(3))*cosd(EA(2));
    g(1,2)=sind(EA(1))*cosd(EA(3))+cosd(EA(1))*sind(EA(3))*cosd(EA(2));
    g(1,3)=sind(EA(3))*sind(EA(2));
    g(2,1)=-cosd(EA(1))*sind(EA(3))-sind(EA(1))*cosd(EA(3))*cosd(EA(2));
    g(2,2)=-sind(EA(1))*sind(EA(3))+cosd(EA(1))*cosd(EA(3))*cosd(EA(2));
    g(2,3)=cosd(EA(3))*sind(EA(2));
    g(3,1)=sind(EA(1))*sind(EA(2));
    g(3,2)=-cosd(EA(1))*sind(EA(2));
    g(3,3)=cosd(EA(2));
end