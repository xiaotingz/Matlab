function en = set110_oneGB(geom110,pars)
% en = set110(geom110,pars)
%
% Dimensionless contribution to energy from <110> rotations
% Very similar to set100; see comments therein for general information.
% Comments in this file will be limited to 110-specific information.


    pwr1 = pars(20); % 110 tilt/twist mix power law:  Twist
    pwr2 = pars(21); % 110 tilt/twist mix power law:  Tilt

    ksi = geom110(1);
    eta = geom110(2);
    phi = geom110(3);

    %
    entwist = twists110(ksi,pars) ;
    entilt  = atgbs110(eta,ksi,pars) ;
    x = phi/(pi/2);
    en = entwist.*(1-x).^pwr1 + entilt.*x.^pwr2 ;
end

function en = twists110(th,pars)
% en = twists110(th,pars)
%
% See comments on set110.

    th1 = pars(22); % 110 twist peak position

    en1 = pars(23); % 110 twist energy peak value
    en2 = pars(24); % Sigma3 energy (110 twist, so not a coherent twin)
    en3 = pars(25); % energy at the symmetry point

    a01 = 0.5;
    a12 = 0.5;
    a23 = 0.5;
    %
    th2 = acos(1/3) ; % Sigma3
    th3 = pi/2 ;  % 110 90-degree boundary is semi-special, although not a CSL

    perio = pi ;  % the twist period
    %
    th = mod(abs(th),perio) ;      % rotation symmetry
    th(th > perio/2) = perio - th(th > perio/2) ;

    en = zeros(size(th));
    %
    select = th <= th1;
    en(select) = en1*rsw(th(select),0,th1,a01) ;

    select = th > th1 & th <= th2;
    en(select) = en2 + (en1-en2)*rsw(th(select),th2,th1,a12) ;

    select = th > th2;
    en(select) = en3 + (en2-en3)*rsw(th(select),th3,th2,a23) ;

end    

function en = stgbs110(th,pars)
% en = stgbs110(th,pars)
% See comments on set110.

    en2 = pars(27); % peak between Sigma1 and Sigma3
    en3 = pars(28); % Coherent Sigma3 twin relative energy; one of the more important element-dependent parameters
    en4 = pars(29); % energy peak between Sigma3 and Sigma11
    en5 = pars(30); % Sigma11 energy
    en6 = pars(31); % energy peak between Sigma11 and Sigma1

    th2 = pars(32); % peak between Sigma1 and Sigma3
    th4 = pars(33); % peak between Sigma3 and Sigma11
    th6 = pars(34); % peak between Sigma11 and higher Sigma1

    a12 = 0.5;
    a23 = 0.5;
    a34 = 0.5;
    a45 = 0.5;
    a56 = 0.5;
    a67 = 0.5;
    %
    % 
    en1 = 0 ;
    en7 = 0 ;

    th1 = 0 ;
    th3 = acos(1/3) ;  % Sigma3
    th5 = acos(-7/11) ; % Sigma11
    th7 = pi ;
    %
    th = pi - th ;  % This is a legacy of an earlier (ksi,eta) mapping 
    %
    en = zeros(size(th));

    select = th<=th2;
    en(select) = en1 + (en2-en1).*rsw(th(select),th1,th2,a12) ;

    select = th >= th2 & th <= th3;
    en(select) = en3 + (en2-en3).*rsw(th(select),th3,th2,a23) ;

    select = th >= th3 & th <= th4;
    en(select) = en3 + (en4-en3).*rsw(th(select),th3,th4,a34) ;

    select = th >= th4 & th <= th5;
    en(select) = en5 + (en4-en5).*rsw(th(select),th5,th4,a45) ;

    select = th >= th5 & th <= th6;
    en(select) = en5 + (en6-en5).*rsw(th(select),th5,th6,a56) ;

    select = th >= th6 & th <= th7;
    en(select) = en7 + (en6-en7).*rsw(th(select),th7,th6,a67) ;


end

function en = atgbs110(eta,ksi,pars)
% en = atgbs110(eta,ksi,pars)
% See comments on set110.

    a = pars(26); % 110 atgb interpolation rsw shape factor
    %
    period = pi ;
    en1 = stgbs110(ksi,pars) ;
    en2 = stgbs110(period-ksi,pars) ;

    en = zeros(size(eta));

    % Power-law interpolation did not work well in this case.  Did an rsw
    % function instead.
    select = en1>=en2;

    en(select) = en2(select) + (en1(select)-en2(select)).*rsw(eta(select),pi,0,a);
    en(~select) = en1(~select) + (en2(~select)-en1(~select)).*rsw(eta(~select),0,pi,a);

end

function en = rsw(theta,theta1,theta2,a)
% en = rsw(theta,theta1,theta2,a)
%
% This function computes the value of Read-Shockley-Wolf function at theta.
% The RSW function is normalized to be 1.0 at theta2 and 0.0 at theta1.
% 
% theta             angle at which to compute the function
% theta1            the starting angle of the interval
% theta2            the end angle of the interval
% a                 parameter defining the shape of the RSW function
%
    dtheta = theta2 - theta1 ;      % Interval of angles where defined
    theta = (theta-theta1)./dtheta*pi/2 ;    % Normalized angle
    % The rest is the RSW function evaluation
    sins = sin(theta) ;
    xlogx = zeros(size(sins));

    % Cut off at small sins to avoid 0*infinity problem.  The proper limit is 0.
    select = sins >= 0.000001;
    xlogx(select) = sins(select).*log(sins(select));

    en = sins - a*xlogx ;
end