function en = set111_oneGB(geom111,pars)
% en = set111(geom111,pars)
%
% Dimensionless contribution to energy from <111> rotations
% Very similar to set100; see comments therein for general information.
% Comments in this file will be limited to 111-specific information.

    a = pars(35); % linear part of 111 tilt/twist interpolation
    b = a - 1;  % Ensures correct value at x = 1.

    ksi = geom111(1);
    eta = geom111(2);
    phi = geom111(3);

    entwist = twists111(ksi,pars) ;
    entilt  = atgbs111(eta,ksi,pars) ;
    x       = phi/(pi/2); 

    % This one fit well enough with a simple one-parameter parabola that the
    % more complicated power laws in the other sets weren't needed.
    en      = entwist + (entilt - entwist).*(a*x - b*x.^2) ;
end

function en = twists111(theta,pars)
% en = twists111(theta,pars)
%
% See comments on set111

    thd = pars(37); % 111 twist peak position

    enm = pars(38); % 111 twist energy at the peak
    en2 = pars(28); % Coherent sigma3 twin shows up in two distinct places in the code

    a1 = pars(36); % 111 twist rsw shape parameter
    a2 = a1;

    theta(theta > pi/3) = 2*pi/3-theta(theta > pi/3);

    select = (theta<=thd);
    en = zeros(size(theta));
    en(select) = enm*rsw(theta(select),0,thd,a1) ;
    en(~select) = en2 + (enm - en2)*rsw(theta(~select),pi/3,thd,a2);
end

function en = atgbs111(eta,ksi,pars)
% en = atgbs111(eta,ksi,pars)
%
% This function is a fit to the energies of all 111-tilt boundaries

% There's an additional symmetry in 111 atgbs that doesn't exist in 100 or
% 110 atgbs.  This is because a rotation about [111] equal to half the period
% (i.e. 60 degrees) is equivalent to a mirror reflection in the (111)
% plane.  Both are Sigma3 operations.  The same is not true of the
% 45-degree [100] or the 90-degree [110] rotation.
% The following two lines account for this extra symmetry.
    ksi(ksi > pi/3) = 2*pi/3 - ksi(ksi>pi/3);
    eta(eta > pi/3) = 2*pi/3 - eta(eta>pi/3);

    % Below the following value of ksi, we ignore the eta dependence.  This is
    % because there's little evidence that it actually varies.  Above this
    % value, we interpolate on an rsw function that follows the Sigma3 line,
    % which is also a line of symmetry for the function.
    ksim  = pars(39); % 111 atgb ksi break

    enmax = pars(40); % Energy at the peak (ksi == ksim)
    enmin = pars(41); % energy at the minimum (sigma3, eta == 0)
    encnt = pars(42); % energy at the symmetry point (sigma3, eta == pi/3)

    a1    = 0.5;
    a2    = 0.5;

    etascale = pars(43); % eta scaling parameter for 111 atgb rsw function on Sigma3 line
        % This rsw function is unusual in that the change in shape of the
        % function is much better captured by changing the angular scale rather
        % than changing the dimensionless shape factor.

    en = zeros(size(ksi));

    select = (ksi <= ksim);
    en(select)  = enmax*rsw(ksi(select),0,ksim,a1) ;

    % chi is the shape of the function along the sigma3 line.
    chi = enmin + (encnt-enmin)*rsw(eta(~select),0,pi/(2*etascale),0.5);
    en(~select)  = chi   + (enmax - chi).*rsw(ksi(~select),pi/3,ksim,a2) ;
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