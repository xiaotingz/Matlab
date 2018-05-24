function en = twist100(ksi,pars)
% en = twist100(ksi,pars)
%
% Dimensionless 100 twist contribution to the energy

    a = pars(10);   % 100 twist maximum energy
    b = pars(10)*pars(11);  % 100 twist rsw shape factor. The unusual split into two parameters is a holdover from an older version.

    perio = pi/2 ;  % the twist period
    ksi = mod(abs(ksi),perio) ;      % rotation symmetry

    ksi(ksi > perio/2) = perio-ksi(ksi>perio/2);

    % Implement an rsw function of ksi
    sins = sin(2*ksi) ;
    xlogx = sins.*log(sins);
    xlogx(isnan(xlogx))=0; % Force the limit to zero as x -> 0.
    en =  a*sins - b*xlogx ;

end