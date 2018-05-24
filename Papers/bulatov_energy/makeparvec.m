function [par43,AlCuparameter] = makeparvec(AlCuparameter,eRGB,par42Al,par42Cu)
% [par43,AlCuparameter] = makeparvec(AlCuparameter,eRGB,par42Al,par42Cu)
%
% Creates a 43-parameter vector as used by weightedmeanenergy
%
% Arguments are:
% AlCuparameter:  Position on the Al-Cu axis, where 0 is Al and 1 is Cu
% (this parameter is capital Phi in the journal article).
% This is related to eSF/eRGB, where eSF is the stacking-fault energy.
% Optionally, AlCu parameter is a string 'Al', 'Ni', 'Au', or 'Cu', which
% then sets all other parameters automatically.  You can call it with
% just this parameter if you wish.
%
% eRGB:  Energy of a "random" grain boundary in J/m^2
%
% There are some additional options that are still written into the
% function but are currently not used:
% par42Al:  The 42 dimensionless parameters for Al
%
% par42Cu:  The 42 dimensionless parameters for Cu
%
% Note a majority of the entries for par42Al and par42Cu are normally
% identical.
%
% All parameters have default values.  Defaults for par42Al and par42Cu are
% the values found by numerical fitting to the 388*4 boundaries.  
% eRGB and AlCuparameter default to the values for Cu.
%
% Optionally returns the numerical AlCuparameter so the user can read the
% default value for each element.

    if ~exist('eRGB','var') || isempty(eRGB),
        eRGB = 1.03669431227427;    % Value for Cu
    end

    if ~exist('AlCuparameter','var') || isempty(AlCuparameter),
        AlCuparameter = 1;  % Value for Cu
    end

    if ~exist('par42Al','var') || isempty(par42Al),
        par42Al = [0.405204179289160;0.738862004021890;0.351631012630026;2.40065811939667;1.34694439281655;0.352260396651516;0.602137375062785;1.58082498976078;0.596442399566661;1.30981422643602;3.21443408257354;0.893016409093743;0.835332505166333;0.933176738717594;0.896076948651935;0.775053293192055;0.391719619979054;0.782601780600192;0.678572601273508;1.14716256515278;0.529386201144101;0.909044736601838;0.664018011430602;0.597206897283586;0.200371750006251;0.826325891814124;0.111228512469435;0.664039563157148;0.241537262980083;0.736315075146365;0.514591177241156;1.73804335876546;3.04687038671309;1.48989831680317;0.664965104218438;0.495035051289975;0.495402996460658;0.468878130180681;0.836548944799803;0.619285521065571;0.844685390948170;1.02295427618256];
    end

    if ~exist('par42Cu','var') || isempty(par42Cu),
        par42Cu = [0.405204179289160;0.738862004021890;0.351631012630026;2.40065811939667;1.34694439281655;3.37892632736175;0.602137375062785;1.58082498976078;0.710489498577995;0.737834049784765;3.21443408257354;0.893016409093743;0.835332505166333;0.933176738717594;0.896076948651935;0.775053293192055;0.509781056492307;0.782601780600192;0.762160812499734;1.10473084066580;0.529386201144101;0.909044736601838;0.664018011430602;0.597206897283586;0.200371750006251;0.826325891814124;0.0226010533470218;0.664039563157148;0.297920289861751;0.666383447163744;0.514591177241156;1.73804335876546;2.69805148576400;1.95956771207484;0.948894352912787;0.495035051289975;0.301975031994664;0.574050577702240;0.836548944799803;0.619285521065571;0.844685390948170;0.0491040633104212];
    end

    if ischar(AlCuparameter),
        switch AlCuparameter
            case 'Ni'
                eRGB = 1.44532834613925;
                AlCuparameter = 0.767911805073948;
            case 'Al'
                eRGB = 0.547128733614891;
                AlCuparameter = 0;
            case 'Au'
                eRGB = 0.529912885175204;
                AlCuparameter = 0.784289766313152;
            case 'Cu'
                eRGB = 1.03669431227427;
                AlCuparameter = 1;
            otherwise
                error('Undefined element')
        end
    end

    par43 = [eRGB;(par42Al+AlCuparameter*(par42Cu-par42Al))];        
end
