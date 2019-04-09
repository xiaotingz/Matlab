function [xlabelName, ylabelName] = setCorrAxis(objName)
    if strcmp(objName, 'energy_GBCD')
%         xlabelName = 'Grain Boundary Energy (a.u.)';
%         ylabelName = 'Relative Grain Boundary Area (MRD)';
        xlabelName = 'Energy Range (r.u.)';
        ylabelName = 'Average Population (MRD)';
    elseif strcmp(objName, 'energy_ln_GBCD')
        xlabelName = 'Energy Range (r.u.)';
        ylabelName = 'Ln(Average Population, MRD)';
    elseif strcmp(objName, 'GBCD_energy')
        xlabelName = 'Average Population (MRD)';
        ylabelName = 'Energy Range (r.u.)';
    elseif strcmp(objName, 'GBCD_GBCurvD')
        xlabelName = 'Average Population (MRD)';
        ylabelName = 'Average Mean Curvature (\mum^{-1})';
    elseif strcmp(objName, 'GBCD_GBCurvD')
        xlabelName = 'Average Mean Curvature (\mum^{-1})';
        ylabelName = 'Average Population (MRD)';
    elseif strcmp(objName, 'GBCurvD_energy')
        xlabelName = 'Average Mean Curvature (\mum^{-1})';
        ylabelName = 'Energy Range (r.u.)';
    elseif strcmp(objName, 'energy_GBCurvD')
        xlabelName = 'Energy Range (r.u.)';
        ylabelName = 'Average Mean Curvature (\mum^{-1})';
    end
end