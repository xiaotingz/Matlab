function [xlabelName, ylabelName] = setCorrAxis(objName)
    if strcmp(objName, 'energy_GBCD')
        xlabelName = 'Grain Boundary Energy (a.u.)';
        ylabelName = 'Relative Grain Boundary Area (MRD)';
    elseif strcmp(objName, 'GBCD_energy')
        xlabelName = 'Relative Grain Boundary Area (MRD)';
        ylabelName = 'Grain Boundary Energy (a.u.)';
    elseif strcmp(objName, 'GBCD_GBCurvD')
        xlabelName = 'Relative Grain Boundary Area (MRD)';
        ylabelName = 'Grain Boundary Curvature (\mum^{-1})';
    elseif strcmp(objName, 'GBCD_GBCurvD')
        xlabelName = 'Grain Boundary Curvature (\mum^{-1})';
        ylabelName = 'Relative Grain Boundary Area (MRD)';
    elseif strcmp(objName, 'GBCurvD_energy')
        xlabelName = 'Grain Boundary Curvature (\mum^{-1})';
        ylabelName = 'Grain Boundary Energy (a.u.)';
    elseif strcmp(objName, 'energy_GBCurvD')
        xlabelName = 'Grain Boundary Energy (a.u.)';
        ylabelName = 'Grain Boundary Curvature (\mum^{-1})';
    end
end