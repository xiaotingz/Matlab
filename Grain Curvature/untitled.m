centroids = roundn(h5read('/Volumes/RESEARCH/Grain Curvature/Jul.20 Austenite/processing/Jul20_ForCurv_centroid.dream3d','/VoxelDataContainer/FIELD_DATA/Centroids'),-5).';
max(centroids(:,1))
max(centroids(:,2))
max(centroids(:,3))
facelabel = double(h5read('/Volumes/RESEARCH/Grain Curvature/Jul.20 Austenite/processing/Jul20_ForCurv_centroid.dream3d','/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels'));
grain_diameter_raw = roundn(h5read('/Volumes/RESEARCH/Grain Curvature/Jul.20 Austenite/processing/Jul20_ForCurv_centroid.dream3d'),-5);
