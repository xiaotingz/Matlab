file_an4 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixOrigin3_Hsmooth_curvSwapped_forParaview.dream3d';
file_an5 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_cropToAn4_Hsmooth_forParaview.dream3d';


fl_an4 = h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
fl_an5 = h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
track_status_an4 = h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/TrackStatus')';
track_status_an5 = h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/TrackStatus')';
fl_an4 = sort(fl_an4, 2);
fl_an5 = sort(fl_an5, 2);
tracked_uniqueface_an4_full = sort(tracked_uniqueface_an4_full, 2);
tracked_uniqueface_an5_full = sort(tracked_uniqueface_an5_full, 2);
tracked_uniqueface_an4 = sort(tracked_uniqueface_an4, 2);
tracked_uniqueface_an5 = sort(tracked_uniqueface_an5, 2);

mask_an4 = ismember(fl_an4, tracked_uniqueface_an4_full, 'rows');
mask_an5 = ismember(fl_an5, tracked_uniqueface_an5_full, 'rows');
track_status_an4(~mask_an4) = 0;
track_status_an5(~mask_an5) = 0;
track_status_an4(any(fl_an4<0, 2)) = -1;
track_status_an5(any(fl_an5<0, 2)) = -1;
track_status_an4(mask_an4) = 1;
track_status_an5(mask_an5) = 1;


track_status_an4(ismember(fl_an4, tracked_uniqueface_an4, 'rows')) = 2;
track_status_an5(ismember(fl_an5, tracked_uniqueface_an5, 'rows')) = 2;

h5write(file_an4,'/DataContainers/TriangleDataContainer/FaceData/TrackStatus', track_status_an4');
h5write(file_an5,'/DataContainers/TriangleDataContainer/FaceData/TrackStatus', track_status_an5');







