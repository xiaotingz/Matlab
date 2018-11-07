% Example of using CaptureFigVid
% Cheers, Dr. Alan Jennings, Research assistant professor, 
% Department of Aeronautics and Astronautics, Air Force Institute of Technology


% %% Set up 3D plot to record
% figure(171);clf;
% surf(peaks,'EdgeColor','none','FaceColor','interp','FaceLighting','phong')
% daspect([1,1,.3]);axis tight;


% ########################################################
% https://www.mathworks.com/matlabcentral/fileexchange/41093-create-video-of-rotating-3d-plot
% ########################################################
% Set up recording parameters (optional), and record
OptionZ.FrameRate=20; OptionZ.Duration=8; OptionZ.Periodic=true;
file_name = ['subgrpaph_2to1_corresp_pair', num2str(idx), ''];
% file_name = 'good_example_pair_6';
CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], file_name ,OptionZ)

% CaptureFigVid([-20,10;-110,10;-190,80], file_name, OptionZ)
% CaptureFigVid([0,10;-182,10], 'test', OptionZ)



