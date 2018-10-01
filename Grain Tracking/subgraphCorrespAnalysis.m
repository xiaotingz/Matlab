load('180822');

num_pieces = zeros(length(face_piecewise), 2);
for i = 1:length(num_pieces)
    idx = face_piecewise(i);
    obj_facelabel_an4 = tracked_uniqueface_an4(idx, :);
    obj_facelabel_an5 = tracked_uniqueface_an5(idx, :);
    
    mask_objface_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2));
    mask_objface_an5 = (facelabel_an5(:,1) == obj_facelabel_an5(1) & facelabel_an5(:,2) == obj_facelabel_an5(2));

    face_tri_nodeid_an4 = tri_node_an4(mask_objface_an4, :);
    face_unique_nodeid_an4 = unique(face_tri_nodeid_an4);
    face_tri_nodeid_an5 = tri_node_an5(mask_objface_an5, :);
    face_unique_nodeid_an5 = unique(face_tri_nodeid_an5);
    [subgraph_an4, subgraph_an5] = findSubgraph(face_unique_nodeid_an4, face_unique_nodeid_an5, face_tri_nodeid_an4, face_tri_nodeid_an5);
    
    num_pieces(i, 1) = max(subgraph_an4);
    num_pieces(i, 2) = max(subgraph_an5);
end

%%
all_multipiece = face_piecewise(all(num_pieces>1, 2));
tmp = num_pieces(all(num_pieces>1, 2), :);
all_multipiece_asym = all_multipiece(tmp(:,1) ~= tmp(:,2));
all_multipiece_sym = all_multipiece(tmp(:,1) == tmp(:,2));


amp_asym_size = ones(length(all_multipiece_asym), 2);
for i = 1:length(all_multipiece_asym)    
    idx = all_multipiece_asym(i);
    obj_facelabel_an4 = tracked_uniqueface_an4(idx, :);
    obj_facelabel_an5 = tracked_uniqueface_an5(idx, :);
    
    mask_objface_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2));
    mask_objface_an5 = (facelabel_an5(:,1) == obj_facelabel_an5(1) & facelabel_an5(:,2) == obj_facelabel_an5(2));

    face_tri_nodeid_an4 = tri_node_an4(mask_objface_an4, :);
    face_unique_nodeid_an4 = unique(face_tri_nodeid_an4);
    face_tri_nodeid_an5 = tri_node_an5(mask_objface_an5, :);
    face_unique_nodeid_an5 = unique(face_tri_nodeid_an5);
    [subgraph_an4, subgraph_an5] = findSubgraph(face_unique_nodeid_an4, face_unique_nodeid_an5, face_tri_nodeid_an4, face_tri_nodeid_an5);
    
    cnt_subgraph_an4 = histcounts(subgraph_an4);
    cnt_subgraph_an5 = histcounts(subgraph_an5);
    
    amp_asym_size(i, 1) = max(cnt_subgraph_an4);
    amp_asym_size(i, 2) = max(cnt_subgraph_an5);
end

sum(any(amp_asym_size < 10, 2))












