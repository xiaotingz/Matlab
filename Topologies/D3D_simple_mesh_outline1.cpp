  // Quick Mesh key procedures

  // first determining which nodes are actually boundary nodes and
  // count number of nodes and triangles that will be created
  for(int64_t k = 0; k < zP; k++)
  {
    for(int64_t j = 0; j < yP; j++)
    {
      for(int64_t i = 0; i < xP; i++)
      {
        point = (k * xP * yP) + (j * xP) + i;
        neigh1 = point + 1;
        neigh2 = point + xP;
        neigh3 = point + (xP * yP);

        if i = 0 | i = xP-1 | j ...
          nodeId_1,2,3,4 = XXX 
          nodeCount ++ 
          triangleCount ++
          triangleCount ++

 // now create node and triangle arrays knowing the number that will be needed
   for(int64_t k = 0; k < zP; k++)
  {
    for(int64_t j = 0; j < yP; j++)
    {
      for(int64_t i = 0; i < xP; i++)
      {
        point = (k * xP * yP) + (j * xP) + i;
        neigh1 = point + 1;
        neigh2 = point + xP;
        neigh3 = point + (xP * yP);

        // every face is defined by one neighboring grain and 4 points, among which the point 2 and point 3 are shared by the two triangels in this face
        // loop all voxels are there are only 3 cases in which mesh should be generated：
	        // 1. The first voxel, on free surface (0)
	        // 2. The last voxel, on free surface  (end - 1)
	        // 3. Two voxels with different GrainID meet
        if i = 0
          nodeId1 = (k * (xP + 1) * (yP + 1)) + (j * (xP + 1)) + i;
          nodeId2 = (k * (xP + 1) * (yP + 1)) + ((j + 1) * (xP + 1)) + i;
          nodeId3 = ((k + 1) * (xP + 1) * (yP + 1)) + (j * (xP + 1)) + i;
          nodeId4 = ((k + 1) * (xP + 1) * (yP + 1)) + ((j + 1) * (xP + 1)) + i;
          
          triangle[triangleIndex * 3 + 0-2] = m_NodeIds[nodeId 123];
          triangleIndex++;
          triangle[triangleIndex * 3 + 0-2] = m_NodeIds[nodeId 243];
          triangleIndex++;

          ownerLists[m_NodeIds[nodeId 1-4]].insert(m_FeatureId(point))
          ownerLists[m_NodeIds[nodeId 1-4]].insert(-1);

        if j = 0
        if k = 0
        if(i == (xP - 1))
          ownerLists[m_NodeIds[nodeId 1-4]].insert(m_FeatureId(point))
          ownerLists[m_NodeIds[nodeId 1-4]].insert(-1);
        else if(m_FeatureIds[point] != m_FeatureIds[neigh1])
          ownerLists[m_NodeIds[nodeId 1-4]].insert(m_FeatureIds[point]);
          ownerLists[m_NodeIds[nodeId 1-4]].insert(m_FeatureIds[neigh1]);
        if(j == (yP - 1))
        else if(m_FeatureIds[point] != m_FeatureIds[neigh2])
          ownerLists[m_NodeIds[nodeId 1-4]].insert(m_FeatureIds[point]);
          ownerLists[m_NodeIds[nodeId 1-4]].insert(m_FeatureIds[neigh2]);
        if(k == (zP - 1))
        else if(m_FeatureIds[point] != m_FeatureIds[neigh3])
          ownerLists[m_NodeIds[nodeId 1-4]].insert(m_FeatureIds[point]);
          ownerLists[m_NodeIds[nodeId 1-4]].insert(m_FeatureIds[neigh3]);
