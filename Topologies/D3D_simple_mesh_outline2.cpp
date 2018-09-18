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
      if i = 1
      {
        neigh4 = 
        if(m_FeatureIds[point] != m_FeatureIds[neigh4]) 
        {
          nodeId1
          nodeId2
          nodeId3
          nodeId4
        }
      } 
      if j = 1
      {
        neigh5 = 
        if(m_FeatureIds[point] != m_FeatureIds[neigh5]) 
        {

        }
      } 
      if k = 1
      {
        neigh6 = 
        if(m_FeatureIds[point] != m_FeatureIds[neigh6]) 
        {

        }
      } 

      if(m_FeatureIds[point] != m_FeatureIds[neigh1])
      {
        nodeId1 = (k * (xP + 1) * (yP + 1)) + (j * (xP + 1)) + (i + 1);
        nodeId2 = (k * (xP + 1) * (yP + 1)) + ((j + 1) * (xP + 1)) + (i + 1);
        nodeId3 = ((k + 1) * (xP + 1) * (yP + 1)) + (j * (xP + 1)) + (i + 1);
        nodeId4 = ((k + 1) * (xP + 1) * (yP + 1)) + ((j + 1) * (xP + 1)) + (i + 1);
        
        ownerLists[nodeId1].insert(m_FeatureIds[point]);
        ownerLists[nodeId1].insert(m_FeatureIds[neigh1]);
        ownerLists[nodeId2].insert(m_FeatureIds[point]);
        ownerLists[nodeId2].insert(m_FeatureIds[neigh1]);
        ownerLists[nodeId3].insert(m_FeatureIds[point]);
        ownerLists[nodeId3].insert(m_FeatureIds[neigh1]);
        ownerLists[nodeId4].insert(m_FeatureIds[point]);
        ownerLists[nodeId4].insert(m_FeatureIds[neigh1]);
      }
      if(m_FeatureIds[point] != m_FeatureIds[neigh2])
      {
        nodeId1 = (k * (xP + 1) * (yP + 1)) + ((j + 1) * (xP + 1)) + (i + 1);
        nodeId2 = (k * (xP + 1) * (yP + 1)) + ((j + 1) * (xP + 1)) + i;
        nodeId3 = ((k + 1) * (xP + 1) * (yP + 1)) + ((j + 1) * (xP + 1)) + (i + 1);
        nodeId4 = ((k + 1) * (xP + 1) * (yP + 1)) + ((j + 1) * (xP + 1)) + i;
      
        ownerLists[nodeId1].insert(m_FeatureIds[point]);
        ownerLists[nodeId1].insert(m_FeatureIds[neigh2]);
        ...
      }
      if(m_FeatureIds[point] != m_FeatureIds[neigh3])
      {
        nodeId1 = ((k + 1) * (xP + 1) * (yP + 1)) + (j * (xP + 1)) + (i + 1);
        nodeId2 = ((k + 1) * (xP + 1) * (yP + 1)) + (j * (xP + 1)) + i;
        nodeId3 = ((k + 1) * (xP + 1) * (yP + 1)) + ((j + 1) * (xP + 1)) + (i + 1);
        nodeId4 = ((k + 1) * (xP + 1) * (yP + 1)) + ((j + 1) * (xP + 1)) + i;
        
        ownerLists[nodeId1].insert(m_FeatureIds[point]);
        ownerLists[nodeId1].insert(m_FeatureIds[neigh3]);
        ...
      }
    }
  }
}

















