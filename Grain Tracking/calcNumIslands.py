def exploreIsland(mat, new, i, j, k):
    next_pos = [[i, j, k]]
    grain_id = mat[i, j, k]
    
    while next_pos:
        [i, j, k] = next_pos.pop()
        legal_pos = (0 <= i < mat.shape[0] and 0 <= j < mat.shape[1] and 0 <= k < mat.shape[2])
        explore = legal_pos and new[i, j, k] and mat[i, j, k] == grain_id
        if explore:
            new[i, j, k] = False
            next_pos.append([i+1, j, k])
            next_pos.append([i, j+1, k])
            next_pos.append([i, j, k+1])
            next_pos.append([i-1, j, k])
            next_pos.append([i, j-1, k])
            next_pos.append([i, j, k-1])
    return new

def calcNumOfIslands(mat):
    new = np.ones(mat.shape).astype(bool)
    cnt = 0 
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            for k in range(mat.shape[2]):
                if new[i, j, k]:
                    cnt += 1
                    new = exploreIsland(mat, new, i, j, k)
    return cnt

print 'an4, #(unique grainIds) = ', np.unique(an4).shape[0], ';   #(islands) = ', calcNumOfIslands(an4)    
print 'an5, #(unique grainIds) = ', np.unique(an5).shape[0], ';   #(islands) = ', calcNumOfIslands(an5)       
# num_grains_an4 = calcNumOfIslands(an4)