import numpy as np
from  RotRep import GetSymRotMat
O=GetSymRotMat(symtype='Cubic')

def EAtoG(EA):
    """
    Input: a set of Euler Angle
                size=[3,]
    Output: the corresponding orientation matrix g
                size = [3, 3]
    """
    g = np.zeros((3,3))
    EA = np.radians(EA)
    
    g[0,0]=np.cos(EA[0])*np.cos(EA[2])-np.sin(EA[0])*np.sin(EA[2])*np.cos(EA[1])
    g[0,1]=np.sin(EA[0])*np.cos(EA[2])+np.cos(EA[0])*np.sin(EA[2])*np.cos(EA[1])
    g[0,2]=np.sin(EA[2])*np.sin(EA[1])
    g[1,0]=-np.cos(EA[0])*np.sin(EA[2])-np.sin(EA[0])*np.cos(EA[2])*np.cos(EA[1])
    g[1,1]=-np.sin(EA[0])*np.sin(EA[2])+np.cos(EA[0])*np.cos(EA[2])*np.cos(EA[1])
    g[1,2]=np.cos(EA[2])*np.sin(EA[1])
    g[2,0]=np.sin(EA[0])*np.sin(EA[1])
    g[2,1]=-np.cos(EA[0])*np.sin(EA[1])
    g[2,2]=np.cos(EA[1])
    return g


def dgInFZ(g1, g2, i):
    """
    Input: g1, g2, orientation matrixes
                size = [3, 3]
           i, i+1 = the ID for current triple junction, pass in for error check
    Output: RFvex, the variant of the misorientation, as Rodrigues vector, that lies in the Fundamental Zone
                size = [3,]
            g_FZ, the variant of the misorientation, as orientation matrix, that lies in the Fundamental Zone
                size = [3, 3]
    """
    for j in range(24):
        for k in range(24):
            gg1 = np.dot(O[j,:,:], g1)
            gg2 = np.dot(O[k,:,:], g2)
            
            dg = np.dot(gg1, gg2.T)
            misA = np.arccos(0.5*(np.trace(dg)-1))
# # misA error check-----------------------------------------------------------------
#             if ((0.5*(np.trace(dg)-1) <= 1) & (0.5*(np.trace(dg)-1) >= -1)):
#                 misA = np.arccos(0.5*(np.trace(dg)-1))
#             else:
#                 print 0.5*(np.trace(dg)-1
#                 print 'misA problem at i =', i
#                 return 
# # misA error check-----------------------------------------------------------------
            
            RFvec = [dg[1,2]-dg[2,1], dg[2,0]-dg[0,2], dg[0,1]-dg[1,0]] /(2*np.sin(misA)) * np.tan(misA/2)
            inFZ = ((all(RFvec >= 0)) & (RFvec[0] <= (np.sqrt(2) - 1)) & \
                    (RFvec[0] >= RFvec[1]) & (RFvec[1] >= RFvec[2]) & (sum(RFvec) <= 1))
            if (inFZ == True):
                return (RFvec, dg)
                
            dg = np.dot(gg2, gg1.T)
            misA = np.arccos(0.5*(np.trace(dg)-1))
# # misA error check-----------------------------------------------------------------
#             if ((0.5*(np.trace(dg)-1) <= 1) & (0.5*(np.trace(dg)-1) >= -1)):
#                 misA = np.arccos(0.5*(np.trace(dg)-1))
#             else:
#                 print 0.5*(np.trace(dg)-1
#                 print 'misA problem at i =', i
#                 return
# # misA error check-----------------------------------------------------------------                
            
            RFvec = [dg[1,2]-dg[2,1], dg[2,0]-dg[0,2], dg[0,1]-dg[1,0]] /(2*np.sin(misA)) * np.tan(misA/2)
            inFZ = ((all(RFvec >= 0)) & (RFvec[0] <= (np.sqrt(2) - 1)) & \
                    (RFvec[0] >= RFvec[1]) & (RFvec[1] >= RFvec[2]) & (sum(RFvec) <= 1))
            if (inFZ == True):
                return (RFvec, dg)
            
    print 'error: no copy in FZ. TJ_id =', i+1 
    return ([0, 0, 0], np.zeros((3,3)))


def convertInFZ(EAs, numTJ):
    """
    Input: EAs, the Euler Angles of the 3 grains at a TJ
                size = [numTJ, 3, 3]
                In each group, the data is [TJ directon, EA1, GB1, EA2, GB2, EA3, GB3]
    Output: RFvecs, the equivalent misorientation in Fundamental Zone as rodrigues vector for the 3 GBs at a TJ
                size = [numTJ, 3, 3]
            gs_FZ,  the equivalent misorientation in Fundamental Zone as orientation matrix for the 3 GBs at a TJ
                size = [numTJ, 3, 3, 3]
    Parameters: gs, the misorientation from EA, not 
    """
    RFvecs = np.zeros((numTJ, 3, 3))
    dgs_FZ = np.zeros((numTJ, 3, 3, 3))
    for i in range(numTJ):
        g1 = EAtoG(EAs[i, 0, :])
        g2 = EAtoG(EAs[i, 1, :])
        g3 = EAtoG(EAs[i, 2, :])
        (RFvec1, dg1_FZ) = dgInFZ(g1, g2, i)
        (RFvec2, dg2_FZ) = dgInFZ(g2, g3, i)
        (RFvec3, dg3_FZ) = dgInFZ(g3, g1, i)
        RFvecs[i, 0, :] = RFvec1
        RFvecs[i, 1, :] = RFvec2
        RFvecs[i, 2, :] = RFvec3
        dgs_FZ[i, 0, :, :] = dg1_FZ
        dgs_FZ[i, 1, :, :] = dg2_FZ
        dgs_FZ[i, 2, :, :] = dg3_FZ
    return (RFvecs, dgs_FZ)
    
    

# (RFvecs, dgs_FZ) = convertInFZ(EAs, numTJ)


def diffMisorientations(dg1, dg2):
    """
    Input: dg1 & dg2, the misorientation of the first grain boundary 
                size = [3, 3]
    Output: ddg, the misorientation between the two misorientations
                size = [3, 3]
    Parameters: disA, the minimum misorientation anlge
                
    !!NOTICE!! dg1 and dg2 needs to be in the Fundamental Zone.
    """
    disA = 1000
    for i in range(24):
        ddg = np.dot(np.dot(O[i,:,:], dg1), dg2.T)
        misA = np.arccos(0.5*(np.trace(ddg)-1))
        misA = np.degrees(misA)
        if misA < disA:
            disA = misA
        ddg = np.dot(np.dot(O[i,:,:], dg2), dg1.T)
        misA = np.arccos(0.5*(np.trace(ddg)-1))
        misA = np.degrees(misA)
        if misA < disA:
            disA = misA
    return disA
    
def AAtoG(phi, n):
    """
    Input: (phi, n), the orientation as axis-angle pair 
                phis is scalar, in degrees
                n.size = [3,]
    Output: dg, the orientation as matrix
                
    !!NOTICE!! phi is in degrees
    """
    phi = np.radians(phi)
    n = n/np.sqrt(n[0]**2 + n[1]**2 + n[2]**2)
    dg = np.zeros((3,3))
    dg[0,0] = np.cos(phi)+(1.0-np.cos(phi))*(n[0]**2)
    dg[0,1] = n[0]*n[1]*(1.0-np.cos(phi))-n[2]*np.sin(phi)
    dg[0,2] = n[0]*n[2]*(1.0-np.cos(phi))+n[1]*np.sin(phi)
    dg[1,0] = n[0]*n[1]*(1.0-np.cos(phi))+n[2]*np.sin(phi)
    dg[1,1] = np.cos(phi)+(1.0-np.cos(phi))*(n[1]**2)
    dg[1,2] = n[2]*n[1]*(1.0-np.cos(phi))-n[0]*np.sin(phi)
    dg[2,0] = n[0]*n[2]*(1.0-np.cos(phi))-n[1]*np.sin(phi)
    dg[2,1] = n[1]*n[2]*(1.0-np.cos(phi))+n[0]*np.sin(phi)
    dg[2,2] = np.cos(phi)+(1.0-np.cos(phi))*(n[2]**2)
    return dg