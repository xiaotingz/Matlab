L_EA = np.zeros((3*numTJ, 3))
R_EA = np.zeros((3*numTJ, 3))
n  = np.zeros((3*numTJ, 3))
for i in range(numTJ):
        L_EA[3*i,:] = EAs[i,1,:]
        R_EA[3*i,:] = EAs[i,2,:]
        n[3*i,:] = norms[i,0,:]
        L_EA[3*i + 1,:] = EAs[i,2,:]
        R_EA[3*i + 1,:]  = EAs[i,0,:]
        n[3*i + 1,:] =norms[i,1,:]
        L_EA[3*i + 2,:] = EAs[i,0,:]
        R_EA[3*i + 2,:] = EAs[i,1,:]
        n[3*i + 2,:] = norms[i,2,:]   

PHI=np.degrees(np.arccos(n[:,2])).reshape((n.shape[0],1))
THETA=np.degrees(np.arctan2(n[:,1],n[:,0])).reshape((n.shape[0],1))
THETA[THETA<0] += 360

doomy = np.ones((len(PHI),1))*1200
MorE = MorE.reshape((len(PHI),1))

toPrint = np.hstack((L_EA, R_EA, PHI, THETA, doomy, MorE))
np.savetxt('60000_181116_THETA360.gbdat', toPrint, delimiter=',', fmt='%3.4f '*8+'%4d '+'%3.7f')
