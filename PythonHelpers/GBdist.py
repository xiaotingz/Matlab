import numpy as np

O = GetSymRotMat(symtype='Cubic')

bsym=np.zeros((24,4,4))
for i in range(24):
    bsym[i][:3,:3]=O[i]
    bsym[i][3,3]=1
    
def spv4dist(b1,b2,moreinfo=False):
    v4dist=5.0
    bb=np.empty((4,4))
    truei=-1
    truej=-1
    tran=0
    for ii in range(24):
        for jj in range(24):
            ttt=bsym[ii].dot(b2).dot(bsym[jj])
            tt=b1-ttt
            compar=np.trace(tt.dot(tt.T))
            if compar < v4dist:
                v4dist = compar
                bb=ttt.copy()
                truei=ii
                truej=jj
                tran=1
            tt=b1-ttt.T
            compar=np.trace(tt.dot(tt.T))
            if compar < v4dist:
                v4dist = compar
                bb=ttt.T.copy()
                truei=ii
                truej=jj
                tran=2
            
            ttt[3,0:3]*=-1
            ttt[0:3,3]*=-1
            tt=b1-ttt
            compar=np.trace(tt.dot(tt.T))
            if compar < v4dist:
                v4dist = compar
                bb=ttt.copy()
                truei=ii
                truej=jj
                tran=3
            tt=b1-ttt.T
            compar=np.trace(tt.dot(tt.T))
            if compar < v4dist:
                v4dist = compar
                bb=ttt.T.copy()
                truei=ii
                truej=jj
                tran=4
    if moreinfo:
        return v4dist,bb,truei,truej,tran
    return v4dist,bb

def convertB(gx,cn):
    gx=gx.reshape((3,3))
    cn=cn.reshape((3,1))
    b=np.zeros((4,4))
    b[:3,:3]=gx
    b[:3,3]=cn.ravel()
    b[3,:3]=-cn.T.dot(gx).ravel()
    return b

def GetSymRotMat(symtype='Cubic'):
    """
    return an array of active rotation matrices of the input crystal symmetry

    Parameters
    ----------
    symtype: string
            Symmetry type of crystal. For now only 'Cubic' and 'Hexagonal' are implemented. 

    Returns
    ----------
    m:  ndarray
        A three dimensional numpy array, which has the shape (n,3,3). 
    """
    if symtype == 'Cubic':
        m=np.zeros((24,3,3))
        m[0][0,1]=1
        m[0][1,0]=-1
        m[0][2,2]=1

        m[1][0,0]=-1
        m[1][1,1]=-1
        m[1][2,2]=1

        m[2][0,1]=-1
        m[2][1,0]=1
        m[2][2,2]=1

        m[3][0,2]=-1
        m[3][1,1]=1
        m[3][2,0]=1

        m[4][0,0]=-1
        m[4][1,1]=1
        m[4][2,2]=-1

        m[5][0,2]=1
        m[5][1,1]=1
        m[5][2,0]=-1

        m[6][0,0]=1
        m[6][1,2]=1
        m[6][2,1]=-1

        m[7][0,0]=1
        m[7][1,1]=-1
        m[7][2,2]=-1

        m[8][0,0]=1
        m[8][1,2]=-1
        m[8][2,1]=1

        m[9][0,1]=1
        m[9][1,2]=1
        m[9][2,0]=1

        m[10][0,2]=1
        m[10][1,0]=1
        m[10][2,1]=1

        m[11][0,2]=-1
        m[11][1,0]=1
        m[11][2,1]=-1

        m[12][0,1]=1
        m[12][1,2]=-1
        m[12][2,0]=-1

        m[13][0,2]=1
        m[13][1,0]=-1
        m[13][2,1]=-1

        m[14][0,1]=-1
        m[14][1,2]=-1
        m[14][2,0]=1

        m[15][0,2]=-1
        m[15][1,0]=-1
        m[15][2,1]=1

        m[16][0,1]=-1
        m[16][1,2]=1
        m[16][2,0]=-1

        m[17][0,0]=-1
        m[17][1,2]=1
        m[17][2,1]=1

        m[18][0,2]=1
        m[18][1,1]=-1
        m[18][2,0]=1

        m[19][0,1]=1
        m[19][1,0]=1
        m[19][2,2]=-1

        m[20][0,0]=-1
        m[20][1,2]=-1
        m[20][2,1]=-1

        m[21][0,2]=-1
        m[21][1,1]=-1
        m[21][2,0]=-1

        m[22][0,1]=-1
        m[22][1,0]=-1
        m[22][2,2]=-1

        m[23][0,0]=1
        m[23][1,1]=1
        m[23][2,2]=1
        
        return m
    elif symtype == 'Hexagonal':
        m=np.zeros((12,3,3))
        m[0][0,0]=0.5
        m[0][1,1]=0.5
        m[0][2,2]=1
        m[0][0,1]=-np.sqrt(3)*0.5
        m[0][1,0]=np.sqrt(3)*0.5

        m[1]=m[0].dot(m[0])
        m[2]=m[1].dot(m[0])
        m[3]=m[2].dot(m[0])
        m[4]=m[3].dot(m[0])
        m[5]=np.eye(3)

        m[6][0,0]=1
        m[6][1,1]=-1
        m[6][2,2]=-1

        m[7]=m[0].dot(m[6])
        m[8]=m[1].dot(m[6])
        m[9]=m[2].dot(m[6])
        m[10]=m[3].dot(m[6])
        m[11]=m[4].dot(m[6])

        return m
    else:
        print "not implemented yet"
        return 0
