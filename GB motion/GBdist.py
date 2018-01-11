import numpy as np
import RotRep
O=RotRep.GetSymRotMat(symtype='Cubic')

bsym=np.zeros((24,4,4))
for i in range(24):
    bsym[i][:3,:3]=O[i]
    bsym[i][3,3]=1
    
def spv4dist(b1,b2,moreinfo=False):
"""
Calculate the "distance" between two GB.
Note!:	It is not a good definition when GBs are close to no misorientation. e.g.
	If two GBs both have 0 misorientation but different 'plane normal', this
	definition will give you a nonzero distance, however these two GBs should
	be the same in fact.

Input:	b1,b2
	numpy array, shape=(4,4)
Return:	v4dist
	real,0<v4dist<5, the distance of two boundaries
	bb
	ndarray, shape=(4,4). bb is symmetrically equivalent
	to b2, and it has smallest ||b1-bb|| value

	truei,truej
	integer, indicate which symmetry operators bsym have 
	been used to get bb
	tran
	1 or 2 or 3 or 4, indicates which case has happened, please see
	the source code for more information
"""
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
"""
Get (4,4) matrix from misorientation gx and plane normal cn
"""
    gx=gx.reshape((3,3))
    cn=cn.reshape((3,1))
    b=np.zeros((4,4))
    b[:3,:3]=gx
    b[:3,3]=cn.ravel()
    b[3,:3]=-cn.T.dot(gx).ravel()
    return b

