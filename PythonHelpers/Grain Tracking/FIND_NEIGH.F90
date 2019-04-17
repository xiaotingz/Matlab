  MODULE globals

    REAL*4, SAVE    :: pi, eps

  END MODULE globals



  MODULE symmetry_operations

    INTEGER*4, SAVE   :: nse, nselau

    REAL*4, SAVE    :: gcsym(24,3,3), bsym(24,4,4)

  END MODULE symmetry_operations


  MODULE param_limits

!***********************************************************************

! ffmax - upper limits on boundary parameters (angles) - constant.

! afmax - upper limits on boundary parameters (cos of angles phi and thata) - constant.

! ndiv  - number of cells in each of 5 directions of boundary parameters space - constant.

! af    - 5 dimensional specification of a cell (cos of angles phi and thata).

!***********************************************************************

    INTEGER*4, SAVE, DIMENSION(5) :: ndiv

    REAL*4, SAVE, DIMENSION(5)  :: afmax, ffmax
    
  END MODULE 


  MODULE cell_and_huge

    INTEGER*4, PARAMETER :: np=9, ksym=36, n_max=4*np**5


  END MODULE 

  MODULE triple_junction_module

!***********************************************************************

!  'sn_i' normals to the boundaries (sample coordinate system)

!  'sl' direction of the triple line (sample coordinate system)

!  'g_i' grain orientations.

!  'gx_i' grain misorientations.

!  'cn_i' normals to the boundaries (crystal coordinate system)

!***********************************************************************

    REAL*4, SAVE, DIMENSION(3,3) :: g1,g2,g3, gx1,gx2,gx3

    REAL*4, SAVE, DIMENSION(3) :: cn1,cn2,cn3, sn1,sn2,sn3, sl

  END MODULE triple_junction_module


  MODULE all_bs
    INTEGER*4, PARAMETER :: nJ=120000
    REAL*4                bs(nJ,3,4,4)
    INTEGER*4, SAVE, DIMENSION(nJ,3,36) :: idxcell
  END MODULE all_bs




!**********************************************************************

!**********************************************************************

  PROGRAM transformator_w_torque

  USE cell_and_huge

  USE globals

  USE param_limits


  WRITE(*,'(30(/),32X,''Yu-Feng Shen, A.Morawiec'',5(/))')

! Initialization of constants.

  pi=ACOS(-1.)

  eps=0.000001

  ndiv=(/np,np,np,np,4*np/)

  f1max=pi/2.;    fmax=pi/2.;  f2max=pi/2.;    xtmax=pi/2.; xpmax=2.*pi;  

  afmax(1)=f1max; afmax(2)=1.; afmax(3)=f2max; afmax(4)=1.; afmax(5)=xpmax; 

  ffmax(1)=f1max; ffmax(2)=fmax;   

  ffmax(3)=f2max; ffmax(4)=xtmax; ffmax(5)=xpmax

  CALL symbs_gb

  OPEN(1,FILE='triples.dat_120000',STATUS='OLD')



  CALL read_bs(ksym)
!  CALL cal_idxcell(ksym)  
  CALL read_idxcell(ksym)
  CALL find_nearest

  CLOSE(1)

  STOP

  END





!***********************************************************************

  SUBROUTINE read_bs(ksym)

!***********************************************************************

! Assignes cell numbers 'num' (= beta) to boundaries.

! Assignes coefficients 'coeffo' (=A) to boundaries.

! Writes both to a file.

!***********************************************************************

  USE triple_junction_module
  USE globals
  USE all_bs

! The loop over triple junctions

  readloop : DO imain=1,nJ
    CALL read_dat
    CALL bconvert(gx1,cn1,bs(imain,1,:,:))
    CALL bconvert(gx2,cn2,bs(imain,2,:,:))
    CALL bconvert(gx3,cn3,bs(imain,3,:,:))
!    CALL read_dat
  ENDDO readloop 

  RETURN


  END

!***********************************************************************

  SUBROUTINE cal_idxcell(ksym)

!***********************************************************************

! Assignes cell numbers 'num' (= beta) to boundaries.

! Assignes coefficients 'coeffo' (=A) to boundaries.

! Writes both to a file.

!***********************************************************************

  USE triple_junction_module
  USE globals
  USE all_bs

  OPEN(2,FILE='idxcell_120000.txt',STATUS='NEW')
! The loop over triple junctions

  !$OMP PARALLEL DO
  DO imain=1,nJ
    CALL find_symeq_num(bs(imain,1,:,:),ksym,imain,1)
    CALL find_symeq_num(bs(imain,2,:,:),ksym,imain,2)
    CALL find_symeq_num(bs(imain,3,:,:),ksym,imain,3)
  ENDDO 
  !$OMP END PARALLEL DO

  DO imain=1,nJ
    WRITE(2,*) idxcell(imain,1,:)
    WRITE(2,*) idxcell(imain,2,:)
    WRITE(2,*) idxcell(imain,3,:)
  ENDDO

  CLOSE(2)
  RETURN


  END

!***********************************************************************

  SUBROUTINE read_idxcell(ksym)

!***********************************************************************

! Assignes cell numbers 'num' (= beta) to boundaries.

! Assignes coefficients 'coeffo' (=A) to boundaries.

! Writes both to a file.

!***********************************************************************

  USE triple_junction_module
  USE globals
  USE all_bs

  OPEN(2,FILE='idxcell_120000.txt',STATUS='OLD')
! The loop over triple junctions

  DO imain=1,nJ
    READ(2,*) idxcell(imain,1,:)
    READ(2,*) idxcell(imain,2,:)
    READ(2,*) idxcell(imain,3,:)
  ENDDO 

  CLOSE(2)
  RETURN


  END
!***********************************************************************

  SUBROUTINE find_nearest

!***********************************************************************

! Assignes cell numbers 'num' (= beta) to boundaries.

! Assignes coefficients 'coeffo' (=A) to boundaries.

! Writes both to a file.

!***********************************************************************
  USE triple_junction_module
  USE all_bs
  USE symmetry_operations
  DIMENSION idx1(36),idx2(36),xM(3,3),xG2(3,3)
  INTEGER*4,ALLOCATABLE::          idxbs(:,:,:)
  REAL*4,ALLOCATABLE::            disbs(:,:,:)
  REAL*4,ALLOCATABLE::            Tranbs(:,:,:,:,:)!last two 3 are (3,3) Matrix
  ALLOCATE(idxbs(nJ,3,50),disbs(nJ,3,50),Tranbs(nJ,3,3,3,50))
  OPEN(3,FILE='idxbs_120000.txt',STATUS='NEW')
  OPEN(4,FILE='disbs_120000.txt',STATUS='NEW')
  OPEN(5,FILE='Tranbs_120000.txt',STATUS='NEW')
  !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(Tranbs,idxbs,disbs,idxcell,bs,gcsym)
  DO imain=1,nJ
!  WRITE(*,*) imain
  DO ib=1,3
    idxbs(imain,ib,:)=0
    disbs(imain,ib,:)=5.
    idx1(:)=idxcell(imain,ib,:)
    ii=1
    DO jmain=1,nJ
    DO jb=1,3
        IF(ii.LE.50.AND.(jmain.NE.imain.OR.jb.NE.ib)) THEN
        idx2(:)=idxcell(jmain,jb,:)
        IF(intersect(idx1,idx2).GT.1) THEN
            CALL spv4dist(bs(imain,ib,:,:),bs(jmain,jb,:,:),x,itrue,jtrue,icase)
            IF(x.LT. 0.023) THEN
                xG2=bs(jmain,jb,1:3,1:3)
                IF(icase.EQ.1) THEN
                    xM(:,:)=gcsym(itrue,:,:)
                ELSE IF(icase.EQ.2) THEN
                    xM(:,:)=-TRANSPOSE(MATMUL(xG2,gcsym(jtrue,:,:)))
                ELSE IF(icase.EQ.3) THEN
                    xM(:,:)=-gcsym(itrue,:,:)
                ELSE IF(icase.EQ.4) THEN
                    xM(:,:)=TRANSPOSE(MATMUL(xG2,gcsym(jtrue,:,:)))
                END IF
                Tranbs(imain,ib,:,:,ii)=xM(:,:)
                disbs(imain,ib,ii)=x
                idxbs(imain,ib,ii)=jb+3*(jmain-1)
                ii=ii+1
            ENDIF
        ENDIF
        ENDIF
    ENDDO
    ENDDO
    IF(imain.EQ.1.OR.imain.EQ.3)THEN
        WRITE(*,*)idxbs(imain,ib,:)
        WRITE(*,*)disbs(imain,ib,:)
    ENDIF
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO

  DO imain=1,nJ
  DO ib=1,3
    WRITE(3,*) idxbs(imain,ib,:)
    WRITE(4,*) disbs(imain,ib,:)
    WRITE(5,*) Tranbs(imain,ib,:,:,:)
  ENDDO
  ENDDO
  
  CLOSE(3)
  CLOSE(4)
  CLOSE(5)
  DEALLOCATE(idxbs,disbs,Tranbs)
  RETURN


  END


!***********************************************************************

  SUBROUTINE find_symeq_num(b, ksym,i,j)

!***********************************************************************

! Calculates cell numbers 'num' and appropriate  coefficients 'coeffo'

! for a boundary 'b' between grains with orientations 'g1' and 'g2'.

! i is the TJ number, j = 1,2,3 the three boundaries of that TJ

!***********************************************************************

  USE param_limits

  USE symmetry_operations

  USE all_bs

  INTEGER*4 :: i,j
  REAL, DIMENSION(3,3) :: b(4,4),bt(4,4),bx(4,4)
  REAL*4,DIMENSION(5) :: af

  LOGICAL testor


  bt=TRANSPOSE(b)

  nsymeq=0

! The loop over symmetry operations.

  DO ii=1,nselau

  DO jj=1,nselau

!   The boundary between 'g1' and 'g2'.

    CALL equivalox(ii,jj,b,bx)

    CALL matrix_to_param(bx,af,isignum)

!   Test whether 'bx' is in the enlarged asymmetric zone.

    testor=(af(1).le.ffmax(1)).AND.(af(2).le.ffmax(2)).AND.(af(3).le.ffmax(3)).AND.(af(4).le.ffmax(4)).AND.(af(5).le.ffmax(5))

    IF(testor) THEN

      af(2)=COS(af(2)) ; af(4)=COS(af(4)) ;

      nsymeq=nsymeq+1

      IF(nsymeq.gt.ksym) RETURN

!     Cell number

      idxcell(i,j,nsymeq)=neuler_to_cell(af,afmax,ndiv)



    ENDIF

!   The boundary between 'g2' and 'g1'.

    CALL equivalox(ii,jj,bt,bx)

    CALL matrix_to_param(bx,af,isignum)

!   Test whether 'bx' is in the enlarged asymmetric zone.

    testor=(af(1).le.ffmax(1)).AND.(af(2).le.ffmax(2)).AND.(af(3).le.ffmax(3)).AND.(af(4).le.ffmax(4)).AND.(af(5).le.ffmax(5))

    IF(testor) THEN

      af(2)=COS(af(2)) ; af(4)=COS(af(4)) ;

      nsymeq=nsymeq+1

      IF(nsymeq.gt.ksym) RETURN

!     Cell number

      idxcell(i,j,nsymeq)=neuler_to_cell(af,afmax,ndiv)


    ENDIF

  ENDDO

  ENDDO

  END


!***********************************************************************

  FUNCTION neuler_to_cell(af,afmax,ndiv)

!***********************************************************************

! Replaces the 5 dimensional specification of a cell (af)

! by one dimensional (neuler_to_cell).

!***********************************************************************

  DIMENSION af(5),afmax(5),ndiv(5),m(5)

  m=MIN(INT((af*ndiv)/afmax),ndiv-1)

  neuler_to_cell=1+m(1)+m(2)*ndiv(1)+m(3)*ndiv(1)*ndiv(2)+ &

   m(4)*ndiv(1)*ndiv(2)*ndiv(3)+m(5)*ndiv(1)*ndiv(2)*ndiv(3)*ndiv(4)

  RETURN

  END



!***********************************************************************

  SUBROUTINE matrix_to_param(b,fp,isignum)

!***********************************************************************

! Macroscopic boundary parameters fp(5) from boundary matrix b.

!***********************************************************************

  USE globals

  DIMENSION b(4,4), fp(5), g(3,3), xn(3)

  g=b(:3,:3)

  CALL MEul(g,fp(1),fp(2),fp(3))

  xn=b(1:3,4)

  CALL Vec_to_polar(xn,fp(4),fp(5))

  isignum=1

! change the direction to positive 'z'.

  IF(fp(4).gt.pi/2) THEN

     fp(4)=pi-fp(4)

     fp(5)=fp(5)+pi

     IF(fp(5).gt.(2.*pi)) fp(5)=fp(5)-2.*pi

     isignum=-1

  ENDIF

  RETURN

  END



      subroutine EMat(f1,f,f2,g)

!***********************************************************************

!     Euler angles --> rotation matrix                                 *

!***********************************************************************

      dimension g(3,3)

      sf=sin(f)

      cf=cos(f)

      sf1=sin(f1)

      cf1=cos(f1)

      sf2=sin(f2)

      cf2=cos(f2)

      g(1,1)=cf1*cf2-sf1*sf2*cf

      g(1,2)=sf1*cf2+cf1*sf2*cf

      g(1,3)=sf2*sf

!***

      g(2,1)=-cf1*sf2-sf1*cf2*cf

      g(2,2)=-sf1*sf2+cf1*cf2*cf

      g(2,3)=cf2*sf

!***

      g(3,1)=sf1*sf

      g(3,2)=-cf1*sf

      g(3,3)=cf

      end



      subroutine MEul(g,a1,a,a2)

!***********************************************************************

!     Rotation matrix --> Euler angles                                 *

!***********************************************************************

      USE globals

      dimension g(3,3)

      x=g(3,3)

      sf=sqrt(1.-min(x*x,1.))

      if(sf.gt.eps) then

        a1=-g(3,2)/sf

        a1=acoss(a1)

        if(g(3,1).lt.0.) a1=2.*pi-a1

        a=acos(x)

        a2=g(2,3)/sf      

        a2=acoss(a2)

        if(g(1,3).lt.0.) a2=2.*pi-a2

      else 

        a1=g(1,1)

        a1=acoss(a1)

        if(g(1,2).lt.0.) a1=2.*pi-a1

        a=acoss(x)

        a2=0.

      endif

      end



      function acoss(ang)

      angl=ang

      if(angl.gt.1) angl=1

      if(angl.lt.-1) angl=-1

      acoss=acos(angl)

      end



!***********************************************************************

      subroutine Vec_to_polar(x,theta,psi)

!***********************************************************************

!     Polar coordinates (theta,psi) of a vector x.

!***********************************************************************

      USE globals

      dimension x(3)

      theta=acoss(x(3))

      sf=sqrt(1.-min(x(3)*x(3),1.))

      if(sf.gt.eps) then

        psi=acoss(x(1)/sf)

        if(x(2).lt.0.) psi=2.*pi-psi

      else 

        psi=acoss(x(1))

        if(x(2).lt.0.) psi=2.*pi-psi

      endif

      end



!***********************************************************************

      SUBROUTINE symbs_gb

!***********************************************************************

!     Reads symmetry operations

!***********************************************************************

      USE symmetry_operations



      OPEN(7,file='symmetry.prp',status='old',err=11)

      READ(7,*,err=11,END=11)

      READ(7,*,err=11,END=11)

      READ(7,*,err=11,END=11) is,nse,nselau

      bsym=0.

      DO 1 nn=1,nselau

         READ(7,*,err=11,END=11)

         DO 2 i=1,3

   2        READ(7,*,err=11,END=11) (gcsym(nn,i,j),j=1,3)

         bsym(nn,1:3,1:3)=gcsym(nn,:,:)

         bsym(nn,4,4)=1.

  1   CONTINUE

      CLOSE(7,err=11)

      RETURN

  11  WRITE(*,*) 'Cannot open (read) SYMMETRY.PRP.'

      STOP

      END



!***********************************************************************

      SUBROUTINE equivalox(ii,jj,bin,bout)

!***********************************************************************

! Boundary matrix 'bout' equivalent to boundary matrix 'bin' by 

! symmetry operations with numbers 'ii' and 'jj'.

!***********************************************************************

      USE symmetry_operations

      dimension bin(4,4),bout(4,4)

      do 2 i=1,4

      do 2 j=1,4

         bout(i,j)=0.

         do 2 k=1,4

         do 2 l=1,4

  2            bout(i,j)=bout(i,j)+bsym(ii,i,k)*bin(k,l)*bsym(jj,l,j)

      end





!***********************************************************************

  SUBROUTINE read_dat

!***********************************************************************

! Reads 'triples.dat' file. 

! Calculates normals in crystal coordinate system.

! Calculates misorientations.

!***********************************************************************

  USE triple_junction_module

  USE globals



  deg_rad=pi/180.



  READ(1,*)

  READ(1,*) sl(:)

  CALL fnormv(sl)



  READ(1,*) f1,f,f2

  READ(1,*) sn1(:)

  CALL fnormv(sn1)

  f1=deg_rad*f1 ; f=deg_rad*f ; f2=deg_rad*f2 ; 

  CALL EMat(f1,f,f2,g1)



  READ(1,*) f1,f,f2

  READ(1,*) sn2(:)

  CALL fnormv(sn2)

  f1=deg_rad*f1 ; f=deg_rad*f ; f2=deg_rad*f2 ; 

  CALL EMat(f1,f,f2,g2)



  READ(1,*) f1,f,f2

  READ(1,*) sn3(:)

  CALL fnormv(sn3)

  f1=deg_rad*f1 ; f=deg_rad*f ; f2=deg_rad*f2 ; 

  CALL EMat(f1,f,f2,g3)



! Calculates normals in crystal coordinate system

  cn1=MATMUL(g2,sn1)

  cn2=MATMUL(g3,sn2)

  cn3=MATMUL(g1,sn3)



! Calculates misorientations

  gx1=MATMUL(g2,TRANSPOSE(g3))

  gx2=MATMUL(g3,TRANSPOSE(g1))

  gx3=MATMUL(g1,TRANSPOSE(g2))



  END



!***********************************************************************

  SUBROUTINE bconvert(g,xn,b)

!***********************************************************************

! Boundary matrix b from misorientation matrix g and inclination xn.

!***********************************************************************

  DIMENSION g(3,3),xn(3),b(4,4), xn2(3)

  b(1:3,1:3)=g

  b(1:3,4)=xn

  xn2=-MATMUL(TRANSPOSE(g),xn)

  b(4,1:3)=xn2

  b(4,4)=0.

  RETURN

  END



SUBROUTINE fnormv(x)

REAL*4  :: x(3)



xnorm=SQRT(DOT_PRODUCT(x,x))

IF(xnorm < 0.00001) THEN

 write(*,*)' ZERO, ZERO, ZERO'

 RETURN

ENDIF

x=x/xnorm

END SUBROUTINE fnormv

      FUNCTION xmdist(g1,g2)
!***********************************************************************
!     Distance between two misorientations 'g1' and 'g2'.
!     The table gcsym(,,) contains matrices of 'nse' crystal symmetry
!     elements. 
!***********************************************************************

      USE symmetry_operations
      DIMENSION g1(3,3),g2(3,3)
      xmdist=-1.
      DO 1 ii=1,nselau
      DO 1 jj=1,nselau
         omeg=0.
         DO 2 i=1,3
         DO 2 j=1,3
         DO 2 k=1,3
         DO 2 l=1,3
  2         omeg=omeg+gcsym(ii,i,j)*g1(i,k)*gcsym(jj,k,l)*g2(j,l)
         IF(omeg.gt.xmdist) xmdist=omeg 
         omeg=0.
         DO 3 i=1,3
         DO 3 j=1,3
         DO 3 k=1,3
         DO 3 l=1,3
  3         omeg=omeg+gcsym(ii,i,j)*g1(i,k)*gcsym(jj,k,l)*g2(l,j)
         IF(omeg.gt.xmdist) xmdist=omeg 
  1   CONTINUE
!      xmdist=acoss((xmdist-1.)/2.)
      xmdist=3.-xmdist
      END


      FUNCTION x4dist(b1,b2)
!***********************************************************************
!     Distance between two 'boundaries' 'b1' and 'b2'.
!     The table bsym(,,) contains matrices of 'nse' crystal symmetry
!     elements. 
!***********************************************************************
      USE symmetry_operations
      DIMENSION b1(4,4),b2(4,4),t(4,4),tt(4,4),ttt(4,4)
      x4dist=5.

      DO 1 ii=1,nselau
      DO 1 jj=1,nselau
         ttt=MATMUL(MATMUL(bsym(ii,:,:),b2),bsym(jj,:,:))

         tt=(b1-ttt)
         t=MATMUL(tt,TRANSPOSE(tt))
         compar=t(1,1)+t(2,2)+t(3,3)+t(4,4)
         IF(compar.LT.x4dist) x4dist=compar

         tt=(TRANSPOSE(b1)-ttt)
         t=MATMUL(tt,TRANSPOSE(tt))
         compar=t(1,1)+t(2,2)+t(3,3)+t(4,4)
         IF(compar.LT.x4dist) x4dist=compar

         ttt(4,1:3)=-ttt(4,1:3)
         ttt(1:3,4)=-ttt(1:3,4)

         tt=(b1-ttt)
         t=MATMUL(tt,TRANSPOSE(tt))
         compar=t(1,1)+t(2,2)+t(3,3)+t(4,4)
         IF(compar.LT.x4dist) x4dist=compar

         tt=(TRANSPOSE(b1)-ttt)
         t=MATMUL(tt,TRANSPOSE(tt))
         compar=t(1,1)+t(2,2)+t(3,3)+t(4,4)
         IF(compar.LT.x4dist) x4dist=compar
  1   CONTINUE
      END

      FUNCTION intersect(iar1,iar2)
          DIMENSION iar1(36),iar2(36)
          intersect=0
          DO ii=1,36
          DO jj=1,36
            IF(iar1(ii).EQ.iar2(jj)) intersect=intersect+1
          ENDDO
          ENDDO
      END


      SUBROUTINE spv4dist(b1,b2,v4dist,itrue,jtrue,icase)
!***********************************************************************
!     Distance between two 'boundaries' 'b1' and 'b2'.
!     The table bsym(,,) contains matrices of 'nse' crystal symmetry
!     elements. 
!     It returns 'v4dist=x4dist' and 'bb' - symmetrically equivalent 
!     to 'b2' and nearest to 'b1'.
!***********************************************************************
      USE symmetry_operations
      REAL*4, DIMENSION(4,4) :: b1,b2,t,tt,ttt
 
      v4dist=5.
      DO 1 ii=1,nselau
      DO 1 jj=1,nselau
         ttt=MATMUL(MATMUL(bsym(ii,:,:),b2),bsym(jj,:,:))

         tt=(b1-ttt)
         t=MATMUL(tt,TRANSPOSE(tt))
         compar=t(1,1)+t(2,2)+t(3,3)+t(4,4)
         IF(compar.LT.v4dist) THEN 
           v4dist=compar
           itrue=ii
           jtrue=jj
           icase=1
         ENDIF

         tt=(b1-TRANSPOSE(ttt))
         t=MATMUL(tt,TRANSPOSE(tt))
         compar=t(1,1)+t(2,2)+t(3,3)+t(4,4)
         IF(compar.LT.v4dist) THEN 
           v4dist=compar
           itrue=ii
           jtrue=jj
           icase=2
         ENDIF

         ttt(4,1:3)=-ttt(4,1:3)
         ttt(1:3,4)=-ttt(1:3,4)

         tt=(b1-ttt)
         t=MATMUL(tt,TRANSPOSE(tt))
         compar=t(1,1)+t(2,2)+t(3,3)+t(4,4)
         IF(compar.LT.v4dist) THEN 
           v4dist=compar
           itrue=ii
           jtrue=jj
           icase=3
         ENDIF

         tt=(b1-TRANSPOSE(ttt))
         t=MATMUL(tt,TRANSPOSE(tt))
         compar=t(1,1)+t(2,2)+t(3,3)+t(4,4)
         IF(compar.LT.v4dist) THEN 
           v4dist=compar
           itrue=ii
           jtrue=jj
           icase=4
         ENDIF

  1   CONTINUE
      END
