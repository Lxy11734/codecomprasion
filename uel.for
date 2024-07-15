! ==========================// UEL //=================================== 
! User Subroutine UEL for Abaqus:
! Element type U1: quadrilateral plane stress (strain) element
! DOF 1-2 : Displacement
!
! Material properties to be given through the input file (*.inp), are:
! PROPS(1) = E     --- Young's modulus
! PROPS(2) = nu    --- Possion's ratio
! PROPS(3) = thck  --- thickness
!
! -----------------------------------------------------------------------
!         Date                Programmer
!      2024/05/30               L. Min
!________________________________________________________________________    
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,ndofel,NRHS,nsvars,
     *     PROPS,npropse,COORDS,mcrd,nnode,U,DU,V,A,JTYPE,TIME,DTIME,
     *     KSTEP,kinc,jelem,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     *     NPREDF,LFLAGS,mlvarx,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     *     PERIOD)
      INCLUDE 'ABA_PARAM.INC'
      INTEGER nsvars,nrhs,ndofel,mdload,mlvarx,npredf,mcrd,npropse,
     * nnode,jtype,kinc,ndload,njprop
      REAL(8) SVARS(nsvars),ENERGY(8),PROPS(npropse),
     * COORDS(mcrd,nnode),U(ndofel),DU(mlvarx,1),V(ndofel),A(ndofel),
     * TIME(2),PARAMS(3),JDLTYP(mdload,*),ADLMAG(16,1),
     * DDLMAG(mdload,*),PREDEF(2,npredf,nnode),LFLAGS(*),JPROPS(*),
     * AMATRX(ndofel,ndofel),RHS(mlvarx,1),DTIME,PNEWDT,PERIOD
      
      ! Calculate RHS & AMATRX
      CALL ISOTROPIC(SVARS,PROPS,COORDS,U,AMATRX,RHS,nsvars,
     *                    npropse,mcrd,nnode,ndofel,mlvarx,kinc,jelem)
      
      
      END SUBROUTINE UEL
! ////////////////////////////////////////////////////////////////////////    
      
! Calculating AMATRX and RHS
      SUBROUTINE ISOTROPIC(SVARS,PROPS,COORDS,U,AMATRX,RHS,nsvars,
     *      npropse,mcrd,nnode,ndofel,mlvarx,kinc,jelem)
      IMPLICIT NONE
      INTEGER nsvars,npropse,mcrd,nnode,ndofel,mlvarx,kinc,jelem,
     * ipt,i,j,k,l
      REAL(8) SVARS(nsvars),PROPS(npropse),COORDS(mcrd,nnode),U(ndofel),
     * AMATRX(ndofel,ndofel),RHS(8,1),
     * xi, eta, GAUSS(nnode,2),det_J,WEIGHT(nnode),
     * Bu(3,2*nnode),Bu_T(2*nnode,3),D(3,3),UU(2*nnode,1),STRAIN(3,1),
     * STRESS(3,1),Ru(2*nnode),Kuu(8,8),Kuu11(8,3),Kuu12(8,8),Kuu1(8,8),
     * Ru1(8,1), Ex, Ey, G, nu, thck
  
      ! Initialization
      AMATRX=0.d0; RHS=0.d0; det_J=0.d0; Bu=0.d0; Bu_T=0.d0; D=0.d0
      UU=0.d0; STRAIN=0.d0; STRESS=0.d0; Ru=0.d0; Kuu=0.d0; Kuu11=0.d0
      Kuu12=0.d0; Kuu1=0.d0
! ______________________________________________________________________
! **************************************************************
! *  Constructing elemet TYPE U1( Displacement & Phase-field ) *
! **************************************************************
      
      Ex=PROPS(1);
      Ey=PROPS(2);
      G=PROPS(3);
      nu=PROPS(4);
      thck=PROPS(5)                 ! Material parameters
      
      
      CALL CONSTITUTIVE_MATRX(Ex,Ey,G,nu,D)                        ! Materials stiffness matrix
      
      
      DO i =1,nnode
          UU(2*i-1,1)=U(2*i-1); UU(2*i,1)=U(2*i)              ! DOF of displacement
      END DO
      
      CALL GAUSS_IPT(GAUSS,WEIGHT)                            ! Caussian point and its weights

      ! Calculating properties at each integral point
      LOOP_IPT :DO ipt=1,4
      ! ------------------------------< LOOP:Begin >----------------------------------
          xi  = GAUSS(ipt,1); eta = GAUSS(ipt,2)              ! Local coordinates of ipt

          CALL BASIC_MATRX(Bu,det_J,xi,eta,COORDS)      ! Sshape functions and its local derivatives
          
          CALL UMATMUL(Bu,UU,STRAIN,3,1,8)                    ! Local Strains & Stresses
          CALL UMATMUL(D,STRAIN,STRESS,3,1,3)
          
          ! Element stiffness matrix
          ! --------------------- Kuu ----------------------
          CALL UTRANSPOSE(Bu,Bu_T,3,8)
          CALL UMATMUL(Bu_T,D,Kuu11,8,3,3)
          CALL UMATMUL(Kuu11,Bu,Kuu12,8,8,3)
          Kuu1 = Kuu12*det_J*WEIGHT(ipt)*WEIGHT(ipt)*thck
          Kuu = Kuu+Kuu1

          ! Internal forces
          ! --------------------- Ru -----------------------
          CALL UMATMUL(Kuu1,UU,Ru1,8,1,8)
          DO i=1,8
              Ru(i)=Ru(i)-Ru1(i,1)
          END DO
          
      ! ------------------------------< LOOP:End >------------------------------------ 
      END DO LOOP_IPT
      ! Updatding AMATRX & RHS
      DO i=1,4
        DO j=1,4
          AMATRX(2*i-1,2*j-1)=Kuu(2*i-1,2*j-1);
          AMATRX(2*i,2*j-1)=Kuu(2*i,2*j-1)
          AMATRX(2*i-1,2*j)=Kuu(2*i-1,2*j)
          AMATRX(2*i,2*j)=Kuu(2*i,2*j)
        END DO
      END DO
      
      DO i=1,nnode
          RHS(2*i-1,1)=Ru(2*i-1)
          RHS(2*i,1)=Ru(2*i)
      END DO
      
      END SUBROUTINE ISOTROPIC

      
      
      
      
! Other subroutines might be used:
! < 1 >
!______________________________________________________________________
      SUBROUTINE GAUSS_IPT(GAUSS,WEIGHT)
      IMPLICIT NONE
      REAL(8) GAUSS(4,2),WEIGHT(4), x
      x = 0.577350269189626
      GAUSS(1,1)=-x; GAUSS(1,2)=-x
      GAUSS(2,1)= x; GAUSS(2,2)=-x
      GAUSS(3,1)= x; GAUSS(3,2)= x
      GAUSS(4,1)=-x; GAUSS(4,2)= x
      WEIGHT=1.d0
      END SUBROUTINE GAUSS_IPT
!______________________________________________________________________
      
! < 2 >
!______________________________________________________________________
      SUBROUTINE BASIC_MATRX(Bu,det_J,xi,eta,COORDS)
! Calculate the shape function and its derivatives
      IMPLICIT NONE
      INTEGER i,j
      REAL(8) N(1,4),dNdx(2,4),dNdxi(2,4),COORDS(2,4)
      REAL(8) JACOBI(2,2),INV_JACOBI(2,2), det_J
      REAL(8) Bu(3,8),Bp(2,4)
      REAL(8) xi,eta
      ! Initialization
      N=0.d0; dNdx=0.d0; dNdxi=0.d0; JACOBI=0.d0
      INV_JACOBI=0.d0; det_J=0.d0; Bu=0.d0; Bp=0.d0
      ! Values of shape functions as a function of local coord.
      N(1,1) = 1.d0/4.d0*(1.d0-xi)*(1.d0-eta)
      N(1,2) = 1.d0/4.d0*(1.d0+xi)*(1.d0-eta)
      N(1,3) = 1.d0/4.d0*(1.d0+xi)*(1.d0+eta)
      N(1,4) = 1.d0/4.d0*(1.d0-xi)*(1.d0+eta)
      ! Derivatives of shape functions respect to local coordinates
      dNdxi(1,1) = 1.d0/4.d0*(-1.d0+eta)
      dNdxi(1,2) = 1.d0/4.d0*(1.d0-eta)
      dNdxi(1,3) = 1.d0/4.d0*(1.d0+eta)
      dNdxi(1,4) = 1.d0/4.d0*(-1.d0-eta)
      dNdxi(2,1) = 1.d0/4.d0*(-1.d0+xi)
      dNdxi(2,2) = 1.d0/4.d0*(-1.d0-xi)
      dNdxi(2,3) = 1.d0/4.d0*(1.d0+xi)
      dNdxi(2,4) = 1.d0/4.d0*(1.d0-xi)
      ! Jacobi matrix
      DO j=1,4
          JACOBI(1,1)=JACOBI(1,1)+dNdxi(1,j)*COORDS(1,j)
          JACOBI(1,2)=JACOBI(1,2)+dNdxi(1,j)*COORDS(2,j)
          JACOBI(2,1)=JACOBI(2,1)+dNdxi(2,j)*COORDS(1,j)
          JACOBI(2,2)=JACOBI(2,2)+dNdxi(2,j)*COORDS(2,j)
      END DO
      ! Jacobi determinant
      det_J=JACOBI(1,1)*JACOBI(2,2)-JACOBI(1,2)*JACOBI(2,1)
      ! Inverse of JACOBI
      INV_JACOBI(1,1)= JACOBI(2,2)/det_J
      INV_JACOBI(1,2)=-JACOBI(1,2)/det_J
      INV_JACOBI(2,1)=-JACOBI(2,1)/det_J
      INV_JACOBI(2,2)= JACOBI(1,1)/det_J
      ! dNdx
      CALL UMATMUL(INV_JACOBI,dNdxi,dNdx,2,4,2)
      ! Bp & Bu matrix( Bp=dNdx )
      DO i=1,4
          Bp(1,i)=dNdx(1,i);      Bp(2,i)=dNdx(2,i)
          Bu(1,i*2-1)=dNdx(1,i);  Bu(2,i*2)  =dNdx(2,i)
          Bu(3,i*2-1)=dNdx(2,i);  Bu(3,i*2)  =dNdx(1,i)
      END DO
      END SUBROUTINE BASIC_MATRX
!______________________________________________________________________
  
! < 3 >
!______________________________________________________________________
      SUBROUTINE CONSTITUTIVE_MATRX(Ex,Ey,G,nu,D)
      IMPLICIT NONE
      REAL(8) Ex,Ey, G, nu, D(3,3)
      D=0.d0
      ! Planae stress
      D(1,1)=Ex/(1.d0-nu**2)
      D(2,2)=Ey/(1.d0-nu**2)
      D(3,3)=G
      D(1,2)=Ex*nu/(1.d0-nu**2)
      D(2,1)=Ey*nu/(1.d0-nu**2)
      
      ! Plane strain
      !D(1,1)=Ex/((1.d0+nu)*(1.d0-2.d0*nu))*(1.d0-nu)
      !D(2,2)=Ey/((1.d0+nu)*(1.d0-2.d0*nu))*(1.d0-nu)
      !D(3,3)=G
      !D(1,2)=Ex/((1.d0+nu)*(1.d0-2.d0*nu))*nu
      !D(2,1)=Ey/((1.d0+nu)*(1.d0-2.d0*nu))*nu
      END SUBROUTINE CONSTITUTIVE_MATRX
!______________________________________________________________________ 

      
! < 4 >
!______________________________________________________________________      
      SUBROUTINE UMATMUL(A,B,C,m,n,l)
! User defined matrix multiplication
      IMPLICIT NONE
      INTEGER i,j,k,m,n,l
      REAL(8) A(m,l),B(l,n),C(m,n)
      DO i=1,m
        DO j=1,n
            C(i,j)=0.d0
            DO k=1,l
                C(i,j)=C(i,j)+A(i,k)*B(k,j)
            END DO
        END DO
      END DO
      END SUBROUTINE UMATMUL
!______________________________________________________________________

! < 8 >
!______________________________________________________________________
      SUBROUTINE UTRANSPOSE(A,A_T,m,n)                                
! User defined matrix transpose
      IMPLICIT NONE
      INTEGER i,j,m,n
      REAL(8) A(m,n),A_T(n,m)
      DO i=1,m
          DO j=1,n
              A_T(j,i)=A(i,j)
          END DO
      END DO
      END SUBROUTINE UTRANSPOSE
!______________________________________________________________________