      subroutine strl81ld_matlab(PROPS,STRAN,STRSG,
     1                EQPST,QFUNR,
     2                FLOAD,FACTM,FACTP,FACTQ,
     3                HPACT,
     4                CPLAS,VPTRIA,VP,DIRP1,DIRP2,FG33)
C    SUBROUTINE STRL81_matlab(PROPS,STRAN,STRSG,
C   .                  FGRAD,EQPST,QFUNR,EQBIF,
C   .                  FLOAD, FACTM, FACTP,FACTQ,
C   .                  HPCON, 
C   .                  CPLAS, VPTRIA,VP,DIRP1, DIRP2,FG33)    
C************************************************************************
C
C     THIS ROUTINE EVALUATES STRESSES FOR ELASTO-PLASTIC  MODEL (J2)
C      (STRONG DISCONTINUITIES with variable Bandwidth) 
C                       -  PLANE STRESS  (ntype=1) -
C                       -  PLANE STRAIN  (ntype=2) -
C     ===========      LARGE  DEFORMATION    ===========
C                            NCRIT=81  
C************************************************************************
C
!    INPUT VARIABLES :
C
C       PROPS = Properties     
C       STRAN = Strain tensor, (Right Cauchy-Green Tensor b)     
C       ANGLN = Strong discontinuity angles   
C
!    OUTPUT VARIABLES :
C
C       STRAN = Strain tensor
C       EHOLD = Historical variables at time t_n
C       EHIST = Historical variables at time t_n+1
C       ESOLD = Strong discont. variables at time t_n
C       ESCUR = Strong discont. variables at time t_n+1
C       ANGLN = Strong discont. angles at time t_n+1
C       PROPS = Material properties      
C       PRNFL = Printing flag      
C
!    OUTPUT VARIABLES :
C
C       STRSG = Stress tensor
C       DMATX = Mechanical tangent matrix
C       EHIST = Historical variables at time t_n+1
C       ESCUR = Strong discont. variables at time t_n+1
C
C       STOLD = previous strain (internal forces -DIM: 4)
C             = auxiliary tangent matrix (tangent stiffness - DIM: 4x4)
C       ICRIB  = criteria for evaluating the bandwidth variation
C                = 1 : h(q) linear
C                = 2 : h = H_crit / H_bar
C
C************************************************************************
C
      IMPLICIT NONE
     
C     Constantes de dimensionamiento de arrays
      include 'ArraysDim.om'
C
      REAL*8  PROPS(NPROP),STRAN(2,2),STRSG(NSTR1)
C      REAL*8  FGRAD(2,*)
      REAL*8  FGRAD(2,2)
      REAL*8  DMATX(NSTR1,NSTR1),BEFIN(2,2),FGBIF(2,2)
      REAL*8  DERQ,DERH,AAUX(6),HCPRE
C
      REAL*8 YOUNG,POISS
      REAL*8 YIELD,GFVAL,HPOUT,
     .       EQPST,QFUNR,EQPRE,QFPRE,FLOAD,EQBIF,QBIFU,
     .       FACTM,FACTP,FACTQ,QFSTD,HPACT,FG33,FGT33,
     .       BETA,HB,MULAM,LAMBL,STRIA(4),CPPRE(4)
      REAL*8 CPLAS(4), VPTRIA(3),VP(3),DIRP1(3), DIRP2(3),CPINV(2,2),
     .       FTRAN(2,2),BFTRI(2,2),TRADE,TAUME,DEVST(3),SDTRE,DETCP
      REAL*8 COE1,COE2,DTRES,FTRIAL,KAPPA,TRZTS,TSPTS,DETFG,YFUNC,
     .       tcauc,HPPRE,CSD(2),DESS1,DESS2,T1,T2     
C
      DTRES = 2.D0/3.D0
      SDTRE = DSQRT(DTRES)
C
C***Recover properties 
C
      YOUNG = PROPS(1)                      ! YOUNG MODULUS
      POISS = PROPS(2)                      ! POISSON RATIO
      YIELD = PROPS(3)                      ! ULTIMATE STRESS (sigma_ult)
      HB = PROPS(4)                      ! SOFTENING PARAMETER (H_y)
C     MULAN: SHEAR MODULUS 
      MULAM=YOUNG/2.D0/(1.D0+POISS)
C     LAMBL: LAMÉ'S FIRST PARAMETER
      LAMBL=YOUNG*POISS/(1.D0+POISS)/(1.D0-2.D0*POISS)
C      IF(NTYPE.EQ.1)LAMBL=2.D0*MULAM*LAMBL/(LAMBL+2.D0*MULAM)
C     KAPPA: BULK MODULUS
      KAPPA=LAMBL+DTRES*MULAM
C     En HPACT se guarda el módulo de endurecimiento regularizado (todavía no calculado) 
      HPACT = HB

C
C***  INITIALIZE VARIABLES OF THE PRESENT STEP for bandwidth integration
C-----------------------------------------------------------------------
C
      IF(FACTQ.LE.0.D0)FACTQ=1.D0
C
      EQPRE=EQPST
      QFPRE=QFUNR
C      BWPRE=BANDW
      FGT33=FG33
      CALL VECASI(NSTR1,CPLAS,CPPRE)        !   CPPRE = CPLAS
C
C***  Update defor. gradient
C----------------------------------------------------------------------
      CALL VECASI(NSTR1,STRAN,FGRAD)        !   FGRAD = STRAN
C 
      DETFG = (FGRAD(1,1)*FGRAD(2,2)-FGRAD(1,2)*FGRAD(2,1))*FG33
C test necesario para saber si DETFG > 0
C***  evaluate the trial Right-Cauchy-Green tensor (b_trial)
C
      DETCP=(CPLAS(1)*CPLAS(2)-CPLAS(3)*CPLAS(3))

      CPINV(1,1)= CPLAS(2)/DETCP
      CPINV(2,2)= CPLAS(1)/DETCP
      CPINV(1,2)=-CPLAS(3)/DETCP
      CPINV(2,1)=-CPLAS(3)/DETCP
C      
      FTRAN(1,1)= FGRAD(1,1)
      FTRAN(1,2)= FGRAD(2,1)
      FTRAN(2,1)= FGRAD(1,2)
      FTRAN(2,2)= FGRAD(2,2)
C
      VPTRIA(3) = FG33*FG33/CPLAS(4)
C
C**** Evaluation of b^{trial} = F C_p^{-1} F^T
      CALL FCNVFT(FGRAD,CPINV,BFTRI)
      CALL DPPAL0(BFTRI,VPTRIA, DIRP1,DIRP2 )
      CALL VECASI(3,VPTRIA,VP)  !  VP    = VPTRIA
C
C***Evaluation of SIGMA_trial
C
      TRADE=(DLOG(VPTRIA(1)*VPTRIA(2)*VPTRIA(3)))/2.D0
      TAUME=KAPPA*TRADE
C     Direcciones del desviador del Tau de trial      
      DEVST(1)=MULAM*(DLOG(VPTRIA(1)) - DTRES*TRADE)
      DEVST(2)=MULAM*(DLOG(VPTRIA(2)) - DTRES*TRADE)
      DEVST(3)=MULAM*(DLOG(VPTRIA(3)) - DTRES*TRADE)
C     Sigma de truak
      STRSG(1)=DEVST(1)*DIRP1(1)+DEVST(2)*DIRP2(1)+TAUME
      STRSG(2)=DEVST(1)*DIRP1(2)+DEVST(2)*DIRP2(2)+TAUME
      STRSG(3)=DEVST(1)*DIRP1(3)+DEVST(2)*DIRP2(3)
      STRSG(4)=DEVST(3)+TAUME
C  Evaluation of yield function (trial state)
C ------------------------------------------
      FTRIAL=DSQRT(DEVST(1)*DEVST(1)+DEVST(2)*DEVST(2)+
     .     DEVST(3)*DEVST(3)) - SDTRE*(YIELD-QFUNR)
C
C***  assume elastic  or unloading condition 
C     --------------------------------------
      FLOAD=0.D0
      FACTM=0.D0
      FACTQ=0.D0
C
      CALL VECASI(NSTR1,STRSG,STRIA)        !   STRIA = STRSG
C
      IF(FTRIAL.GT.0.D0)THEN
C***  loading condition
C     -----------------
         FLOAD=1.D0
C        evaluate the stresses and internal variables
C        --------------------------------------------
         CALL PCJ81(STRSG,TAUME,DEVST,VP,CPLAS,
     .     PROPS,EQPST,QFUNR,FTRIAL,TRADE,FTRAN,
     .     HPACT,FLOAD,FACTM,FACTP,FACTQ,DIRP1,DIRP2,FG33)
C
      ENDIF
C
      RETURN
      END
C
C
C************************************************************************
      SUBROUTINE FCNVFT(F,C,B)
C*********************************************************************
C    CALCULA   B= F C F^T   MATRICES DE [2X2]
C********************************************************************
      IMPLICIT NONE
C
      REAL*8 B(2,2),F(2,2),C(2,2), R11,R12,R21,R22
C
      R11=F(1,1)*C(1,1) + F(1,2)*C(2,1)
      R12=F(1,1)*C(1,2) + F(1,2)*C(2,2)
      R21=F(2,1)*C(1,1) + F(2,2)*C(2,1)
      R22=F(2,1)*C(1,2) + F(2,2)*C(2,2)
C
      B(1,1)= R11*F(1,1)+R12*F(1,2)
      B(1,2)= R11*F(2,1)+R12*F(2,2)
      B(2,1)= R21*F(1,1)+R22*F(1,2)
      B(2,2)= R21*F(2,1)+R22*F(2,2)
C
      RETURN
      END
C
C
C
C************************************************************************
      SUBROUTINE DPPAL0(T,VP, M1,M2 )
C************************************************************************
C
C    THIS ROUTINE EVALUATE THE PRINCIPAL DIRECTION AND EIGENVALUES
C              OF A SYMMETRIC TENSOR (2,2)
C************************************************************************
C
!    INPUT VARIABLES :
C
C       T  =  SYMMETRIC TENSOR (2,2)    
C       VP =  VECTOR OF EIGENVALUES (3) 
C       M1 =  TENSOR OF PPAL. DIRECT.(N1 TIMES N1) STORED AS VECT(3)
C       M2 =  TENSOR OF PPAL. DIRECT.(N2 TIMES N2) STORED AS VECT(3) 
C************************************************************************
C
      IMPLICIT NONE
C
      REAL*8 T(2,2),VP(3),M1(3),M2(3),THETA,N1(2),N2(2),R11,R12,R21,R22 
C
      IF(T(1,2).eq.0.d0.and.(T(1,1)-T(2,2)).eq.0.d0) then
         THETA = 0.D0
      ELSE
         THETA = ATAN2(2*T(1,2) , (T(1,1)-T(2,2))) / 2.D0
      ENDIF
C
      N1(1) = DCOS(theta)
      N1(2) = DSIN(theta)
      N2(1) = -N1(2)
      N2(2) =  N1(1)

      R11 = N1(1)
      R12 = N2(1)
      R21 = N1(2)
      R22 = N2(2)

      VP(1)=R11*(T(1,1)*R11+T(1,2)*R21)+R21*(T(2,1)*R11+T(2,2)*R21)
      VP(2)=R12*(T(1,1)*R12+T(1,2)*R22)+R22*(T(2,1)*R12+T(2,2)*R22)

      M1(1)= N1(1)*N1(1)
      M1(2)= N1(2)*N1(2)
      M1(3)= N1(1)*N1(2)

      M2(1)= N2(1)*N2(1)
      M2(2)= N2(2)*N2(2)
      M2(3)= N2(1)*N2(2)
C
      RETURN
      END
C
C************************************************************************
      SUBROUTINE PCJ81(STRSG,TAUME,DEVST,VP,CPLAS,
     .           PROPS,EQPST,QFUNR,FTRIAL,TRADE,FTRAN,
     .           HPACT,FLOAD,FACTM,FACTP,FACTQ,DIRP1,DIRP2,FG33)
C************************************************************************
C    THIS FUNCTION EVALUATES THE STRES-STRAIN RELATION.
C    CONSTITUTIVE MODEL: Elasto-Plastic J2 Theory
C************************************************************************
C
C!    INPUT VARIABLES :
C
C       TAUME = Mean stress 
C       DEVST = Principal Trial deviatoric Stress (3)
C       EQPST = Int.var.like strain(Theor.J2)-Equiv.plast.strain::alpha)
C       QFUNR = Scalar variable like stress  (q(alpha))
C       EQBIF = Variable (alpha) at  bifurcation.
C       FLOAD = Load/Unload flag   
C       FACTM = A factor to be used in the tangent matrix evaluation
C       FACTP = A factor to be used in the tangent matrix evaluation
C       FACTQ = A factor to be used in the tangent matrix evaluation
C
C    Output variables: 
C
C       STRSG(NSTR1)  Stresses at Gauss point NGAUL
C       CPLAS(NSTR1)  Plastic strains at Gauss point NGAUL
C       LESLT(3)      Elastic principal stretches (tensor b)
C************************************************************************
C
      IMPLICIT NONE
C
C    INCLUDE 'propert.om'  ! -->  PH,PA

      REAL*8 STRSG(*),TAUME,DEVST(*),VP(*),CPLAS(*),
     .       PROPS(*),EQPST,QFUNR,HPACT,FLOAD,
     .       FACTM,FACTP,FACTQ,TRADE,DEPSI(3),BFIN(2,2),FTRAN(2,2),
     .       YIELD,DGAMA,FTRIAL,
     .       DTRES,SDTRE,MULAM,LAMBL,DOSMU,KAPPA,MSTRI,COE
      REAL*8 DETBE,BEINV(2,2),CPAUX(2,2),DIRP1(*),DIRP2(*)   !,QFAUX
      REAL*8 YOUNG,POISS
      REAL*8 HB,TITA,KINFB,KCEROB,DELTA
      
      REAL*8 alpha_new,norm_gfunction,toldg
      REAL*8 k_new,kp_new,gfunction,dgfunction 
cccccccccccccccccccccccccccccccccc
      REAL*8 DETFG, DETCP, CPINV(2,2), FG33
      REAL*8 BFTRI(2,2),VPTRIA(3), STRG_aux(4),norma_dif
c
      INTEGER  iterco,lineal,iter
C
      DTRES = 2.D0/3.D0
      SDTRE = DSQRT(DTRES)      
C
      YOUNG = PROPS(1)                      ! YOUNG MODULUS
      POISS = PROPS(2)                      ! POISSON RATIO
      YIELD = PROPS(3)                      ! ULTIMATE STRESS (sigma_ult)
C     Se usa el HB que viene como argumento para considerar que puede venir regularizado con el espesor del elemento      
C     HB = PROPS(4)
      TITA = PROPS(5)
      KINFB = PROPS(6)
      KCEROB = PROPS(7)
      DELTA = PROPS(8)
C     MULAN: SHEAR MODULUS 
      MULAM=YOUNG/(2.D0*(1.D0+POISS))
C     LAMBL: LAMÉ'S FIRST PARAMETER
      LAMBL=YOUNG*POISS/(1.D0+POISS)/(1.D0-2.D0*POISS)
      DOSMU=2.D0*MULAM
C     KAPPA: BULK MODULUS
      KAPPA=LAMBL+DTRES*MULAM
      
C     Para simplificar el cálculo se determina si el modelo es lineal.  
      IF(KINFB.NE.KCEROB .AND. DELTA.NE.0.D0) THEN
C        Exponential-Linear Hardening
         lineal = 0
      ELSE
C        Linear Hardenning      
         lineal = 1
      ENDIF
      
      IF(QFUNR.GE.YIELD)THEN
         QFUNR=YIELD
         STRSG(1)=TAUME                         
         STRSG(2)=TAUME
         STRSG(4)=TAUME
         STRSG(3)=0.d0
         
         VP(1)=DEXP(DTRES*TRADE)
         VP(2)=VP(1)
         VP(3)=VP(1)
         BFIN(1,1)=VP(1)*(DIRP1(1)+DIRP2(1))
         BFIN(2,2)=VP(1)*(DIRP1(2)+DIRP2(2))
C     
         BEINV(1,1)= 1.D0/BFIN(1,1)
         BEINV(2,2)= 1.D0/BFIN(2,2)
         BEINV(1,2)= 0.D0
         BEINV(2,1)= 0.D0
C     
C     Evaluation of C_p  = F^T B^{-1} F
C     ---------------------------------
         CALL FCNVFT(FTRAN,BEINV,CPAUX)
C     
         CPLAS(1)=CPAUX(1,1)
         CPLAS(2)=CPAUX(2,2)
         CPLAS(3)=CPAUX(1,2)
         CPLAS(4)=FG33*FG33 /VP(3)
         
         factm=0.d0
         factq=0.d0
         RETURN
      ENDIF
C     
C     Norma del tensor desviador
      MSTRI=DSQRT(DEVST(1)*DEVST(1)+DEVST(2)*DEVST(2)+
     .      DEVST(3)*DEVST(3))
      IF(MSTRI.EQ.0.D0) THEN
         MSTRI=1.D0
      ENDIF
C    
      if(lineal.EQ.1) then
      
         DGAMA=FTRIAL/(DTRES*HPACT+DOSMU) ! Incr. of gamma
C     Parece que no se usa para nada QFAUX
C        QFAUX=QFUNR
         QFUNR=QFUNR-SDTRE*HPACT*DGAMA                 
      
C        FACTM = 1.D0/(1.D0+HPACT/3.D0/MULAM)

      else   

         norm_gfunction = 1.d00 
         toldg  = 1.d-6*YIELD    
      
         dgama = 0.d0 
         alpha_new = EQPST
         iter  = 0

         do while (norm_gfunction .gt. toldg .and. iter.lt. 100)
            iter=iter+1
            QFUNR = -((TITA*HPACT*alpha_new)+(KINFB-KCEROB)*
     1         (1.-dexp(-DELTA*alpha_new)))
            k_new = YIELD - QFUNR
            kp_new = (TITA*HPACT)+(KINFB-KCEROB)*
     1         (DELTA*dexp(-DELTA*alpha_new))
            gfunction = -SDTRE*k_new + MSTRI - 2*MULAM*dgama
            dgfunction = -2*MULAM*(1+kp_new/3/MULAM)
            dgama = dgama - gfunction/dgfunction
            alpha_new = EQPST + SDTRE*dgama
            norm_gfunction = dabs(gfunction)
         end do
         if (iter. ge.99) stop 'iteraciones dgamma > 100'
         HPACT = TITA*HPACT+(KINFB-KCEROB)*DELTA*dexp(-DELTA*alpha_new)
          
      endif
      FACTQ = DOSMU*DGAMA/MSTRI
      FACTM = 1.D0/(1.D0+HPACT/3.D0/MULAM)
c     ----------------------
C     contraction of S_trial
C     ----------------------
      COE = 1.D0-FACTQ      
C
      STRSG(1)=COE*(DEVST(1)*DIRP1(1)+DEVST(2)*DIRP2(1)) + TAUME
      STRSG(2)=COE*(DEVST(1)*DIRP1(2)+DEVST(2)*DIRP2(2)) + TAUME
      STRSG(3)=COE*(DEVST(1)*DIRP1(3)+DEVST(2)*DIRP2(3))
      STRSG(4)=COE*(DEVST(3)) + TAUME
C     -------------------------
C     update equiv.plast.strain
C     -------------------------
      EQPST=EQPST+SDTRE*DGAMA
C     ---------------------------
C     update plastic strain CPLAS
C     ---------------------------
      DEPSI(1)=(COE*DEVST(1))/MULAM
      DEPSI(2)=(COE*DEVST(2))/MULAM
      DEPSI(3)=(COE*DEVST(3))/MULAM
C     
      VP(1)=DEXP((DEPSI(1)+DTRES*TRADE))
      VP(2)=DEXP((DEPSI(2)+DTRES*TRADE))
      VP(3)=DEXP((DEPSI(3)+DTRES*TRADE))
C     
      BFIN(1,1)=VP(1)*DIRP1(1)+VP(2)*DIRP2(1)
      BFIN(2,2)=VP(1)*DIRP1(2)+VP(2)*DIRP2(2)
      BFIN(1,2)=VP(1)*DIRP1(3)+VP(2)*DIRP2(3)
      BFIN(2,1)=BFIN(1,2)
C     c
      DETBE=BFIN(1,1)*BFIN(2,2)-BFIN(1,2)*BFIN(2,1)
      BEINV(1,1)= BFIN(2,2)/DETBE
      BEINV(2,2)= BFIN(1,1)/DETBE
      BEINV(1,2)=-BFIN(1,2)/DETBE
      BEINV(2,1)= BEINV(1,2)
C     ---------------------------------
C     Evaluation of C_p  = F^T B^{-1} F
C     ---------------------------------
      CALL FCNVFT(FTRAN,BEINV,CPAUX)
C     
      CPLAS(1)=CPAUX(1,1)
      CPLAS(2)=CPAUX(2,2)
      CPLAS(3)=CPAUX(1,2)
      CPLAS(4)= FG33*FG33 / VP(3)
C     
      RETURN
      END
C
C
C************************************************************************
      SUBROUTINE MTDL81(STRAN,PROPS,DMATX,FLOAD,FACTM,FACTP,FACTQ,
     .                  VPTRIA,VP,MA,MB)
C************************************************************************
C
C    THIS ROUTINE EVALUATES MATRIX FOR DAMAGE MODELS
C                   (STRONG DISCONTINUITIES)
C
C************************************************************************
C
!    INPUT VARIABLES :
C
C       STRAN = Strain tensor     
C       RINTV = Internal varible describing the size of the 
C               elastic region   
C       FLOAD = Load/Unload flag   
C       FACTM = A factor to be used in the tangent matrix evaluation
C
!    OUTPUT VARIABLES :
C
C       DMATX = Constitutive matrix
C
C       DGAT2 = Geometric additional matrix (present constitutive model)
C************************************************************************
C
      IMPLICIT NONE
C
C     INCLUDE 'element.om'  ! -->  NHIST, NSTR1
C     INCLUDE 'propert.om'  ! -->  PH,PA
C     INCLUDE 'interval.om' ! -->  DTIME
      include 'ArraysDim.om'
C
      REAL*8  DMATX(NSTR1,NSTR1), STRAN(2,2), PROPS(*)
C      
      REAL*8  FLOAD,FACTM,FACTP,FACTQ,FP,E1,E2,E3,
     .        VPTRIA(*),VP(*),MA(*),MB(*),SIG1,SIG2,SIG3,SIGI(2,2),
     .        MULAM,LAMBL,DOSMU,KAPPA,TP(3),TREPS,DTRES,UTRES,
     .        TE1,TE2,
     .        H1,H2,H3,H4,H5,H6,H11,H22,H21,H31,H51,B1,B2,B3,
     .        MA1,MA2,MA3,MB1,MB2,MB3, MODDES,FT1,FT2,FT3, SN1,DN1,
     .        DSD,DSU,DD,V1,V2,V3,VP1,VP2,VP3,PERT,TKMED,
     .        C1(6),C2(6),C3(6),D(6),DETB,COE1,COE2,COE3,COE4,
     .        THETA,N11,N12,N21,N22,POIPS,PSIG(2),VEC(2),SIPSI
      REAL*8 YOUNG,POISS
C
      INTEGER KLOAD,I      
C
      DTRES=2.D0/3.D0
      UTRES=1.D0/3.D0
      
      YOUNG=PROPS(1)                      ! YOUNG MODULUS
      POISS=PROPS(2)                      ! POISSON RATIO

C     MULAN: SHEAR MODULUS      
      MULAM=YOUNG/(2.D0*(1.D0+POISS))
C     LAMBL: LAMÉ'S FIRST PARAMETER  
      LAMBL=YOUNG*POISS/(1.D0+POISS)/(1.D0-2.D0*POISS)
C     Término que se utiliza en las operaciones siguientes
      DOSMU=2.D0*MULAM
C     KAPPA: BULK MODULUS
      KAPPA=LAMBL+2.D0/3.D0*MULAM
C
C------------------------   perturbar VPTRIAL  -----------------
      pert = 1.D-7*DSQRT(VPTRIA(1)*VPTRIA(1)+VPTRIA(2)*VPTRIA(2)+
     .                   VPTRIA(3)*VPTRIA(3))  
      V1 = VPTRIA(1)
      V2 = VPTRIA(2)
      V3 = VPTRIA(3)  
      IF(Dabs(V1-V2).LT.pert)THEN
         V1= (V1+V2)/2.D0+pert 
         V2= (V1+V2)/2.D0-pert
      ENDIF
      IF(Dabs(V1-V3).LT.pert)THEN
         V1 = (V1+V3)/2.D0+pert
      ENDIF
      IF(Dabs(V3-V2).LT.pert)THEN
         V2 = (V3+V2)/2.D0-pert
      ENDIF
C C------------------------   perturbar VP  -----  -----------------
C     Se perturba VP para el caso que los autovalores son nulos, y el logaritmo no se haga -Inf.
      pert = 1.D-7*DSQRT(VP(1)*VP(1)+VP(2)*VP(2)+
     .                   VP(3)*VP(3))  
      VP1 = VP(1)
      VP2 = VP(2)
      VP3 = VP(3)  
      IF(Dabs(VP1-VP2).LT. pert)THEN
         VP1= (VP1+VP2)/2.D0+pert 
         VP2= (VP1+VP2)/2.D0-pert
      ENDIF
      if(Dabs(VP1-VP3).LT. pert)THEN
          VP1= (VP1+VP3)/2.D0+pert
      ENDIF
      if(Dabs(VP3-VP2).LT. pert)THEN
          VP2= (VP3+VP2)/2.D0-pert
      ENDIF
C C--------------------------------------------------------------- 
      E1=DLOG(VP1)
      E2=DLOG(VP2)
      E3=DLOG(VP3)
C     Se está calculando las tensiones desviatóricas (no trial).
      TREPS = E1+E2+E3
      TKMED = KAPPA*TREPS/2.D0
C     Valores principales del desviador de las tensiones Tau (no trial).
      TP(1) = MULAM*(E1-UTRES*TREPS)
      TP(2) = MULAM*(E2-UTRES*TREPS)
      TP(3) = MULAM*(E3-UTRES*TREPS)
C     Direcciones principales (corresponden al tensor V)
      MA1 = MA(1)
      MA2 = MA(3)
      MA3 = MA(2)
      MB1 = MB(1)
      MB2 = MB(3)
      MB3 = MB(2)
C------------------------   perturbar VP  -----  -----------------
      IF(MA1.EQ.1D0.OR.MA3.EQ.1.D0) THEN
         THETA=1.D-6
         N11 = DCOS(theta)
         N12 = DSIN(theta)
         N21 = -N12
         N22 =  N11
         MA1= N11*N11
         MA3= N12*N12
         MA2= N11*N12
         MB1= N21*N21
         MB3= N22*N22
         MB2= N21*N22
      ENDIF
C-----------------------------------------------------------------
      B1 = V1*MA1+V2*MB1
      B3 = V1*MA3+V2*MB3
      B2 = V1*MA2+V2*MB2

C      KLOAD=INT(FLOAD) 
C       IF(KLOAD .EQ.0.D0) THEN
C C        IF(NTYPE.EQ.1)THEN
C C           H1=LAMBL+DOSMU
C C           H2=LAMBL
C C           H4=H1
C C        ELSE
C             H1 = DOSMU 
C             H2 = 0.D0
C             H3 = 0.D0
C             H4 = DOSMU
C             H5 = 0.D0    
C             H6 = DOSMU
C C            
C             H11 = ( 2.D0 * H1 - H2 - H3 ) / 3.D0  
C             H21 = ( 2.D0 * H2 - H1 - H3 ) / 3.D0  
C             H22 = ( 2.D0 * H4 - H5 - H2 ) / 3.D0
C             H1 = H11+KAPPA      
C             H2 = H21+KAPPA       
C             H4 = H22+KAPPA
C C           Se agrega los directamente los valores del tensor tangente del modelo de Hencky.
C             H3 = KAPPA-DOSMU/3.D0
C             H5 = H3
C             H6 = KAPPA+2.D0*DOSMU/3.D0
C C        ENDIF
C       ELSE
C C        IF(NTYPE.EQ.1)THEN
C Ccc            poips=poiss/(1.d0+poiss)
C C           poips=poiss
C C           SIG1=3.D0*YOUNG/(3.D0*(1.D0-POIPS)+FACTQ*YOUNG)
C C           SIG2=DOSMU/(1.D0+DOSMU*FACTQ)
C C           SIGI(1,1)=0.5*(sig1+sig2)
C C           SIGI(2,1)=0.5D0*(sig1-sig2)
C C           SIGI(1,2)=SIGI(2,1)
C C           SIGI(2,2)=SIGI(1,1)
C C           PSIG(1)= DTRES*TP(1)-UTRES*TP(2)
C C           PSIG(2)=-UTRES*TP(1)+DTRES*TP(2)
C C           VEC(1)= SIGI(1,1)*PSIG(1)+SIGI(1,2)*PSIG(2)
C C           VEC(2)= SIGI(2,1)*PSIG(1)+SIGI(2,2)*PSIG(2)
C C           SIPSI=VEC(1)*PSIG(1)+VEC(2)*PSIG(2)
C C           COE1=1.D0/(SIPSI+FACTM)
C C
C C           H1=SIGI(1,1)-COE1*VEC(1)*VEC(1)
C C           H2=SIGI(1,2)-COE1*VEC(1)*VEC(2)
C C           H4=SIGI(2,2)-COE1*VEC(2)*VEC(2)
C C        ELSE
C     Dirección normal a la superficie de fluencia, como la norma del tensor desviador
C     de tensiones (en este caso se podría usar del tensor desviador de tensiones trial,
C     ya que tiene las mismas direcciones). 
      moddes = DSQRT(TP(1)*TP(1)+TP(2)*TP(2)+TP(3)*TP(3))   
      ft1=TP(1)/moddes
      ft2=TP(2)/moddes
      ft3=TP(3)/moddes
      TE1 = 1-FACTQ
      TE2 = FACTM-FACTQ
C             
      H1 = KAPPA+DOSMU*( TE1*DTRES-TE2*ft1*ft1)
      H2 = KAPPA+DOSMU*(-TE1*UTRES-TE2*ft1*ft2)
      H3 = KAPPA+DOSMU*(-TE1*UTRES-TE2*ft1*ft3)
      H4 = KAPPA+DOSMU*( TE1*DTRES-TE2*ft2*ft2)
      H5 = KAPPA+DOSMU*(-TE1*UTRES-TE2*ft2*ft3)
      H6 = KAPPA+DOSMU*( TE1*DTRES-TE2*ft3*ft3)
      
C     Se determina los valores principales de las tensiones totales (no trial)
      TP(1) = TP(1)+TKMED
      TP(2) = TP(2)+TKMED
      TP(3) = TP(3)+TKMED
C
C*** 
C
      C1(1) = H1*MA1 *MA1 + 2.D0*H2*MA1*MB1 + H4*MB1*MB1  

      C1(2) = H1*MA3*MA1 + H2*(MA1*MB3+MA3*MB1) +  H4*MB3*MB1  

      C1(3) = H1*MA1*MA2 + H2*(MA1*MB2+MA2*MB1) +  H4*MB2*MB1  

      C1(4) = H1*MA3*MA3 + 2.D0*H2*MA3*MB3 + H4*MB3*MB3 

      C1(5) = H1*MA3*MA2 + H2*(MA3*MB2 + MA2*MB3) +  H4*MB2*MB3  

      C1(6) = H1*MA2*MA2 + 2.D0*H2*MA2*MB2 + H4*MB2*MB2
C  
C***
C
      detB = B1*B3-B2*B2 

      COE1 = 1.D0/( ( V1 - V2)*(V1-V3) ) 
      COE2 = V2*V3 
      COE3 = V1 * COE1 
      COE4 = V2+V3-3.D0*V1 

      C2(1) = COE1* ( - COE2*(1.D0-(1.D0-MA1)*(1.D0-MA1) )) + 
     .         COE3* ( 2.D0*B1*MA1+ COE4*MA1*MA1 ) 

      C2(2) = COE1*(-detB-COE2*(-(1.D0-MA1)*(1.D0-MA3)))
     .       + COE3*((B1*MA3+B3*MA1) + COE4*(MA1*MA3) ) 

      C2(3) = COE1* ( - COE2*((1.D0-MA1)*MA2) ) + COE3*
     .     ((B1*MA2+B2*MA1) + COE4*(MA1*MA2) ) 

      C2(4) = COE1* ( - COE2*((1.D0- (1.-MA3)*(1.-MA3) ) )) + 
     .         COE3*(2.D0*B3*MA3 + COE4*MA3*MA3) 

      C2(5) = COE1* ( - COE2*((1.D0-MA3)*MA2) ) + COE3*
     .         ((B3*MA2+B2*MA3) + COE4*(MA3*MA2) ) 

      C2(6) = COE1* (0.5D0*detB - COE2*(0.5D0-MA2*MA2) ) 
     .     + COE3*(2.D0*B2*MA2 + COE4*MA2*MA2 ) 
C
C***
C
      COE1 = 1.D0/((V2-V1)*(V2-V3) ) 
      COE2 = V1*V3 
      COE3 = V2*COE1 
      COE4 = V1+V3-3.D0*V2 

      C3(1) = COE1* ( - COE2*(1.D0-(1.D0-MB1)*(1.D0-MB1)) ) + 
     .         COE3* ( 2.D0*B1*MB1+ COE4*MB1*MB1 )  

      C3(2) = COE1*(-detB-COE2*(-(1.D0-MB1)*(1.D0-MB3)))
     .       + COE3*((B1*MB3+B3*MB1) + COE4*(MB1*MB3) )  

      C3(3) = COE1* ( - COE2*((1.D0-MB1)*MB2) ) + COE3*
     .         ((B1*MB2+B2*MB1) + COE4*(MB1*MB2) )  

      C3(4) = COE1* ( - COE2*(1.D0-(1.-MB3)*(1.-MB3)) ) + 
     .         COE3* (2.D0*B3*MB3 + COE4*(MB3*MB3) )  

      C3(5) = COE1* ( - COE2*((1.D0-MB3)*MB2)) + COE3*
     .         ((B3*MB2+B2*MB3) + COE4*(MB3*MB2) )  

      C3(6) = COE1* (0.5D0*detB - COE2*(0.5D0-MB2*MB2)) 
     .       + COE3*(2.D0*B2*MB2 + COE4*(MB2*MB2) )  
C
      DO i=1,6
         D(i) = C1(i)+2.D0*(TP(1)*C2(i)+TP(2)*C3(i))
      ENDDO
      
      DMATX(1,1) = D(1)         
      DMATX(2,1) = D(2)        
      DMATX(3,1) = D(3)  
      DMATX(1,2) = DMATX(2,1)
      DMATX(2,2) = D(4)        
      DMATX(3,2) = D(5) 
      DMATX(1,3) = DMATX(3,1)  
      DMATX(2,3) = DMATX(3,2)  
      DMATX(3,3) = D(6)     
                
      DMATX(4,4)= H6-2.D0*TP(3)  
      DMATX(1,4)= H3*MA1+H5*MB1      
      DMATX(2,4)= H3*MA3+H5*MB3
      DMATX(3,4)= H3*MA2+H5*MB2 
      DMATX(4,1)= DMATX(1,4)      
      DMATX(4,2)= DMATX(2,4)  
      DMATX(4,3)= DMATX(3,4) 
 
C     Hay que corregir lo de abajo para considerar que Fzz no es nulo?
      !IF(KLOAD .EQ.1.D0) THEN
      !   if((factm.eq.0.d0.and.factq.eq.0.d0)) then
      !     DMATX(1,1) = kappa + dosmu*0.6666666666666d-5    
      !     DMATX(2,1) = kappa - dosmu*0.3333333333333d-5
      !     DMATX(3,1) = 0.d0
      !     DMATX(1,2) = kappa- dosmu*0.3333333333333d-5
      !     DMATX(2,2) = kappa + dosmu*0.6666666666666d-5
      !     DMATX(3,2) = 0.d0
      !     DMATX(1,3) = 0.d0 
      !     DMATX(2,3) = 0.d0 
      !     DMATX(3,3) = 0.d0 + dosmu*0.5d-5
      !                           
      !     DMATX(4,4)=  kappa 
      !     DMATX(1,4)=   0.d0
      !     DMATX(2,4)=   0.d0
      !     DMATX(3,4)=   0.d0
      !     DMATX(4,1)=   0.d0 
      !     DMATX(4,2)=   0.d0
      !     DMATX(4,3)=   0.d0
      !   endif
      !endif
C
      RETURN
C
      END
C
