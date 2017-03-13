      SUBROUTINE STRDVB_DO(PROPS,STRAN,STRSG,DMATX,
     .                  STREF,RINTV,QFUNR,MOD_BETA,QBIFU,dq_beta,
     .                  beta_s,FLOAD, FACTM, FACTP,
     .                  FACTQ,ESCUR,ESOLD,ANGLN,HPCRI,STOLD,
     .                  HPCON,
     .                  angl1,angl2,ltflg,angle0,VTRAC,tenspg,CARTD,first_t,
     .                  eps_el, ang_hist, var_g )
C************************************************************************
C
C     THIS ROUTINE EVALUATES STRESSES FOR DAMAGE  MODELS
C      (STRONG DISCONTINUITIES with variable Bandwidth)  
C         Modelo cohesivo discreto de danio oliver
C         NCRIT = 78    Cohesive Model Projected of the SD method
C 
C************************************************************************
C
!    INPUT VARIABLES : 
C
C       PROPS = Properties     
C       STRAN = Strain tensor      
C       ESCUR = LTRAK,NSOLI,BANDW,QBIFU(makes the funct.of a bif.time)
C       ANGLN = Strong discontinuity angles   
C
!    OUTPUT VARIABLES :
C
C       STRSG    = Stress tensor
C       STREF    = Efective Stress  
C       MOD_BETA = equivalent displacement of the crack 
C       QFUNR    = Scalar damage variable  (q(r))
C       QBIFU    = equivalent traction at the bifurcation time.
C       RBIFU    = Variable (r) at  bifurcation.
C       BWBIF    = Bandwidth (h) at bifurcation. 
C       FLOAD    = Load/Unload flag   
C       FACTM    = A factor to be used in the tangent matrix evaluation
C       FACTP    = A factor to be used in the tangent matrix evaluation
C       FACTQ    = A factor to be used in the tangent matrix evaluation
C       ESCUR    = LTRAK,NSOLI,BANDW
C       ANGLN    = Strong discontinuity angles   
C
!    MAIN VARIBLES USED INTO THIS ROUTINE
C
C       RTRIAL = trial value of r
C       RIPRE  = value of r at previous converged step
C       QFPRE  = value of q at previous converged step 
C
C************************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'problem.om'  ! -->  NPROP
      INCLUDE 'element.om'  ! -->  NHIST, NSTR1, NGAUL,NCRIT
      INCLUDE 'propert.om'  ! -->  PH,PA
      INCLUDE 'auxiliar.om' ! -->  IGAUL
      INCLUDE 'interval.om' ! -->  IITER,DTIME
      INCLUDE 'logical.om' ! -->  IITER,DTIME
      INCLUDE 'cstrong.om'                                   ! STDC
      INCLUDE 'pararcl.om'
      INCLUDE 'exfem.om'
      INCLUDE 'dynamic.om'
C
      INTEGER LTRAK,ICRIB,ILARG, implicito,istr,jstr,kstr,iter
      INTEGER LTFLG(NELEM),I,J,IANG,ltra2,LFLAG,nsol_b
      REAL*8  PROPS(NPROP), STRAN(NSTR1), STRSG(NSTR1)
      REAL*8  STREF(NSTR1), RINTV, DAMAG, HPBIF, TDAMG
      REAL*8  ESCUR(*), ESOLD(*),     ANGLN(NDIME,NDIME)
      REAL*8  DMATX(NSTR1,NSTR1), STPOS(6),BWBIF,QFUNR
      REAL*8  STOLD(NSTR1),AAUX(6),CARTD(NDIME,NNODL),VECDOT
C
      REAL*8  FTULT,GFVAL,HPCRI,DGPRE, dete,eps_salto(4),RIPRE,QFPRE,rtrial,
     .        FLOAD,RBIFU,QBIFU,MOD_BETA_PRE, MOD_BETA_TILDE, MOD_BETA,
     .        FACTM,FACTP,FACTQ,HPOUT,RINT0,ZEROQ,rtilde,qtilde
      REAL*8  hpaux,dummy,pimed,dq_beta,dr_pre,deltr
      REAL*8  ANGL1(NELEM,*),ANGL2(NELEM,5),HPCON
      REAL*8  daux,angle,anglax(3,3),coe, angle0,DSTRa(6),DSTRE(6),beta_pre(2)
      REAL*8  strsg_p,TER(2),tracc(2),MAXQI,tenspg(nelem,NSTR1),factor_i,beta_impli(2)
      real*8  NX,NY,KBETA,FACTOR,sigmacri,der_chi,normt,COEFI,qe_inv(2,2),
     .        qfi_0(2,2),dtracc(2), qe(2,2),eps_el(nstr1),nc(4,2),MULAM,LAMBL, angax(2),
     .        mbeta,nu(2),nu1(2),mod_nu,ang_hist(6)
	
	REAL*8  VTRAC(NDIME),mod_beta_ax,tn,ts,bn,bs,nxx,nyy,first_t,damage
	REAL*8  qfi(2,2),qfi_inv(2,2),cartd_b(4,2),n_m(4,2),beta_S(2),b(4,2),sig_tilde
      real*8  gamma_pre,mod_gamma_pre,gamma_pre_szz,mx,my,tzz,tmm,gamma_s,mod_gamma,
     .        mod_gamma_pre_zz,qbifu_g_zz,gamma_szz,mod_gamma_zz,var_g(6),sigma_mm,
     .        gfval_el,qbifu_g,COTK,hintr,EXPOH
      PARAMETER (PIMED=1.57079632679490D0)
C
C***Recover properties  
C
      YOUNG=PROPS(11)                      ! YOUNG MODULUS
      POISS=PROPS(12)                      ! POISSON RATIO
      MULAM=YOUNG/2.D0/(1.D0+POISS)
      LAMBL=YOUNG*POISS/(1.D0+POISS)/(1.D0-2.D0*POISS)
      FTULT=PROPS(13)                      ! ULTIMATE STRESS (sigma_ult)
      GFVAL=PROPS(14)                      ! FRACTURE ENERGY
	IF (KGEOM.EQ.1) THEN
         LAMBL= 2.D0*LAMBL*MULAM/(2.D0*MULAM+LAMBL)
      ENDIF
      HPOUT=PROPS(22)                      ! SOFTENING PARAMETER (H_y)
	EXPOH=PROPS(15) 
      KBETA= 1E-10/young

 	gfval_el= gfval/dsqrt(dvolt/THICK)


      implicito=impli
c
cmayo23      LTRAK=int(ESOLD(1))
cmayo23      LFLAG=LTFLG(IELEM) 
C
      MOD_BETA_PRE=MOD_BETA
      mod_beta_tilde=mod_beta 
      gamma_pre       = var_g(1) 
      mod_gamma_pre   = var_g(2)       
      qbifu_g         = var_g(3)
c
	gamma_s   = gamma_pre
      mod_gamma = mod_gamma_pre
c
      gamma_pre_szz      = var_g(4) 
      mod_gamma_pre_zz   = var_g(5)       
      qbifu_g_zz         = var_g(6) 
c
	gamma_szz          = gamma_pre_szz 
      mod_gamma_zz       = mod_gamma_pre_zz
C
C***  Evaluation of SIGMA_efective
C
      call veczer(4,qfi_inv)
      call veczer(8,cartd_b)
      call veczer(8,n_m)

      qfi_inv(1,1)=kbeta
      qfi_inv(2,2)=kbeta

      beta_pre(1)=beta_s(1)
      beta_pre(2)=beta_s(2)
c
c***  matriz N
c
      NX=DCOS(ESOLD(5))
      NY=DSIN(ESOLD(5))
      call nsolitario(nx,ny,CARTD,nsol_b,ndime)

      MX=-NY
      MY= NX

c
      if(simet) then
	  coe=dabs(nx*CARTD(1,nsol_b)+ny*CARTD(2,nsol_b))
	  if(coe.eq.0.d0) coe=1.d0
	  nxx=CARTD(1,nsol_b)/coe
	  nyy=CARTD(2,nsol_b)/coe
	  n_m(1,1)=nxx
	  n_m(3,1)=nyy
	  n_m(2,2)=nyy
	  n_m(3,2)=nxx
      else
	  n_m(1,1)=nx
	  n_m(3,1)=ny
	  n_m(2,2)=ny
	  n_m(3,2)=nx
	endif
c***                 matrix B= C*b+
      cartd_b(1,1)=CARTD(1,nsol_b)
      cartd_b(3,1)=CARTD(2,nsol_b)
      cartd_b(2,2)=CARTD(2,nsol_b)
      cartd_b(3,2)=CARTD(1,nsol_b)
C
C***  Evaluation of SIGMA_efective
C
      CALL ELSTR(stref,stran,DMATX)           ! STREF= dmatx*stran
c
      do istr=1,NSTR1
	 do jstr=1,2
   	   B (istr,jstr)=0.d0
   	   nc(istr,jstr)=0.d0
	   do kstr=1,NSTR1
	       B(istr,jstr) = B(istr,jstr)  + DMATX(istr, kstr)*cartd_b(kstr,jstr)
	       nc(istr,jstr)= nc(istr,jstr) + DMATX(istr, kstr)*n_m(kstr,jstr)
	   enddo
	 enddo
      enddo
	


	qe(1,1)= (mulam+lambl)*nx*nx +mulam
	qe(2,2)= (mulam+lambl)*ny*ny +mulam
	qe(1,2)= (mulam+lambl)*ny*nx  
	qe(2,1)= qe(1,2) 
c
c***   varibles de danio estandar
c
      RINT0 = FTULT/DSQRT(YOUNG)
      ZEROQ=1.D-8*RINT0
      if (istep.eq.1.and.itime.eq.1) dr_pre=0.0d0

      IF(QFUNR.LE.0.D0)QFUNR=RINT0
      IF(RINTV.LE.0.D0)RINTV=RINT0
      RIPRE=RINTV
      QFPRE=QFUNR
	IF (ITIME.EQ.1.AND.ISTEP.EQ.1.AND.IITER.EQ.0.
     .              AND.IGAUL.EQ.1.AND.IELEM.EQ.1) THEN
ccc  este valor hay que modificarlo
	  QINCR=dabs(RINT0/HPOUT)*QINCR
	ENDIF
	rtilde=ripre
	qtilde=qfunr
	IF (IITER.EQ.0) THEN
	  CALL ACTUALIZAR(dr_pre,%VAL(P_INCRA),NELEM,IELEM,1,igaul)
	ENDIF

      if (istep.gt.5.or.itime.gt.1) THEN
	     call getrtil(ripre,%VAL(P_INCRA),%VAL(P_INCRF),NELEM,IELEM,rtilde,igaul)
cc           IF(LFLAG.NE.20000) THEN 
                   qtilde= qfpre+ HPOUT* (rtilde-ripre)
cc	     else
cc                   qtilde= qfpre+ HPOUT* (-1)* (rtilde-ripre)
cc	     endif
	   if (qtilde.le.zeroq) qtilde=zeroq
         if (qtilde.gt.rint0) qtilde=rint0
	endif
c.......................... 
      CALL POSI2D(NSTR1,STREF,STPOS)
      RTRIAL=DSQRT(VECDOT(STPOS,STRAN,NSTR1))

      IF (IMPLICITO.EQ.1)  COEFI=QFUNR/RINTV
      IF (IMPLICITO.EQ.0)  COEFI=QTILDE/RTILDE    !totalmente explicito
      CALL VECSCA(NSTR1,COEFI,STREF,STRSG) 
c
c***                  vector  tracc= N*C*epsilon
c
      tracc(1)= n_M(1,1)*STRSG(1)+ n_M(3,1)*STRSG(3)
      tracc(2)= n_M(2,2)*STRSG(2)+ n_M(3,2)*STRSG(3)
	tmm= mx* (mx*STRSG(1)+ my*STRSG(3)) + 
     .          my* (my*STRSG(2)+ mx*STRSG(3) )
	tzz=strsg(4)
c.......................................................
      FLOAD=1.D0
      FACTM=0.D0
      FACTP=0.D0
      der_chi=0.d0
	factor=1.d0
C
c***
cmayo23      IF(LFLAG.NE.20000) THEN 
c
        IF(qbifu.eq.0.d0)THEN
           IF(RTRIAL.GT.RIPRE)THEN
C***         loading condition 
C            -----------------
             FLOAD=-1.D0 
             DELTR=RTRIAL-RIPRE
             RINTV=RTRIAL
CCCCCC evaluar el hpout (qfpre)
             HINTR=hpout
             COTK=-RINT0/GFVAL
             IF(EXPOH.EQ.1.0D0) HINTR=COTK*QFPRE 


cccccccccccccccccccccccc
             QFUNR=QFPRE+hintr*DELTR
             IF (QFUNR.LT.ZEROQ) QFUNR=ZEROQ
	       if (qtilde.le.zeroq) qtilde=zeroq
             if (qtilde.gt.rint0) qtilde=rint0  
             IF (IMPLICITO.EQ.1)  COEFI=QFUNR/RINTV
	       IF (IMPLICITO.EQ.0)  COEFI=QTILDE/RTILDE    !totalmente explicito
	       factor= coefi
             CALL CRIBIS(STREF,DUMMY,PROPS,QFUNR,RINTV,angax,HPCRI,ltrak,RTRIAL)

             IF (HPCRI.GE.HPOUT) THEN   
                 QBIFU= DSQRT(tracc(1)**2+tracc(2)**2)
                 LTRAK = -1 
	           first_t=1.d0
                 FLOAD=-1.D0 
	           qfunr=qfpre
	           rintv=ripre
                 IF (IMPLICITO.EQ.1)  COEFI=QFUNR/RINTV
                 IF (IMPLICITO.EQ.0)  COEFI=QTILDE/RTILDE    !totalmente explicito
             ELSE
                IF (IMPLICITO.EQ.1)  FACTM=HPOUT/RINTV**2 - QFUNR/RINTV**3
	          IF (IMPLICITO.EQ.0)  factm=0.0d0                      !totalmente explicito
             ENDIF
	     ENDIF
           fload=0.d0
           beta_s(1)= qfi_inv(1,1) * tracc(1)+ qfi_inv(1,2) * tracc(2)
           beta_s(2)= qfi_inv(2,1) * tracc(1)+ qfi_inv(2,2) * tracc(2)
	     mod_beta=dsqrt(qe(1,1)*beta_s(1)**2+
     .           (qe(1,2)+qe(2,1)) *beta_s(1)*beta_s(2)+qe(2,2) * beta_s(2)**2)
           beta_impli(1)=beta_s(1)
           beta_impli(2)=beta_s(2)		 
	     factor=COEFI
        ELSE
C***       bifurcation condition has been set TRUE IN PREVIOUS STEPS
C          ---------------------------------------------------------
           QFUNR=QFPRE
           IF(qtilde.LE.ZEROQ) qtilde=zeroq
c
           do istr=1,2
	      do jstr=1,2
   	        qfi_0(istr,jstr)=0.d0
	        do kstr=1,NSTR1
	  	      qfi_0(istr,jstr)=qfi_0(istr,jstr)+ n_m(kstr, istr)*b(kstr,jstr)
	        enddo
	  	    qfi_0(istr,jstr)=qfi_0(istr,jstr)*coefi
	      enddo
           enddo
c***       inicializacion beta
c          -------------------
	     if(first_t.eq.1.d0) then
	       first_t=0.d0
             mod_beta =  MOD_BETA_PRE
	       MOD_BETA_TILDE=mod_beta
             dete= 1/(qe(1,1)*qe(2,2) - qe(1,2)*qe(2,1))
             qe_inv(1,1) =  qe(2,2)* dete
             qe_inv(1,2) = -qe(1,2)* dete
             qe_inv(2,1) = -qe(2,1)* dete
             qe_inv(2,2) =  qe(1,1)* dete
             nu(1)=  qe_inv(1,1) * tracc(1)+ qe_inv(1,2) * tracc(2)
             nu(2)=  qe_inv(2,1) * tracc(1)+ qe_inv(2,2) * tracc(2)
             dete= dsqrt(nu(1)*nu(1)+ nu(2)*nu(2))
             nu(1)=nu(1)/dete
             nu(2)=nu(2)/dete

	       nu1(1)=  qe (1,1) * nu(1)+ qe(1,2) * nu(2)
	       nu1(2)=  qe (2,1) * nu(1)+ qe(2,2) * nu(2)
	       mod_nu=dsqrt(nu1(1)*nu(1) + nu1(2)*nu(2))

	       qbifu= (tracc(1)*nu(1)+tracc(2)*nu(2))/mod_nu
	        
		   mbeta= mod_beta/mod_nu
	       beta_s(1)=  mbeta* nu(1)
	       beta_s(2)=  mbeta* nu(2)
	       beta_impli(1)  = beta_s(1) 
             beta_impli(2)  = beta_s(2)  

             fload=-1
           else
c***
             factor= qbifu/MOD_BETA_TILDE*dexp(-qbifu*MOD_BETA_TILDE/GFVAL)
             qfi(1,1)= factor*qe(1,1)+qfi_0(1,1)
             qfi(2,2)= factor*qe(2,2)+qfi_0(2,2)
             qfi(1,2)= factor*qe(1,2)+qfi_0(1,2)
             qfi(2,1)= factor*qe(2,1)+qfi_0(2,1)
             dete= 1/(qfi(1,1)*qfi(2,2) - qfi(1,2)*qfi(2,1))
             qfi_inv(1,1) =  qfi(2,2)* dete
             qfi_inv(1,2) = -qfi(1,2)* dete
             qfi_inv(2,1) = -qfi(2,1)* dete
             qfi_inv(2,2) =  qfi(1,1)* dete
             beta_impli(1)=  qfi_inv(1,1) * tracc(1)+ qfi_inv(1,2) * tracc(2)
             beta_impli(2)=  qfi_inv(2,1) * tracc(1)+ qfi_inv(2,2) * tracc(2)
      
             beta_s(1)= beta_impli(1)  
             beta_s(2)= beta_impli(2)  
      
	       mod_beta=dsqrt(qe(1,1)*beta_impli(1)**2             +
     .                  (qe(1,2)+qe(2,1)) *beta_impli(1)*beta_impli(2) +
     .                   qe(2,2) * beta_impli(2)**2   )
		   fload=-2
           endif
C  
           if(MOD_BETA .gt. MOD_BETA_pre)   then   
c***          loading in strong discontinuity
c             -------------------------------
              fload=1
c***          NEWTON
C		    ------
              call newton_do(qbifu,mod_beta,gfval,qfi_0,qe,tracc, beta_impli,
     .                     factor_i,der_chi,lambl,mulam,ncrit,nx,ny)
c
	        if(mod_beta.lt. mod_beta_pre) then
	           beta_impli(1)=beta_s(1)
	           beta_impli(2)=beta_s(2)
	           mod_beta=mod_beta_pre
	           fload=-2
                 factor= qbifu/mod_beta*dexp(-qbifu*mod_beta/GFVAL)
	           fload= 1.d0
	           der_chi=0.d0
              else
c             ----------
                if (impli) then
                   beta_s(1)=beta_impli(1)
                   beta_s(2)=beta_impli(2)
                   factor=factor_i
	          endif
	        endif
		    ltrak = -1
	     else     ! 
c***          unloading in strong discontinuity
c             ---------------------------------
              MOD_BETA=MOD_BETA_pre		  
              if(fload.eq.-2) then
			   	
	        else
                 fload=-1
	           der_chi=0.d0
              endif

 	     endif
c**************************************************************
c          % MODO GAMMA
c
           eps_salto(1)= stran(1) - (cartd_b(1,1)* beta_pre(1)+ cartd_b(1,2)* beta_pre(2) )
           eps_salto(2)= stran(2) - (cartd_b(2,1)* beta_pre(1)+ cartd_b(2,2)* beta_pre(2) )
           eps_salto(3)= stran(3) - (cartd_b(3,1)* beta_pre(1)+ cartd_b(3,2)* beta_pre(2) )
           eps_salto(4)= stran(4) - (cartd_b(4,1)* beta_pre(1)+ cartd_b(4,2)* beta_pre(2) )
           CALL ELSTR(STREF,eps_salto,DMATX)           ! STREF= dmatx*stran
           CALL VECSCA(NSTR1,COEFI,STREF,STRSG) 
           sigma_mm  = mx*(strsg(1)*mx+strsg(3)*my)+my*(mx* strsg(3)+my* strsg(2))
c***       tension mm
c          ----------
           if (qbifu_g .eq. 0.d0 ) then

	        tmm = sigma_mm
	  
c	        if( dabs(sigma_mm).gt.ftult ) then  !MODELO SIMETRICO
	        if( (sigma_mm).gt.ftult ) then     !MODELO NO-SIMETRICO
			    qbifu_g =dabs(sigma_mm)
                  mod_gamma   = (dabs(sigma_mm) -ftult)/YOUNG
                  gamma_s=mod_gamma
ccc				tmm=dsign(1.d0,sigma_mm)*ftult 
	        endif

           else
ccccc	        gamma_s     =  (sigma_mm -ftult)/YOUNG

	        gamma_s     =  (sigma_mm -ftult)/YOUNG

c              mod_gamma   = (dabs(sigma_mm) -ftult)/YOUNG   !MODELO SIMETRICO
               mod_gamma   = ((sigma_mm) -ftult)/YOUNG    !MODELO NO-SIMETRICO


	        if(mod_gamma.lt.mod_gamma_pre) mod_gamma=mod_gamma_pre
              write(*,*)     ' paso por aqui, elem= ', ielem 
              write(LULOG,*) ' Modo gamma activo, elem= ', ielem
              if(impli) then
cc			  sig_tilde= Young* mod_gamma+ftult
			  sig_tilde= YOUNG* mod_gamma+ftult



                tmm = qbifu_g*dexp(-qbifu_g*mod_gamma/gfval_el)*sigma_mm/sig_tilde
              else
			  sig_tilde= Young* mod_gamma_pre+ftult
                tmm = qbifu_g*dexp(-qbifu_g*mod_gamma_pre/gfval_el)*sigma_mm/sig_tilde
	        endif


           endif
c***       tension zz
c          ----------
           if (qbifu_g_zz .eq. 0.d0.and. ntype.ne.1 ) then
	        tzz = strsg(4)
c 	        if( dabs(strsg(4)).gt.ftult ) then  !MODELO SIMETRICO
	        if( (strsg(4)).gt.ftult ) then     !MODELO NO-SIMETRICO
			    qbifu_g_zz =dabs(strsg(4))
                  mod_gamma_zz   = (dabs(strsg(4)) -ftult)/YOUNG
	            gamma_szz=mod_gamma_zz
	        endif
           else
	        gamma_szz     =  (strsg(4) -ftult)/YOUNG
c              mod_gamma_zz   = (dabs(strsg(4)) -ftult)/YOUNG   !MODELO SIMETRICO
              mod_gamma_zz   = ((strsg(4)) -ftult)/YOUNG    !MODELO NO-SIMETRICO
	        if(mod_gamma_zz.lt.mod_gamma_pre_zz) mod_gamma_zz=mod_gamma_pre_zz
              if(impli) then
			  sig_tilde= Young* mod_gamma_zz+ftult
                tzz = qbifu_g_zz*dexp(-qbifu_g_zz*mod_gamma_zz/gfval_el)*strsg(4)/sig_tilde
              else
			  sig_tilde= Young* mod_gamma_pre_zz+ftult
                tzz = qbifu_g_zz*dexp(-qbifu_g_zz*mod_gamma_pre_zz/gfval_el)*strsg(4)/sig_tilde
	        endif
           endif

        endif

cmayo23	else
cmayo23c***  elementos fuera del tracking
cmayo23c
cmayo23           fload=0.d0
cmayo23	     factor=1.d0
cmayo23           beta_impli(1)=beta_pre(1)
cmayo23           beta_impli(2)=beta_pre(2)
cmayo23           beta_s(1)=beta_pre(1)
cmayo23           beta_s(2)=beta_pre(2)
cmayo23           MOD_BETA=MOD_BETA_pre		 
cmayo23        if(0) then		  
cmayo23ccccc           IF(RTRIAL.GT.RIPRE)THEN
cmayo23C***         loading condition 
cmayo23C            -----------------
cmayo23             FLOAD=-1.D0 
cmayo23             DELTR=RTRIAL-RIPRE
cmayo23             RINTV=RTRIAL
cmayo23             QFUNR=QFPRE+1.d-1 * DELTR
cmayo23             IF (QFUNR.LT.ZEROQ) QFUNR=ZEROQ
cmayo23	       if (qtilde.le.zeroq) qtilde=zeroq
cmayo23             if (qtilde.gt.rint0) qtilde=rint0  
cmayo23             IF (IMPLICITO.EQ.1)  COEFI=QFUNR/RINTV
cmayo23	       IF (IMPLICITO.EQ.0)  COEFI=QTILDE/RTILDE    !totalmente explicito
cmayo23	       factor= coefi
cmayo23	     endif
cmayo23
cmayo23      ENDIF
c
C***Stress evaluation
c
 200  CONTINUE
      eps_salto(1)= stran(1) - (cartd_b(1,1)* beta_s(1)+ cartd_b(1,2)* beta_s(2) )
      eps_salto(2)= stran(2) - (cartd_b(2,1)* beta_s(1)+ cartd_b(2,2)* beta_s(2) )
      eps_salto(3)= stran(3) - (cartd_b(3,1)* beta_s(1)+ cartd_b(3,2)* beta_s(2) )
      eps_salto(4)= stran(4) - (cartd_b(4,1)* beta_s(1)+ cartd_b(4,2)* beta_s(2) )

	inc_eps_el(1)=eps_salto(1)- eps_el(1)
	inc_eps_el(2)=eps_salto(2)- eps_el(2)
	inc_eps_el(3)=eps_salto(3)- eps_el(3)
	inc_eps_el(4)=eps_salto(4)- eps_el(4)

	eps_el(1)=eps_salto(1)
	eps_el(2)=eps_salto(2)
	eps_el(3)=eps_salto(3)
	eps_el(4)=eps_salto(4)

      CALL ELSTR(strsg,eps_salto,DMATX)           ! STREF= dmatx*stran
      IF (IMPLICITO.EQ.1)  FACTQ=QFUNR/RINTV
	IF (IMPLICITO.EQ.0)  FACTQ=QTILDE/RTILDE         !totalmente explicito

	damage=1.d0-factq

      CALL VECSCA(NSTR1,FACTQ,strsg,STRSG)

      if(qbifu_g.gt.0.d0) CALL ten_max(strsg,PROPS,nx,ny,mx,my,tmm,tzz)
	 
 	dr_pre=rintv-ripre
	CALL ACTUALIZAR(dr_pre,%VAL(P_INCRA),NELEM,IELEM,2,igaul)

      beta_s(1)=beta_impli(1)
      beta_s(2)=beta_impli(2)

      dq_beta= MOD_BETA - MOD_BETA_pre
      if(dq_beta.lt. 0.d0) dq_beta=0.d0

      factp = factor
      factm = der_chi
      
      tracc(1)= n_M(1,1)*strsg(1)+ n_M(2,1)*strsg(2)+ n_M(3,1)*strsg(3)
      tracc(2)= n_M(1,2)*strsg(1)+ n_M(2,2)*strsg(2)+ n_M(3,2)*strsg(3)
      normt=dsqrt(tracc(1)**2 + tracc(2)**2)
c
      CALL CRIBIS(strsg,DUMMY,PROPS,QFUNR,RINTV,angax,HPCRI,ltrak,RTRIAL)
c
	ANGLN(1,1) = angax(1)
	ANGLN(2,1) = angax(2)
      if (kangl.eq.2)then
	     daux=angln(1,1) 
	     angln(1,1) = angln(2,1)
           angln(2,1) = daux
	endif


c	if(qbifu.eq.0.d0) then
c	   ANGLN(1,1) = angax(1)
c	   ANGLN(2,1) = angax(2)
c         if (kangl.eq.2)then
c	     daux=angln(1,1) 
c	     angln(1,1) = angln(2,1)
c           angln(2,1) = daux
c	   endif
c	endif

C

      ANGL2(IELEM,1)= fload
      ANGL2(IELEM,2)= damage 
	if(qbifu.gt.0.d0) then 
       ANGL2(IELEM,3)=  dabs((beta_s(1)*nx+beta_s(2)*ny)) /  ! modo I -> angl2=1
     .          (dsqrt(beta_s(1)**2+beta_s(2)**2))    ! modoII -> angl2=0
       ANGL2(IELEM,4)= gamma_s                
       ANGL2(IELEM,5)= dsqrt(beta_s(1)**2+beta_s(2)**2)                    
	else
       ANGL2(IELEM,3)= 1.d0 
       ANGL2(IELEM,4)= 0.d0 
       ANGL2(IELEM,5)= 0.d0                    
	endif

	do istr=1,nstr1
	  tenspg(ielem,istr)= strsg(istr)
	enddo   

	if(ntype.ne.4) then
        ANGL1(IELEM,1)=esold(5)+pimed
	else
        ANGL1(IELEM,1)=ANGLN(1,1)
        ANGL1(IELEM,2)=ANGLN(2,1)
        ANGL1(IELEM,3)=ANGLN(3,1)
      endif
c
C***Store the variables in the global arrays
C
      ESCUR(1)=FLOAT(LTRAK)
      ESCUR(4)=QBIFU
      ESCUR(3)=normt
c
      var_g(1)= gamma_s  
      var_g(2)= mod_gamma       
      var_g(3)= qbifu_g
c
      var_g(4)= gamma_szz  
      var_g(5)= mod_gamma_zz       
      var_g(6)= qbifu_g_zz
      RETURN
      END
C   
C************************************************************************
      SUBROUTINE newton_DO(qbifu,mod_beta,gfval,qfi_0,qe,tracc, beta_impli,
     .                     factor_i,der_chi,lambl,mulam,ncrit,nx,ny)  
C************************************************************************
C
      IMPLICIT NONE
C
      real*8  qbifu,mod_beta,gfval,qfi_0(2,2),tracc(2),beta_impli(2),
     .        factor_i,der_factor,der_chi,qfi(2,2),res_b(2),mod_res,
     .        normt,k(2,2), k_inv(2,2), dete, delta_b(2),qe(2,2),lambl,mulam,
     .        coe1,coe2,qeb(2),nx,ny,bn,bs,derbn,derbs,expone
	integer iter,ncrit
c     
c***  Newton      ****
c     ------
      factor_i = qbifu/mod_beta*dexp(-qbifu*mod_beta/GFVAL)
      der_factor= (-qbifu**2/GFVAL* dexp(-qbifu*mod_beta/GFVAL)) 
      der_chi = (der_factor  - factor_i ) / mod_beta**2
	coe1= (lambl+2.d0*mulam)
	coe2= mulam 
	 
	bn =   beta_impli(1) * nx + beta_impli(2)* ny
	bs = - beta_impli(1) * ny + beta_impli(2)* nx

      qfi(1,1)= factor_i*qe(1,1)+qfi_0(1,1)
      qfi(2,2)= factor_i*qe(2,2)+qfi_0(2,2)
      qfi(1,2)= factor_i*qe(1,2)+qfi_0(1,2)
      qfi(2,1)= factor_i*qe(2,1)+qfi_0(2,1)

      res_b(1)= qfi(1,1)*beta_impli(1) + qfi(1,2)*beta_impli(2) -tracc(1)
      res_b(2)= qfi(2,1)*beta_impli(1) + qfi(2,2)*beta_impli(2) -tracc(2) 
      mod_res= dsqrt(res_b(1)*res_b(1)+res_b(2)*res_b(2))   
      normt= dsqrt(tracc(1)*tracc(1)+tracc(2)*tracc(2))   
      iter=0

      do while((mod_res.gt.1.e-10*normt).and.iter.lt.100)
          iter=iter+1

          qeb(1)= qe(1,1)*beta_impli(1) +  qe(1,2)*beta_impli(2)
          qeb(2)= qe(2,1)*beta_impli(1) +  qe(2,2)*beta_impli(2)

          derbn= coe1* bn*nx +coe2*bs*(-ny)
	    derbs= coe1* bn*ny +coe2*bs*( nx)

          K(1,1)= qfi(1,1) + der_chi * qeb(1) *derbn 
          K(1,2)= qfi(1,2) + der_chi * qeb(1) *derbs
          K(2,1)= qfi(2,1) + der_chi * qeb(2) *derbn
          K(2,2)= qfi(2,2) + der_chi * qeb(2) *derbs

	    dete= 1.d0/(k(1,1)*k(2,2) - k(1,2)*k(2,1))
          k_inv(1,1) =  k(2,2)* dete
          k_inv(1,2) = -k(1,2)* dete
          k_inv(2,1) = -k(2,1)* dete
          k_inv(2,2) =  k(1,1)* dete
			
		delta_b(1)= -(k_inv(1,1)*res_b(1)+k_inv(1,2)*res_b(2) )
		delta_b(2)= -(k_inv(2,1)*res_b(1)+k_inv(2,2)*res_b(2) )
c***
          beta_impli(1) = beta_impli(1) + delta_b(1)
          beta_impli(2) = beta_impli(2) + delta_b(2)

 	    mod_beta=dsqrt(qe(1,1)*beta_impli(1)**2             +
     .                (qe(1,2)+qe(2,1)) *beta_impli(1)*beta_impli(2) +
     .                 qe(2,2) * beta_impli(2)**2   )

	    bn =   beta_impli(1) * nx + beta_impli(2)* ny
	    bs = - beta_impli(1) * ny + beta_impli(2)* nx
          expone    = dexp(-qbifu*MOD_BETA/GFVAL)
          factor_i  = (qbifu/MOD_BETA)  * expone
          der_factor= (-qbifu**2/GFVAL) * expone
		der_chi   = (der_factor  - factor_i  ) / MOD_BETA**2
          qfi(1,1)= factor_i * qe(1,1) + qfi_0(1,1)
          qfi(2,2)= factor_i * qe(2,2) + qfi_0(2,2)
          qfi(1,2)= factor_i * qe(1,2) + qfi_0(1,2)
          qfi(2,1)= factor_i * qe(2,1) + qfi_0(2,1)

	    res_b(1)= qfi(1,1)*beta_impli(1) + qfi(1,2)*beta_impli(2) -tracc(1) 
	    res_b(2)= qfi(2,1)*beta_impli(1) + qfi(2,2)*beta_impli(2) -tracc(2) 
          mod_res= dsqrt(res_b(1)*res_b(1)+res_b(2)*res_b(2))
 
      enddo
	if(iter.ge.99) stop ' modulo constitutivo no converge'
c
c***  end Newton    *****
      return
	end
c   
C************************************************************************
      SUBROUTINE MTXDVB_do(PROPS,STRAN,DMATX,STRSG, qbifu,beta_s,
     .      STREF,FLOAD,FACTM,FACTP,FACTQ,DEFIC,ESCUR,esold,CARTD )
C************************************************************************
C
C    THIS ROUTINE EVALUATES MATRIX FOR DAMAGE MODELS
C                   (STRONG DISCONTINUITIES)
C
C************************************************************************
C
!    INPUT VARIABLES :
C
C       PROPS = Properties     
C       STRAN = Strain tensor     
C       STRSG = Stress tensor
C       STREF = Efective Stress  
C       RINTV = Internal varible describing the size of the 
C               elastic region   
C       DAMAG = Scalar damage variable  (q)
C       FLOAD = Load/Unload flag   
C       FACTM = A factor to be used in the tangent matrix evaluation
C       HPRIM = Softening parameter
C       DEFIC = Incompatible deformation
C
!    OUTPUT VARIABLES :
C
C       DMATX = Constitutive matrix
C
C************************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'element.om'  ! -->  NHIST, NSTR1
      INCLUDE 'propert.om'  ! -->  PH,PA
      INCLUDE 'interval.om' ! -->  DTIME
      INCLUDE 'auxiliar.om' ! -->  IGAUL
      INCLUDE 'pararcl.om' ! -->  impli
      INCLUDE 'cstrong.om'                                   ! STDC
C
      REAL*8  PROPS(*), STRSG(NSTR1), STREF(NSTR1),STPOS(6)
      REAL*8  DMATX(NSTR1,NSTR1), STRAN(NSTR1),DEFIC(NSTR1)
C      
      REAL*8  RINTV,DAMAG,FLOAD,FACTM,FACTP,FACTQ,HPRIM,COEF1,PSIGMA,
     .        THETA,VECDOT,STRAX(9),ESCUR(*),esold(*),MULAM,lambl
	real*8 CARTD(2,*),KBETA,qfi_inv(2,2),qfi(2,2),cartd_b(4,2),
     .       nx,ny,n_m(4,2),b(4,2),dete,be(4,4),ben(4,4),cen(4,4),
     .       qbifu,beta_s(2),coe,nc(4,2),qe(2,2),nxx,nyy,qeb(2),
     .       derbn,derbs,coe1,coe2,bn,bs,POISS_TP,YOUNG_TP
 
      INTEGER KLOAD,ISTR,JSTR,KSTR,LSTRE,nsol_b,implicito
C
C***Recover elastic properties 
C 
      YOUNG=PROPS(11)
      POISS=PROPS(12)
      MULAM=YOUNG/2.D0/(1.D0+POISS)
      LAMBL=YOUNG*POISS/(1.D0+POISS)/(1.D0-2.D0*POISS)
	KBETA= 1E-8/young
      QBIFU=ESCUR(4)
	IF (ntype.EQ.1) THEN
		POISS_TP = POISS/(POISS+1)
	    YOUNG_TP = (1-(POISS/(POISS+1))**2)*YOUNG
          MULAM=YOUNG_TP/2.D0/(1.D0+POISS_TP)
          LAMBL=YOUNG_TP*POISS_TP/(1.D0+POISS_TP)/(1.D0-2.D0*POISS_TP)
	ENDIF
C
      implicito=impli
C
C***Compute elastic matrix 
C
      CALL MTELAS(DMATX) 
c
c***
c
      DMATX(1,1)=FACTQ*DMATX(1,1)
      DMATX(1,2)=FACTQ*DMATX(1,2)
      DMATX(1,4)=FACTQ*DMATX(1,4)
      DMATX(2,1)=FACTQ*DMATX(2,1)
      DMATX(2,2)=FACTQ*DMATX(2,2)
      DMATX(2,4)=FACTQ*DMATX(2,4)
      DMATX(4,1)=FACTQ*DMATX(4,1)
      DMATX(4,2)=FACTQ*DMATX(4,2)
      DMATX(4,4)=FACTQ*DMATX(4,4)
      DMATX(3,3)=FACTQ*DMATX(3,3)

      if(fload.eq.0.d0) return


	if (fload.eq.-1) then
         CALL POSI2D(NSTR1,STREF,STPOS)
         DMATX(1,1)=DMATX(1,1) + FACTM*STREF(1)*STPOS(1)
         DMATX(1,2)=DMATX(1,2) + FACTM*STREF(1)*STPOS(2)
         DMATX(1,3)=DMATX(1,3) + FACTM*STREF(1)*STPOS(3)
         DMATX(1,4)=DMATX(1,4) + FACTM*STREF(1)*STPOS(4)
         DMATX(2,1)=DMATX(2,1) + FACTM*STREF(2)*STPOS(1)
         DMATX(2,2)=DMATX(2,2) + FACTM*STREF(2)*STPOS(2)
         DMATX(2,3)=DMATX(2,3) + FACTM*STREF(2)*STPOS(3)
         DMATX(2,4)=DMATX(2,4) + FACTM*STREF(2)*STPOS(4)
         DMATX(3,1)=DMATX(3,1) + FACTM*STREF(3)*STPOS(1)
         DMATX(3,2)=DMATX(3,2) + FACTM*STREF(3)*STPOS(2)
         DMATX(3,3)=DMATX(3,3) + FACTM*STREF(3)*STPOS(3)
         DMATX(3,4)=DMATX(3,4) + FACTM*STREF(3)*STPOS(4)
         DMATX(4,1)=DMATX(4,1) + FACTM*STREF(4)*STPOS(1)
         DMATX(4,2)=DMATX(4,2) + FACTM*STREF(4)*STPOS(2)
         DMATX(4,3)=DMATX(4,3) + FACTM*STREF(4)*STPOS(3)
         DMATX(4,4)=DMATX(4,4) + FACTM*STREF(4)*STPOS(4)
	   return
	endif


      call veczer(4,qfi_inv)

	qfi_inv(1,1)=kbeta
	qfi_inv(2,2)=kbeta

      call veczer(8,cartd_b)
	call veczer(8,n_m)


      NX=DCOS(ESOLD(5))
      NY=DSIN(ESOLD(5))

      call nsolitario(nx,ny,CARTD,nsol_b,2)   !  ndime=2

c      nsol1=INT(ESOLD(8))
c	IF(NSOL1.EQ.0) NSOL1=1
c      NX=DCOS(ESOLD(5))
c      NY=DSIN(ESOLD(5))
	cartd_b(1,1)=CARTD(1,nsol_b)
	cartd_b(3,1)=CARTD(2,nsol_b)
	cartd_b(2,2)=CARTD(2,nsol_b)
	cartd_b(3,2)=CARTD(1,nsol_b)
      if(simet) then
	  coe=dabs(nx*CARTD(1,nsol_b)+ny*CARTD(2,nsol_b))
	  if(coe.eq.0.d0) coe=1.d0
	  nxx=CARTD(1,nsol_b)/coe
	  nyy=CARTD(2,nsol_b)/coe
	  n_m(1,1)=nxx
	  n_m(3,1)=nyy
	  n_m(2,2)=nyy
	  n_m(3,2)=nxx
      else
	  n_m(1,1)=nx
	  n_m(3,1)=ny
	  n_m(2,2)=ny
	  n_m(3,2)=nx
	endif
c
      do istr=1,NSTR1
	do jstr=1,2
   	  B(istr,jstr)=0.d0
   	  nc(istr,jstr)=0.d0
	  do kstr=1,NSTR1
	       B(istr,jstr)=B(istr,jstr)+ DMATX(istr, kstr)*cartd_b(kstr,jstr)
	       nc(istr,jstr)=nc(istr,jstr)+ DMATX(istr, kstr)*n_m(kstr,jstr)
	  enddo
	enddo
	enddo
C
      if (qbifu .ne.0  )then
            do istr=1,2
	      do jstr=1,2
   		    QFI(istr,jstr)=0.d0
	        do kstr=1,NSTR1
		       qfi(istr,jstr)=qfi(istr,jstr)+ n_m(kstr, istr)*b(kstr,jstr)
	        enddo
	      enddo
            enddo
	      qe(1,1)= (mulam+lambl)*nx*nx +mulam
	      qe(2,2)= (mulam+lambl)*ny*ny +mulam
	      qe(1,2)= (mulam+lambl)*ny*nx  
	      qe(2,1)= qe(1,2) 

            qfi(1,1)= factp*qe(1,1)+QFI(1,1)
            qfi(2,2)= factp*qe(2,2)+QFI(2,2)
            qfi(1,2)= factp*qe(1,2)+QFI(1,2)
            qfi(2,1)= factp*qe(2,1)+QFI(2,1)
c***************************
	      coe1= (lambl+2.d0*mulam)
	      coe2= mulam 
	      bn =   beta_s(1) * nx + beta_s(2)* ny
	      bs = - beta_s(1) * ny + beta_s(2)* nx
            qeb(1)= qe(1,1)*beta_s(1) +  qe(1,2)*beta_s(2)
            qeb(2)= qe(2,1)*beta_s(1) +  qe(2,2)*beta_s(2)

            derbn= coe1* bn*nx +coe2*bs*(-ny)
	      derbs= coe1* bn*ny +coe2*bs*( nx)
c***************************

        	  if (impli. and.fload.eq.1) then
              qfi(1,1)= qfi(1,1)+ factm* qeb(1) *derbn 
              qfi(1,2)= qfi(1,2)+ factm* qeb(1) *derbs
              qfi(2,1)= qfi(2,1)+ factm* qeb(2) *derbn
              qfi(2,2)= qfi(2,2)+ factm* qeb(2) *derbs
            endif

	      dete= 1/(qfi(1,1)*qfi(2,2) - qfi(1,2)*qfi(2,1))
            qfi_inv(1,1) =  qfi(2,2)* dete
            qfi_inv(1,2) = -qfi(1,2)* dete
            qfi_inv(2,1) = -qfi(2,1)* dete
            qfi_inv(2,2) =  qfi(1,1)* dete
      endif

c....................................................
c***  CE= B* qfi_inv * N * ce;
c....................................................
       do istr=1,4
	 do jstr=1,2
	   be(istr,jstr)=0.d0
	   do kstr=1,2
            be(istr,jstr)= be(istr,jstr)+ B(istr,kstr)*qfi_inv(kstr,jstr) 
	   enddo
	 enddo
	 enddo

       do istr=1,4
	 do jstr=1,4
	   ben(istr,jstr)=0.d0
	   do kstr=1,2
            ben(istr,jstr)= ben(istr,jstr)+ be(istr,kstr)*n_m(jstr,kstr) 
	   enddo
	 enddo
	 enddo

       do istr=1,4
	 do jstr=1,4
	   cen(istr,jstr)=0.d0
	   do kstr=1,4
            cen(istr,jstr)= cen(istr,jstr)+ ben(istr,kstr)*dmatx(kstr,jstr) 
	   enddo
	 enddo
	 enddo
 
       do istr=1,4
	 do jstr=1,4
	   dmatx(istr,jstr)=dmatx(istr,jstr)-cen(istr,jstr)
	 enddo
	 enddo
 
       RETURN
C
      END
C
cc*********************************************************************
	SUBROUTINE getrtil(ripre,INCRA,INCRF,NELEM,IELEM,rtilde,igaul)
cc*********************************************************************
	INTEGER  NELEM,IELEM,igaul
	REAL*8   ripre,INCRA(NELEM,5,*),INCRF(*),rtilde
	REAL*8   FACTOR, tauep, inc_tauep

      if(INCRF(2) .eq.0.d0) then
	    factor=1.d0
	else
      	FACTOR=INCRF(1)/INCRF(2) 
	endif
	rtilde=ripre+INCRA(IELEM,2,igaul)*factor 
	END
