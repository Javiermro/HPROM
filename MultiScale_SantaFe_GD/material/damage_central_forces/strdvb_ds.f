      SUBROUTINE STRDVB_DS(PROPS,STRAN,STRSG,DMATX,
     .                  STREF,giro_normal,QFUNR,MOD_BETA,QBIFU,dq_beta,
     .                  BETA_BIFU,beta_s, FLOAD, FACTM, FACTP,
     .                  FACTQ,ESCUR,ESOLD,ANGLN,HPCRI,STOLD,
     .                  dummy,
     .                  angl1,angl2,ltflg, dummy2 ,dummy3,tenspg,CARTD,first_t,
     .                  eps_el, ang_hist, var_g )
C************************************************************************
C
C     THIS ROUTINE EVALUATES STRESSES FOR DAMAGE  MODELS
C      (STRONG DISCONTINUITIES with variable Bandwidth)  
C    modelo cohesivo discreto de fuerzas centrales (Sancho y planas)
C       NCRIT = 79    Cohesive Model of Sancho-Planas
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
      INTEGER LTRAK, implicito,istr,jstr,kstr
      INTEGER LTFLG(NELEM),I,J,nsol_b
      REAL*8  PROPS(NPROP), STRAN(NSTR1), STRSG(NSTR1),MULAM,LAMBL,YOUNG_TP, POISS_TP
      REAL*8  STREF(NSTR1), DAMAG, ESCUR(*), ESOLD(*),ANGLN(NDIME,NDIME)
      REAL*8  DMATX(NSTR1,NSTR1), QFUNR, STOLD(NSTR1),CARTD(NDIME,NNODL),VECDOT
C
      REAL*8  FTULT,GFVAL,HPCRI, dete,eps_salto(4),FLOAD,QBIFU,MOD_BETA_PRE, 
     .        MOD_BETA_TILDE, MOD_BETA, FACTM,FACTP,FACTQ,dummy,dummy1,dummy2,
     .        dummy3,pimed,dq_beta,
     .        var_g(6),gamma_s, mod_gamma, qbifu_g, ANGL1(NELEM,*),ANGL2(NELEM,3)
      REAL*8  anglax(3,3),beta_pre(2),tracc(2),tenspg(nelem,NSTR1),
     .        factor_i,beta_impli(2), NX,NY,KBETA,FACTOR,sigmacri,der_chi,normt,
     .        COEFI, qfi_0(2,2),  angax(2),
     .        mbeta,nu(2),nu1(2),mod_nu,beta_bifu,eps_el(nstr1),ang_hist(6),tmm
	
	REAL*8  first_t,damage, qfi(2,2),qfi_inv(2,2),cartd_b(4,2),n_m(4,2),beta_S(2),
     .        b(4,2),mx,my, gamma_pre , mod_gamma_pre , mod_gamma_tilde, m_m(4,2),
     .        sigma_mm,gamma_impli,gfval_el,sig_tilde,gamma_pre_szz,tzz,
     .        mod_gamma_pre_zz, qbifu_g_zz, gamma_szz, mod_gamma_zz,t1_t2,giro_normal

      PARAMETER (PIMED=1.57079632679490D0)
C
C***Recover properties  
C
      YOUNG=PROPS(11)                      ! YOUNG MODULUS
      POISS=PROPS(12)                      ! POISSON RATIO
      MULAM=YOUNG/2.D0/(1.D0+POISS)
      LAMBL=YOUNG*POISS/(1.D0+POISS)/(1.D0-2.D0*POISS)
      FTULT=PROPS(13)                       ! ULTIMATE STRESS (sigma_ult)
      GFVAL=PROPS(14)                      ! FRACTURE ENERGY
	IF (KGEOM.EQ.1) THEN
		POISS_TP = POISS/(POISS+1)
	    YOUNG_TP = (1-(POISS/(POISS+1))**2)*YOUNG
          MULAM=YOUNG_TP/2.D0/(1.D0+POISS_TP)
          LAMBL=YOUNG_TP*POISS_TP/(1.D0+POISS_TP)/(1.D0-2.D0*POISS_TP)
	ENDIF
      
 	gfval_el= gfval/dsqrt(dvolt/THICK)
 
      KBETA= 1E-12/young

      implicito=impli
c
      MOD_BETA_PRE=MOD_BETA
      mod_beta_tilde=mod_beta 
      gamma_pre        = var_g(1) 
      mod_gamma_pre    = var_g(2)       
      qbifu_g          = var_g(3)
      mod_gamma_tilde  = mod_gamma_pre
	gamma_s          = gamma_pre
      mod_gamma        = mod_gamma_pre
c
      gamma_pre_szz    = var_g(4) 
      mod_gamma_pre_zz = var_g(5)       
      qbifu_g_zz       = var_g(6) 
	gamma_szz        = gamma_pre_szz 
      mod_gamma_zz     = mod_gamma_pre_zz
C
C***  Evaluation of SIGMA_efective
C
      call veczer(4,qfi_inv)
      call veczer(8,cartd_b)
      call veczer(8,n_m)
      call veczer(8,m_m)

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
	n_m(1,1)=nx
	n_m(3,1)=ny
	n_m(2,2)=ny
	n_m(3,2)=nx
c
      cartd_b(1,1)=CARTD(1,nsol_b)
      cartd_b(3,1)=CARTD(2,nsol_b)
      cartd_b(2,2)=CARTD(2,nsol_b)
      cartd_b(3,2)=CARTD(1,nsol_b)
C
C***  Evaluation of SIGMA_efective
C
      CALL ELSTR(stref,stran,DMATX)           ! STREF= dmatx*stran
c.......................... 
      coefi=1.d0
      CALL VECSCA(NSTR1,COEFI,STREF,STRSG) 
c***                  vector  tracc= N*C*epsilon
      tracc(1)= n_M(1,1)*STRSG(1)+ n_M(3,1)*STRSG(3)
      tracc(2)= n_M(2,2)*STRSG(2)+ n_M(3,2)*STRSG(3)

	tmm= mx*(mx*STRSG(1)+ my*STRSG(3)) + my*(my*STRSG(2)+ mx*STRSG(3))
	tzz=strsg(4)

      do istr=1,NSTR1
	do jstr=1,2
   	   B (istr,jstr)=0.d0
	   do kstr=1,NSTR1
	       B(istr,jstr) = B(istr,jstr) + DMATX(istr,kstr)*cartd_b(kstr,jstr)
	   enddo
	enddo
      enddo
c.......................................................
      FLOAD=0.D0
      FACTM=0.D0
      FACTP=0.D0
      der_chi=0.d0
	factor=1.d0
C
      IF(qbifu.eq.0.d0)THEN
c	
             beta_s(1)= qfi_inv(1,1) * tracc(1)+ qfi_inv(1,2) * tracc(2)
             beta_s(2)= qfi_inv(2,1) * tracc(1)+ qfi_inv(2,2) * tracc(2)
             MOD_BETA=dsqrt(beta_s(1)**2+beta_s(2)**2) 
             beta_impli(1)=beta_s(1)
             beta_impli(2)=beta_s(2)
		   fload=0.d0
		   CALL CRIBIS_sp(STRSG,PROPS,angax,sigmacri, t1_t2 )
cc		   CALL CRIBIS_sp1(STRSG,PROPS,angax,sigmacri)
             IF (sigmacri.gt. FTULT) THEN   
                 QBIFU= DSQRT(tracc(1)**2+tracc(2)**2)
ccc                 QBIFU= sigmacri
	           first_t=1.d0
	           beta_bifu= dsqrt(beta_s(1)*beta_s(1)+beta_s(2)*beta_s(2))
             ENDIF	         

        ELSE
C***      bifurcation condition has been set TRUE IN PREVIOUS STEPS
C         ---------------------------------------------------------
          do istr=1,2
	      do jstr=1,2
   	        qfi_0(istr,jstr)=0.d0
	        do kstr=1,NSTR1
	  	     qfi_0(istr,jstr) =qfi_0 (istr,jstr)+ n_m(kstr, istr)*b(kstr,jstr)
	        enddo
	      enddo  	     
          enddo
c
c***      inicializacion beta
c         -------------------
	    if(first_t.eq.1.d0) then
cc	       first_t=0.d0
	       first_t=2.d0
             mod_beta = MOD_BETA_PRE
	       MOD_BETA_TILDE=mod_beta
             dete= dsqrt(tracc(1)*tracc(1)+ tracc(2)*tracc(2))
             nu(1)=tracc(1)/dete
             nu(2)=tracc(2)/dete
	       beta_s(1)=  mod_beta* nu(1)
	       beta_s(2)=  mod_beta* nu(2)
	       beta_impli(1)  = beta_s(1) 
             beta_impli(2)  = beta_s(2)  
             fload=-1
          else
c***
		   factor= qbifu/MOD_BETA_TILDE*dexp(-qbifu*MOD_BETA_TILDE/GFVAL)
             qfi(1,1)= factor + qfi_0(1,1)
             qfi(2,2)= factor + qfi_0(2,2)
             qfi(1,2)=          qfi_0(1,2)
             qfi(2,1)=          qfi_0(2,1)
             dete= 1/(qfi(1,1)*qfi(2,2) - qfi(1,2)*qfi(2,1))
             qfi_inv(1,1) =  qfi(2,2)* dete
             qfi_inv(1,2) = -qfi(1,2)* dete
             qfi_inv(2,1) = -qfi(2,1)* dete
             qfi_inv(2,2) =  qfi(1,1)* dete
             beta_impli(1)=  qfi_inv(1,1) * tracc(1)+ qfi_inv(1,2) * tracc(2)
             beta_impli(2)=  qfi_inv(2,1) * tracc(1)+ qfi_inv(2,2) * tracc(2)
             beta_s(1)= beta_impli(1)  
             beta_s(2)= beta_impli(2)  
             MOD_BETA=dsqrt(beta_impli(1)**2+beta_impli(2)**2) 
cnewton		   fload=-2
		   fload=1
          endif
C  


cc	   if(1) then     !! modelo con beta del paso previo 
cc
cc             if ((MOD_BETA .lt. MOD_BETA_pre))   then   
cc		         mod_beta=mod_beta_pre
cc                   fload=-2
cc	       endif
cc		   ltrak = -1
cc
cc	   else


            if ((MOD_BETA .gt. MOD_BETA_pre))   then   
c
c***          loading in strong discontinuity
c             -------------------------------
              fload=1
c  
c***          NEWTON
C		    ------
              call newton_ds(qbifu,mod_beta,gfval,qfi_0,tracc, beta_impli,
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

cc        endif


cc        if(1) then   !!!!!!!!!!!!         sin gamma   0    !!!! con gamma = 1
c
c          % MODO GAMMA
c
         if(giro_normal.eq.0.d0) then
c....................................................
           eps_salto(1)= stran(1) - (cartd_b(1,1)*beta_pre(1) + cartd_b(1,2)*beta_pre(2))
           eps_salto(2)= stran(2) - (cartd_b(2,1)*beta_pre(1) + cartd_b(2,2)*beta_pre(2))
           eps_salto(3)= stran(3) - (cartd_b(3,1)*beta_pre(1) + cartd_b(3,2)*beta_pre(2))
           eps_salto(4)= stran(4) - (cartd_b(4,1)*beta_pre(1) + cartd_b(4,2)*beta_pre(2))

           CALL ELSTR(strsg,eps_salto,DMATX)           ! STREF= dmatx*stran

c***  tension mm
           sigma_mm  = mx*(strsg(1)*mx+strsg(3)*my)+my*(mx* strsg(3)+my* strsg(2))
	     tmm       = sigma_mm   
           if (qbifu_g .eq. 0.d0 ) then
	        if( dabs(sigma_mm).gt.ftult ) then  !MODELO SIMETRICO
cc	        if( (sigma_mm).gt.ftult ) then     !MODELO NO-SIMETRICO
			    qbifu_g     =  dabs(sigma_mm)
                  mod_gamma   = (dabs(sigma_mm) -ftult)/YOUNG
                  gamma_s     =  mod_gamma
ccc				tmm=dsign(1.d0,sigma_mm)*ftult 
	        endif
           else
	        gamma_s   = (sigma_mm -ftult)/YOUNG
               mod_gamma = (dabs(sigma_mm) -ftult)/YOUNG  !MODELO SIMETRICO
cc               mod_gamma   = ((sigma_mm) -ftult)/YOUNG  !MODELO NO-SIMETRICO
	        if(mod_gamma.lt.mod_gamma_pre) mod_gamma=mod_gamma_pre
              write(*,*)     ' paso por aqui, elem= ', ielem 
              write(LULOG,*) ' Modo gamma activo, elem= ', ielem
              if(impli) then
			  sig_tilde= YOUNG*mod_gamma + ftult
                tmm = qbifu_g*dexp(-qbifu_g*mod_gamma/gfval_el)*sigma_mm/sig_tilde
              else
			  sig_tilde= YOUNG*mod_gamma_pre + ftult
                tmm = qbifu_g*dexp(-qbifu_g*mod_gamma_pre/gfval_el)*sigma_mm/sig_tilde
	        endif
           endif
c*** tension zz
           if (qbifu_g_zz .eq. 0.d0.and. KGEOM.ne.1  ) then
	        tzz = strsg(4)
	        if( dabs(strsg(4)).gt.ftult ) then  !MODELO SIMETRICO
cc	        if( (strsg(4)).gt.ftult ) then     !MODELO NO-SIMETRICO
			    qbifu_g_zz   =  dabs(strsg(4))
                  mod_gamma_zz = (dabs(strsg(4)) -ftult)/YOUNG
	            gamma_szz    =  mod_gamma_zz
	        endif
           else
	        gamma_szz    =(strsg(4) -ftult)/YOUNG
              mod_gamma_zz =(dabs(strsg(4)) -ftult)/YOUNG  !MODELO SIMETRICO
cc              mod_gamma_zz=((strsg(4)) -ftult)/YOUNG   !MODELO NO-SIMETRICO
	        if(mod_gamma_zz.lt.mod_gamma_pre_zz) mod_gamma_zz=mod_gamma_pre_zz
              if(impli) then
			  sig_tilde= YOUNG*mod_gamma_zz + ftult
                tzz = qbifu_g_zz*dexp(-qbifu_g_zz*mod_gamma_zz/gfval_el)*strsg(4)/sig_tilde
              else
			  sig_tilde= YOUNG*mod_gamma_pre_zz + ftult
                tzz = qbifu_g_zz*dexp(-qbifu_g_zz*mod_gamma_pre_zz/gfval_el)*strsg(4)/sig_tilde
	        endif
           endif
c....................................................
        endif

ccc	 endif    !!!!!! sin gamma

      endif
C***Stress evaluation
c
 200  CONTINUE
      eps_salto(1)= stran(1) - (cartd_b(1,1)* beta_s(1)+ cartd_b(1,2)* beta_s(2) )
      eps_salto(2)= stran(2) - (cartd_b(2,1)* beta_s(1)+ cartd_b(2,2)* beta_s(2) )
      eps_salto(3)= stran(3) - (cartd_b(3,1)* beta_s(1)+ cartd_b(3,2)* beta_s(2) )
      eps_salto(4)= stran(4) - (cartd_b(4,1)* beta_s(1)+ cartd_b(4,2)* beta_s(2) )
c
	inc_eps_el(1)=eps_salto(1)- eps_el(1)
	inc_eps_el(2)=eps_salto(2)- eps_el(2)
	inc_eps_el(3)=eps_salto(3)- eps_el(3)
	inc_eps_el(4)=eps_salto(4)- eps_el(4)

	eps_el(1)=eps_salto(1)
	eps_el(2)=eps_salto(2)
	eps_el(3)=eps_salto(3)
	eps_el(4)=eps_salto(4)

      CALL ELSTR(strsg,eps_salto,DMATX)           ! STREF= dmatx*stran

      factq=1.d0
	damage=1.d0-factq

      CALL VECSCA(NSTR1,FACTQ,strsg,STRSG)

      if(qbifu_g.gt.0.d0) CALL ten_max(strsg,PROPS,nx,ny,mx,my,tmm,tzz)

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
      CALL CRIBIS_sp(strsg,PROPS,angax,sigmacri,t1_t2)
cc      CALL CRIBIS_sp1(strsg,PROPS,angax,sigmacri)
cc      if(first_t.eq.2.and.t1_t2.gt.1.0d0)  first_t=0
      if(first_t.eq.2)  first_t=0

	ANGLN(1,1) = angax(1)
	ANGLN(2,1) = angax(2)

	ang_hist(6)= ang_hist(5)
	ang_hist(5)= ang_hist(4)
	ang_hist(4)= ang_hist(3)
	ang_hist(3)= ang_hist(2)
	ang_hist(2)= ang_hist(1)
	ang_hist(1)= angax(1)
C
      ANGL2(IELEM,1)= fload
      ANGL2(IELEM,2)= damage  
      ANGL2(IELEM,3)=  dabs((beta_s(1)*nx+beta_s(2)*ny)) /  ! modo I -> angl2=1
     .          (dsqrt(beta_s(1)**2+beta_s(2)**2))    ! modoII -> angl2=0
      ANGL2(IELEM,4)=  gamma_s                 
      ANGL2(IELEM,5)= dsqrt(beta_s(1)**2+beta_s(2)**2)                    

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
      ESCUR(2)=BETA_BIFU
      ESCUR(3)=0.d0
c
      var_g(1)= gamma_s  
      var_g(2)= mod_gamma       
      var_g(3)= qbifu_g
c
      var_g(4)= gamma_szz  
      var_g(5)= mod_gamma_zz       
      var_g(6)= qbifu_g_zz
c
      RETURN
      END
C   
C************************************************************************
      SUBROUTINE newton_ds(qbifu,mod_beta,gfval,qfi_0,tracc, beta_impli,
     .                     factor_i,der_chi,lambl,mulam,ncrit,nx,ny)  
C************************************************************************
C
      IMPLICIT NONE
C
      real*8  qbifu,mod_beta,gfval,qfi_0(2,2),tracc(2),beta_impli(2),
     .        factor_i,der_factor,der_chi,qfi(2,2),res_b(2),mod_res,
     .        normt,k(2,2), k_inv(2,2), dete, delta_b(2),lambl,mulam,
     .        coe1,coe2,qeb(2),nx,ny,bn,bs,derbn,derbs,expone,
     .        mod_res_0,beta_impli_0(2), s, mod_beta1 
	integer iter,ncrit,iterbs
c     
c***  Newton      ****
c     ------
      factor_i = qbifu/mod_beta*dexp(-qbifu*mod_beta/GFVAL)
      der_factor= (-qbifu**2/GFVAL* dexp(-qbifu*mod_beta/GFVAL)) 
      der_chi = (der_factor  - factor_i ) / mod_beta**2

      qfi(1,1)= factor_i + qfi_0(1,1)
      qfi(2,2)= factor_i + qfi_0(2,2)
      qfi(1,2)=          + qfi_0(1,2)
      qfi(2,1)=          + qfi_0(2,1)

      res_b(1)= qfi(1,1)*beta_impli(1) + qfi(1,2)*beta_impli(2) -tracc(1)
      res_b(2)= qfi(2,1)*beta_impli(1) + qfi(2,2)*beta_impli(2) -tracc(2) 
      mod_res= dsqrt(res_b(1)*res_b(1)+res_b(2)*res_b(2))   
      normt= dsqrt(tracc(1)*tracc(1)+tracc(2)*tracc(2))   
      iter=0


      do while((mod_res.gt.1.e-10*normt).and.iter.lt.100)
          iter=iter+1

          qeb(1)=  beta_impli(1) 
          qeb(2)=  beta_impli(2)

C         qeb(1)= beta_impli(1)
C         qeb(2)= beta_impli(2)
          derbn= beta_impli(1)
	    derbs= beta_impli(2)

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

          beta_impli(1) = beta_impli(1) + delta_b(1)
          beta_impli(2) = beta_impli(2) + delta_b(2)

          MOD_BETA=dsqrt(beta_impli(1)**2+beta_impli(2)**2) 

          expone    = dexp(-qbifu*MOD_BETA/GFVAL)
          factor_i  = (qbifu/MOD_BETA)  * expone
          der_factor= (-qbifu**2/GFVAL) * expone
		der_chi   = (der_factor  - factor_i  ) / MOD_BETA**2
          qfi(1,1)= factor_i + qfi_0(1,1)
          qfi(2,2)= factor_i + qfi_0(2,2)
          qfi(1,2)=            qfi_0(1,2)
          qfi(2,1)=            qfi_0(2,1)

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
      SUBROUTINE MTXDVB_ds(PROPS,STRAN,DMATX,STRSG, qbifu,beta_s,
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
      REAL*8  DAMAG,FLOAD,FACTM,FACTP,FACTQ,HPRIM,COEF1,PSIGMA,
     .        THETA,VECDOT,STRAX(9),ESCUR(*),esold(*),MULAM,lambl
	real*8 CARTD(2,*),KBETA,qfi_inv(2,2),qfi(2,2),cartd_b(4,2),
     .       nx,ny,n_m(4,2),b(4,2),dete,be(4,4),ben(4,4),cen(4,4),
     .       qbifu,beta_s(2),coe,nc(4,2),nxx,nyy,qeb(2),
     .       derbn,derbs,coe1,coe2,bn,bs
 
      INTEGER KLOAD,ISTR,JSTR,KSTR,LSTRE,nsol_b,implicito
C
C***Recover elastic properties 
C 
      YOUNG=PROPS(11)
      POISS=PROPS(12)
      MULAM=YOUNG/2.D0/(1.D0+POISS)
      LAMBL=YOUNG*POISS/(1.D0+POISS)/(1.D0-2.D0*POISS)
	KBETA= 1E-12/young
      QBIFU=ESCUR(4)
C
      implicito=impli
C
C***Compute elastic matrix 
C
      CALL MTELAS(DMATX) 
c
c***
c

csp c      if(fload.eq.0.d0) then
csp         DMATX(1,1)=FACTQ*DMATX(1,1)
csp         DMATX(1,2)=FACTQ*DMATX(1,2)
csp         DMATX(1,4)=FACTQ*DMATX(1,4)
csp         DMATX(2,1)=FACTQ*DMATX(2,1)
csp         DMATX(2,2)=FACTQ*DMATX(2,2)
csp         DMATX(2,4)=FACTQ*DMATX(2,4)
csp         DMATX(4,1)=FACTQ*DMATX(4,1)
csp         DMATX(4,2)=FACTQ*DMATX(4,2)
csp         DMATX(4,4)=FACTQ*DMATX(4,4)
csp         DMATX(3,3)=FACTQ*DMATX(3,3)
cspc	   return
cspc	endif
csp
csp
      if(fload.eq.0.d0) return
csp
csp
cspcq      if(0) then
	if (fload.eq.-1) then
c         DMATX(1,1)=FACTQ*DMATX(1,1)
c         DMATX(1,2)=FACTQ*DMATX(1,2)
c         DMATX(1,4)=FACTQ*DMATX(1,4)
c         DMATX(2,1)=FACTQ*DMATX(2,1)
c         DMATX(2,2)=FACTQ*DMATX(2,2)
c         DMATX(2,4)=FACTQ*DMATX(2,4)
c         DMATX(4,1)=FACTQ*DMATX(4,1)
c         DMATX(4,2)=FACTQ*DMATX(4,2)
c         DMATX(4,4)=FACTQ*DMATX(4,4)
cspc         DMATX(3,3)=FACTQ*DMATX(3,3)



cc         CALL POSI2D(NSTR1,STREF,STPOS)

             stpos(1)=stref(1)
             stpos(2)=stref(2)
             stpos(3)=stref(3)
             stpos(4)=stref(4)





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

c      nsol_b=INT(ESOLD(8))
c      IF(nsol_b.EQ.0) nsol_b=1

	cartd_b(1,1)=CARTD(1,nsol_b)
	cartd_b(3,1)=CARTD(2,nsol_b)
	cartd_b(2,2)=CARTD(2,nsol_b)
	cartd_b(3,2)=CARTD(1,nsol_b)

	n_m(1,1)=nx
	n_m(3,1)=ny
	n_m(2,2)=ny
	n_m(3,2)=nx
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

            qfi(1,1)= factp +QFI(1,1)
            qfi(2,2)= factp +QFI(2,2)
            qfi(1,2)=        QFI(1,2)
            qfi(2,1)=        QFI(2,1)
c***************************

              qeb(1)= beta_s(1)
              qeb(2)= beta_s(2)
              derbn= beta_s(1)
	        derbs= beta_s(2)

c***************************

        	  if (impli. and.(fload.eq.1.or.fload.eq.-2)) then
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
C***********************************************************************
      SUBROUTINE CRIBIS_sp(STREF,PROPS,ANGLN,sigmacri,t1_t2)
C***********************************************************************
C
C****THIS ROUTINE DETERMINES IF THE CRITICAL CONDITIONS AT BIFURCATION
C            =====   LARGE STRAIN PROBLEM  :: DAMAGE  ========
C
C    Input variables:
C      
C       STREF(NSTR1)       Effective stresses at Gauss point NGAUL
C       STRAN(NSTR1)       Total strains at Gauss point NGAUL
C       PROPS = Properties     
C       RINTV = Internal varible describing the size of the 
C               elastic region   
C       QFUNR = Scalar variable  (q)
C
C    Output variables: 
C
C       ANGLN              Localization angle
C       HPCRI              Critical hardening modulus
C
C***********************************************************************
      IMPLICIT NONE
C
      INCLUDE 'problem.om'  ! -->  NPROP
      INCLUDE 'element.om'  ! -->  NHIST, NSTR1, NGAUL
      INCLUDE 'propert.om'  ! -->  PH,PA
      INCLUDE 'auxiliar.om' ! -->  IGAUL
      INCLUDE 'interval.om' ! -->  IITER,DTIME
C
      REAL*8   STREF(nstr1), DUMMY, angln(*),sigmacri,PROPS(*)
      REAL*8   QFUNR,TAUXX,TAUYY,TAUXY,TAUZZ
      REAL*8   ADDTAU,REDTAU,TAU1,TAU2,VALUE1,VALUE2,ETT,ALFA,TWOPI,
     .         MULAM,LAMBL,THETAP,ANGLCH
C     
      real*8  A,B,N1C,N2C, ANS1,ANS2,NN(2),ALFA2, DENOM,STPOS(4),
     .        TAU1P,TAU2P,t1_t2

      PARAMETER (TWOPI=6.283185307179586D0)
C     
      if (ncrit.eq.78) then
	   YOUNG=PROPS(11)
         POISS=PROPS(12)
         MULAM=YOUNG/2.D0/(1.D0+POISS)
         LAMBL=YOUNG*POISS/(1.D0+POISS)/(1.D0-2.D0*POISS)
         A= 1.D0/MULAM
         B=(LAMBL+MULAM)/(MULAM*(LAMBL+2.D0*MULAM))
	endif

      IF(NCRIT.EQ.79.or.ncrit.eq.78) THEN ! ONLY TRACTION MODEL


ccc         CALL POSI2D(NSTR1,STREF,STPOS)


         CALL VECASI(NSTR1,STREF,STPOS) ! STPOS=STREF



      ELSE
         CALL VECASI(NSTR1,STREF,STPOS) ! STPOS=STREF
      ENDIF 
C     
C     -----------------------
C        POSITIVE STRESS
C        --------------
         TAUXX=STPOS(1)
         TAUYY=STPOS(2)
         TAUXY=STPOS(3)
         TAUZZ=STPOS(4)
         ADDTAU=TAUXX+TAUYY
         REDTAU=TAUXX-TAUYY
         TAU1P=0.5*ADDTAU+DSQRT(0.25D0*REDTAU*REDTAU+TAUXY*TAUXY)
         TAU2P=0.5*ADDTAU-DSQRT(0.25D0*REDTAU*REDTAU+TAUXY*TAUXY)
C        TOTAL STRESS
C        --------------
         TAUXX=STREF(1)
         TAUYY=STREF(2)
         TAUXY=STREF(3)
         TAUZZ=STREF(4)
         ADDTAU=TAUXX+TAUYY
         REDTAU=TAUXX-TAUYY
         TAU1=0.5*ADDTAU+DSQRT(0.25D0*REDTAU*REDTAU+TAUXY*TAUXY)
         TAU2=0.5*ADDTAU-DSQRT(0.25D0*REDTAU*REDTAU+TAUXY*TAUXY)
         IF(TAU1.EQ.0.0D0.AND.TAU2.EQ.0.0D0) THEN !Null stres.=def.values
            ANGLN(1)= 0.D0
            ANGLN(2)= 0.D0   
            sigmacri = TAU1
	      t1_t2=1.d0
            RETURN
         ELSE IF(TAU1.EQ.TAU2) THEN ! Default values
            ANGLN(1)= 0.D0
            ANGLN(2)= 0.D0
            sigmacri = TAU1
	      t1_t2=1.d0
            RETURN   
         ENDIF
         IF(TAU1P.EQ.0.0D0.AND.TAU2P.EQ.0.0D0) THEN !Null stres.=def.values
            ANGLN(1)= 0.D0
            ANGLN(2)= 0.D0   
            sigmacri = TAU1P
            RETURN
         ENDIF

	   t1_t2=tau1/tau2
         THETAP= - 0.5D0*DATAN2(2.0D0*TAUXY,(TAUXX-TAUYY))
C
         if (ncrit.eq.78) then
            N1C = (A*(TAU1*TAU1P-TAU2*TAU2P)
     .           - B*(TAU1*TAU2P+TAU2*TAU1P-2.D0*TAU2*TAU2P))
     .           /(2.D0*B*(TAU1*TAU1P-(TAU1*TAU2P+TAU2*TAU1P)+TAU2*TAU2P))
C
            if ((N1C.LE.1.D0).AND.(N1C.GE.0.D0))THEN
               N2C=1.D0-N1C
               ALFA=DATAN(DSQRT(N2C/N1C))
               ALFA2=DATAN(-DSQRT(N2C/N1C))
            else
               if (TAU1.GT.0.D0)THEN
                  ALFA=0.D0 
                  ALFA2=ALFA
               else
                  ALFA2=TWOPI/4.D0
                  ALFA=ALFA2
               endIF
            endIF
	   else
            alfa=0.0d0        !direcciones principales
 	      alfa2=0.0d0
         endif
	            
         ANS1=ALFA-THETAP
         ANS2=ALFA2-THETAP

         if(ANS1.GT.(TWOPI/4.D0)) THEN
            ANS1=ANS1-TWOPI/2.D0
         endIF
         if(ANS1.LT.-TWOPI/4.D0) THEN
            ANS1=ANS1+TWOPI/2.D0
         endIF
        
         if(ANS2.GT.TWOPI/4.D0) THEN
            ANS2=ANS2-TWOPI/2.D0
         endIF
         if(ANS2.LT.-TWOPI/4.D0) THEN
            ANS2=ANS2+TWOPI/2.D0
         endIF
         
         sigmacri = TAU1
C
         ANGLN(1)=ANS1
         ANGLN(2)=ANS2

      END
C
C***********************************************************************
      SUBROUTINE CRIBIS_sp1(STREF,PROPS,ANGLN,sigmacri)
C***********************************************************************
C
C****THIS ROUTINE DETERMINES IF THE CRITICAL CONDITIONS AT BIFURCATION
C            =====   LARGE STRAIN PROBLEM  :: DAMAGE  ========
C
C    Input variables:
C      
C       STREF(NSTR1)       Effective stresses at Gauss point NGAUL
C       STRAN(NSTR1)       Total strains at Gauss point NGAUL
C       PROPS = Properties     
C       RINTV = Internal varible describing the size of the 
C               elastic region   
C       QFUNR = Scalar variable  (q)
C
C    Output variables: 
C
C       ANGLN              Localization angle
C       HPCRI              Critical hardening modulus
C
C***********************************************************************
      IMPLICIT NONE
C
      INCLUDE 'problem.om'  ! -->  NPROP
      INCLUDE 'element.om'  ! -->  NHIST, NSTR1, NGAUL
      INCLUDE 'propert.om'  ! -->  PH,PA
      INCLUDE 'auxiliar.om' ! -->  IGAUL
      INCLUDE 'interval.om' ! -->  IITER,DTIME
C
      REAL*8   STREF(nstr1), DUMMY, angln(*), PROPS(*)
      REAL*8   TAUXX,TAUYY,TAUXY,TAUZZ
      REAL*8   ADDTAU,REDTAU,TAU1,TAU2,VALUE1,VALUE2,ETT,ALFA,TWOPI,pimed,
     .         MULAM,LAMBL,THETAP,ANGLCH,sigmacri

C     
      real*8  A,B,N1C,N2C,
     .        ANS1,ANS2,NN(2),ALFA2, RCUAD,DENOM,STPOS(4),
     .        TAU1P,TAU2P,N1CUA,N2CUA

      INTEGER  LTRAK
      PARAMETER (TWOPI=6.283185307179586D0,pimed=1.57079632679490d0)
c
      YOUNG=PROPS(11)
      POISS=PROPS(12)
      MULAM=YOUNG/2.D0/(1.D0+POISS)
      LAMBL=YOUNG*POISS/(1.D0+POISS)/(1.D0-2.D0*POISS)

      IF(NTYPE.EQ.1) THEN
         LAMBL= 2.D0*LAMBL*MULAM/(2.D0*MULAM+LAMBL)
      ENDIF
C     
c      IF(NCRIT.EQ.79) THEN ! ONLY TRACTION MODEL
c         CALL POSI2D(NSTR1,STREF,STPOS)
c      ELSE
         CALL VECASI(NSTR1,STREF,STPOS) ! STPOS=STREF
c      ENDIF 
      A= 1.D0/MULAM
      B=(LAMBL+MULAM)/(MULAM*(LAMBL+2.D0*MULAM))
C     
C     COMPUTE INCLINATION OF THE NORMAL ANGLE
C
C        POSITIVE STRESS
C        --------------
         TAUXX=STPOS(1)
         TAUYY=STPOS(2)
         TAUXY=STPOS(3)
         TAUZZ=STPOS(4)
         ADDTAU=TAUXX+TAUYY
         REDTAU=TAUXX-TAUYY
         TAU1P=0.5*ADDTAU+DSQRT(0.25D0*REDTAU*REDTAU+TAUXY*TAUXY)
         TAU2P=0.5*ADDTAU-DSQRT(0.25D0*REDTAU*REDTAU+TAUXY*TAUXY)
C        TOTAL STRESS
C        --------------
         TAUXX=STREF(1)
         TAUYY=STREF(2)
         TAUXY=STREF(3)
         TAUZZ=STREF(4)
         ADDTAU=TAUXX+TAUYY
         REDTAU=TAUXX-TAUYY
         TAU1=0.5*ADDTAU+DSQRT(0.25D0*REDTAU*REDTAU+TAUXY*TAUXY)
         TAU2=0.5*ADDTAU-DSQRT(0.25D0*REDTAU*REDTAU+TAUXY*TAUXY)

         sigmacri = TAU1

         IF(TAU1.EQ.0.0D0.AND.TAU2.EQ.0.0D0) THEN !Null stres.=def.values
c            ANGLN(1)= 0.D0
c            ANGLN(2)= 0.D0
            ANGLN(1)= pimed
            ANGLN(2)= pimed
            RETURN
         ELSE IF(TAU1.EQ.TAU2) THEN ! Default values
c            ANGLN(1)= 0.D0
c            ANGLN(2)= 0.D0
            ANGLN(1)= pimed
            ANGLN(2)= pimed
            RETURN   
         ENDIF
         IF(TAU1P.EQ.0.0D0.AND.TAU2P.EQ.0.0D0) THEN !Null stres.=def.values
c            ANGLN(1)= 0.D0
c            ANGLN(2)= 0.D0   
            ANGLN(1)= pimed
            ANGLN(2)= pimed   
            RETURN
         ENDIF

         THETAP= - 0.5D0*DATAN2(2.0D0*TAUXY,(TAUXX-TAUYY))
C
         N1C = (A*(TAU1*TAU1P-TAU2*TAU2P)
     .        - B*(TAU1*TAU2P+TAU2*TAU1P-2.D0*TAU2*TAU2P))
     .        /(2.D0*B*(TAU1*TAU1P-(TAU1*TAU2P+TAU2*TAU1P)+TAU2*TAU2P))
C
         if ((N1C.LE.1.D0).AND.(N1C.GE.0.D0))THEN
            N2C=1.D0-N1C
            ALFA = DATAN(DSQRT(N2C/N1C))
            ALFA2 = DATAN(-DSQRT(N2C/N1C))
         else
            if (TAU1.GT.0.D0)THEN
               ALFA=0.D0 
               ALFA2=ALFA
            else
               ALFA2=TWOPI/4.D0
               ALFA=ALFA2
            endIF
         endIF

ccc   direcciones principales
c        alfa=0.0d0
c	   alfa2=0.0d0
ccc   
         ANS1=ALFA-THETAP
         ANS2=ALFA2-THETAP

         if(ANS1.GT.(TWOPI/4.D0)) THEN
            ANS1=ANS1-TWOPI/2.D0
         endIF
         if(ANS1.LT.-TWOPI/4.D0) THEN
            ANS1=ANS1+TWOPI/2.D0
         endIF
        
         if(ANS2.GT.TWOPI/4.D0) THEN
            ANS2=ANS2-TWOPI/2.D0
         endIF
         if(ANS2.LT.-TWOPI/4.D0) THEN
            ANS2=ANS2+TWOPI/2.D0
         endIF 
C 
         ANGLN(1)=ANS1
         ANGLN(2)=ANS2
      return
	end
C
C***********************************************************************
      SUBROUTINE ten_max(strsg,PROPS,nx,ny,mx,my,tmm,tzz)
C***********************************************************************
      IMPLICIT NONE
C
      INCLUDE 'problem.om'  ! -->  NPROP
      INCLUDE 'element.om'  ! -->  NHIST, NSTR1, NGAUL
      INCLUDE 'propert.om'  ! -->  PH,PA
      INCLUDE 'auxiliar.om' ! -->  IGAUL
      INCLUDE 'interval.om' ! -->  IITER,DTIME
C
      REAL*8   strsg(nstr1), PROPS(*),nx,ny,mx,my,tracc(2)
      real*8   sig_nn,sig_nm,tmm,tzz
	integer  iban
C     
      tracc(1)= nx*strsg(1)+ ny*strsg(3)
      tracc(2)= nx*strsg(3)+ ny*strsg(2)
	sig_nn  = nx*tracc(1)+ ny*tracc(2)
	sig_nm  = mx*tracc(1)+ my*tracc(2)

      strsg(1) = sig_nn*nx*nx +  sig_nm*(nx*mx+nx*mx) + tmm*mx*mx
      strsg(2) = sig_nn*ny*ny +  sig_nm*(ny*my+ny*my) + tmm*my*my
      strsg(3) = sig_nn*ny*nx +  sig_nm*(nx*my+ny*mx) + tmm*my*mx
      strsg(4) = tzz

      return
	end
c
c**********************************************************************
      subroutine nsolitario(nx,ny,CARTD,nsoli,ndime)
c**********************************************************************
c
      real*8 CARTD(NDIME,*),nx,ny,prod(3),maxprod
	integer nsoli

      prod(1)=nx*CARTD(1,1)+ ny*CARTD(2,1)
      prod(2)=nx*CARTD(1,2)+ ny*CARTD(2,2)
      prod(3)=nx*CARTD(1,3)+ ny*CARTD(2,3)

      prod(1)=prod(1)/dsqrt(CARTD(1,1)**2+ CARTD(2,1)**2)
      prod(2)=prod(2)/dsqrt(CARTD(1,2)**2+ CARTD(2,2)**2)
      prod(3)=prod(3)/dsqrt(CARTD(1,3)**2+ CARTD(2,3)**2)

      nsoli   = 1; 
      maxprod = dabs(prod(1)); 
      if(dabs(prod(2)).gt.maxprod)then
	      nsoli = 2
		  maxprod=dabs(prod(2)) 
	endif
      if(dabs(prod(3)).gt.maxprod)then
	      nsoli = 3
      endif
      if(prod(nsoli).lt.0.d0) then
	   nx=-nx
	   ny=-ny
	endif
      return
	end
