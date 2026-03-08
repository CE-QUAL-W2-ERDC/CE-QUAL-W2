
!***********************************************************************************************************************************
!**                                            S U B R O U T I N E   W I T H D R A W A L                                          **
!***********************************************************************************************************************************

SUBROUTINE WITHDRAWAL
  USE GLOBAL; USE GEOMC; USE TVDC; USE SELWC; USE LOGICC
  USE MAIN, ONLY: DERIVED_CALC,CDN,TDGON,JSG,NNSG,NDO,JWD,EA,SYSTDG,NN2,NDGP,O2DG_DER,TDG_DER, theta_strt, wet_well, STR_FLOW_PROF_EQN
  USE modSYSTDG, ONLY: GTNAME, UPDATE_TDGC, SYSTDG_TDG
  USE SCREENC, ONLY: JDAY
  USE KINETIC, ONLY: CAC
  IMPLICIT NONE
  REAL :: HSWT,HSWB,ELR,WSEL,ELSTR,COEF,RATIO,HT,RHOFT,DLRHOT,HB,RHOFB,DLRHOB,VSUM,DLRHOMAX,HWDT,HWDB,ELWD,TEMPEST,ESTRTEST,QSUMJS
  REAL(R8) :: FRACV,QSUMWD
  REAL(R8) :: dosat, n2sat   ! cb 11/7/17
  INTEGER  :: K,JS,KSTR,KTOP,KBOT,KWD,JJWD  
  LOGICAL  :: jspair
  REAL(R8) :: elmv,dbot,dtop,hwz,z1estr,y1estr
  REAL(R8) :: vupsum, vdnsum, vupavg, vdnavg, Qjsk
  INTEGER  :: kmv, npair, jup, jdn, kcnt, ii, jj 
RETURN

!***********************************************************************************************************************************
!**                                             D O W N S T R E A M   W I T H D R A W A L                                         **
!***********************************************************************************************************************************

ENTRY DOWNSTREAM_WITHDRAWAL (JS)

! Variable initialization

  HSWT = 0.0; HSWB = 0.0; VNORM = 0.0; QNEW = 0.0

! Water surface elevation

  ELR  = SINA(JB)*DLX(ID)*0.5
  WSEL = ELWS(ID)-ELR                   !EL(KT,ID)-Z(ID)*COSA(JB)

! Structure layer

  DO K=KT,KB(ID)
    IF (EL(K,ID)-ELR < ESTR(JS,JB)) EXIT
  END DO
  KSTR = MAX(K-1,KT)
  KSTR = MIN(KSTR,KB(ID))

! Initial withdrawal limits

  KTOP = MAX(KTSW(JS,JB),KT)
  IF (KSTR < KTOP) KTOP = KSTR
  KBOT = MIN(KBSW(JS,JB),KB(ID))
  IF (KBOT <= KT .AND. KBOT /= KB(ID)) KBOT = KT+1
  IF (KBOT > KB(ID)) KBOT = KB(ID)
  ELSTR = ESTR(JS,JB)
  IF (ESTR(JS,JB) <= EL(KB(ID)+1,ID+1)-ELR) THEN
    KSTR  = KB(ID)
    ELSTR = EL(KB(ID),ID)-ELR
  END IF
  IF (ESTR(JS,JB) > EL(KT,ID)-ELR) ELSTR = WSEL
  IF (KBSW(JS,JB) < KSTR) THEN
    KSTR  = KTOP          ! KT
    ELSTR = EL(KTOP,ID)   !WSEL
  END IF

! Boundary interference

  COEF = 1.0
  IF (WSEL-(EL(KBOT,ID)-ELR) /= 0.0) THEN   ! SR 11/2021
    RATIO = (ELSTR-(EL(KBOT,ID)-ELR))/(WSEL-(EL(KBOT,ID)-ELR))
    IF (RATIO < 0.1 .OR. RATIO > 0.9) COEF = 2.0
  END IF

! Withdrawal zone above structure

  DO K=KSTR-1,KTOP,-1

!** Density frequency

    HT    = (EL(K,ID)-ELR)-ELSTR
    RHOFT = MAX(SQRT((ABS(RHO(K,ID)-RHO(KSTR,ID)))/(HT*RHO(KSTR,ID)+NONZERO)*G),NONZERO)

!** Thickness

    IF (POINT_SINK(JS,JB)) THEN
      HSWT = (COEF*QSTR(JS,JB)/RHOFT)**0.333333
    ELSE
      HSWT = SQRT(2.0*COEF*QSTR(JS,JB)/(WSTR(JS,JB)*RHOFT))
    END IF
    IF (HT >= HSWT) THEN
      KTOP = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELSTR+HSWT) < WSEL) THEN
    DLRHOT = ABS(RHO(KSTR,ID)-RHO(KTOP,ID))
  ELSE IF (WSEL == ELSTR) THEN
    DLRHOT = NONZERO
  ELSE
    DLRHOT = ABS(RHO(KSTR,ID)-RHO(KT,ID))*HSWT/(WSEL-ELSTR)
  END IF
  DLRHOT = MAX(DLRHOT,NONZERO)

! Withdrawal zone below structure

  DO K=KSTR+1,KBOT

!** Density frequency

    HB    = ELSTR-(EL(K,ID)-ELR)
    RHOFB = MAX(SQRT((ABS(RHO(K,ID)-RHO(KSTR,ID)))/(HB*RHO(KSTR,ID)+NONZERO)*G),NONZERO)

!** Thickness

    IF (POINT_SINK(JS,JB)) THEN
      HSWB = (COEF*QSTR(JS,JB)/RHOFB)**0.333333
    ELSE
      HSWB = SQRT(2.0*COEF*QSTR(JS,JB)/(WSTR(JS,JB)*RHOFB))
    END IF
    IF (HB >= HSWB) THEN
      KBOT = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELSTR-HSWB) > EL(KBOT+1,ID)) THEN
    DLRHOB = ABS(RHO(KSTR,ID)-RHO(KBOT,ID))
  ELSE IF ((EL(KBOT+1,ID)-ELR) == ELSTR) THEN                                                                          !SR 03/24/13
    DLRHOB = NONZERO                                                                                                   !SR 03/24/13
  ELSE
    DLRHOB = ABS(RHO(KSTR,ID)-RHO(KBOT,ID))*HSWB/(ELSTR-(EL(KBOT+1,ID)-ELR))
  END IF
  DLRHOB = MAX(DLRHOB,NONZERO)

! Velocity profile

  VSUM     = 0.0
!  DLRHOMAX = MAX(DLRHOT,DLRHOB,1.0E-10)                      ! GH 1/31/08
  DO K=KTOP,KBOT
!    VNORM(K) = ABS(1.0-((RHO(K,ID)-RHO(KSTR,ID))/DLRHOMAX)**2)*BHR2(K,ID)
 	   IF(K.GT.KSTR)THEN
       DLRHOMAX = MAX(DLRHOB,1.0E-10)                          !GH 1/31/08
       ELSE
       DLRHOMAX = MAX(DLRHOT,1.0E-10)                          !GH 1/31/08
       ENDIF
     VNORM(K) = 1.0-((RHO(K,ID)-RHO(KSTR,ID))/DLRHOMAX)**2
 	 IF(VNORM(K).GT.1.0) VNORM(K)=1.0                         !GH 1/31/08
	 IF(VNORM(K).LT.0.0) VNORM(K)=0.0                         !GH 1/31/08
	 VNORM(K)=VNORM(K)*BHR2(K,ID)
     VSUM     = VSUM+VNORM(K)
  END DO

! OUTFLOWS
  QSUMJS=0.0                                                  ! SW 7/30/09
  TAVG(JS,JB)=0.0                                                    ! CB 5/12/10
  IF(CONSTITUENTS)CAVG(JS,JB,CN(1:NAC))=0.0
  IF(DERIVED_CALC)CDAVG(JS,JB,CDN(1:NACD(JW),JW))=0.0
  IF(VSUM==0.0)THEN
      WRITE(WRN,'(A,F12.3,A,I5,A,I5,A,I5,A,E12.4,A)')'DOWNSTREAM WITHDRAWAL: VSUM=0.0 on JDAY:',JDAY,' KTOP:',KTOP,' KBOT:',KBOT,' KSTR:',KSTR,' DLRHOMAX:',DLRHOMAX,' SET TO EQUAL WITHDRAWALS WITH DEPTH'
      VSUM=1.0
      DO K=KTOP,KBOT
      VNORM(K)=1.0/(KTOP-KBOT+1)
      ENDDO
  ENDIF
  
  DO K=KTOP,KBOT
    QNEW(K)    = (VNORM(K)/VSUM)*QSTR(JS,JB)
    QOUT(K,JB) =  QOUT(K,JB)+QNEW(K)
    QDSW(K,ID) =  QDSW(K,ID)+QNEW(K)      ! For layer-specific output; ID redefined for SP/PI/PU/GT    !SR 12/19/2022
    TAVG(JS,JB)=TAVG(JS,JB)+QNEW(K)*T2(K,ID)                  ! SW 7/30/09
    IF(CONSTITUENTS)CAVG(JS,JB,CN(1:NAC))=CAVG(JS,JB,CN(1:NAC))+QNEW(K)*C2(K,ID,CN(1:NAC))  
    IF(DERIVED_CALC)CDAVG(JS,JB,CDN(1:NACD(JW),JW))=CDAVG(JS,JB,CDN(1:NACD(JW),JW))+QNEW(K)*CD(K,ID,CDN(1:NACD(JW),JW))
    QSUMJS=QSUMJS+QNEW(K)
  END DO
IF(QSUMJS.GT.0.0)THEN
  TAVG(JS,JB)=TAVG(JS,JB)/QSUMJS
  IF(CONSTITUENTS)then                    ! cb 1/16/13
    CAVG(JS,JB,CN(1:NAC))=CAVG(JS,JB,CN(1:NAC))/QSUMJS
    if(tdgon)then
     IF (nnsg==1) THEN                                ! nnsg==1 is gate flow
           !
           ! systedg
           IF (SYSTDG) THEN                          
               IF(GTNAME(jsg)) THEN                   
                  CALL  UPDATE_TDGC (0,palt(id),jsg,tavg(js,jb),cavg(js,jb,NDO))    
                  CALL  UPDATE_TDGC (1,palt(id),jsg,tavg(js,jb),cavg(js,jb,NN2))  
                  CALL  UPDATE_TDGC (2,palt(id),jsg,tavg(js,jb),cavg(js,jb,NDGP))
               ELSE
                  call total_dissolved_gas (0,palt(id),nnsg,jsg,tavg(js,jb),cavg(js,jb,NDO))    
                  call total_dissolved_gas (1,palt(id),nnsg,jsg,tavg(js,jb),cavg(js,jb,NN2))     ! n2 GAS
                  call total_dissolved_gas (2,palt(id),nnsg,jsg,tavg(js,jb),cavg(js,jb,NDGP)) 
               END IF
           ELSE
               call total_dissolved_gas (0,palt(id),nnsg,jsg,tavg(js,jb),cavg(js,jb,NDO))    
               call total_dissolved_gas (1,palt(id),nnsg,jsg,tavg(js,jb),cavg(js,jb,NN2))     ! n2 GAS
               call total_dissolved_gas (2,palt(id),nnsg,jsg,tavg(js,jb),cavg(js,jb,NDGP))
           END IF
           !
      ELSE
      call total_dissolved_gas (0,palt(id),nnsg,jsg,tavg(js,jb),cavg(js,jb,NDO))    
      call total_dissolved_gas (1,palt(id),nnsg,jsg,tavg(js,jb),cavg(js,jb,NN2))     ! n2 GAS
      call total_dissolved_gas (2,palt(id),nnsg,jsg,tavg(js,jb),cavg(js,jb,NDGP))
      END IF
    end if
  end if
  IF(DERIVED_CALC)then                    ! cb 1/16/13
    CDAVG(JS,JB,CDN(1:NACD(JW),JW))=CDAVG(JS,JB,CDN(1:NACD(JW),JW))/QSUMJS
    !if(tdgon)then                  ! cb 11/6/17
      !cdavg(js,jb,16)  = (cavg(js,jb,ndo)/exp(7.7117-1.31403*(log(tavg(js,jb)+45.93)))*palt(id))*100.0 
      dosat=exp(7.7117-1.31403*(log(tavg(js,jb)+45.93)))*palt(id)
      cdavg(js,jb,O2DG_DER)=(cavg(js,jb,ndo)/dosat)*100.0 
      !IF(ngctdg /= 0)THEN
      If(CAC(NN2)== '      ON') THEN
          EA = DEXP(2.3026D0*(7.5D0*TDEW(JW)/(TDEW(JW)+237.3D0)+0.6609D0))*0.001316   ! in mm Hg   0.0098692atm=7.5006151mmHg  
          !cdavg(js,jb,NDC)  = (cavg(js,jb,NGN2)/(1.5568D06*0.79*(PALT(ID)-EA)*(1.8816D-5 - 4.116D-7 * Tavg(js,jb) + 4.6D-9 * Tavg(js,jb)**2)))*100.0    ! SW 10/27/15      
          n2sat=1.5568D06*0.79*(PALT(ID)-EA)*(1.8816D-5 - 4.116D-7 * Tavg(js,jb) + 4.6D-9 * Tavg(js,jb)**2)
          cdavg(js,jb,TDG_DER)  = 100.*(0.79*(cavg(js,jb,NN2)/n2sat) + 0.21*(cavg(js,jb,ndo)/dosat))
      ELSE IF(CAC(NDGP)== '      ON') THEN
          cdavg(js,jb,TDG_DER)  = cavg(js,jb,NDGP)/palt(id)*100.0
      END IF
      !ENDIF
    !end if
  end if
ELSE
  TAVG(JS,JB)=-99.0
  IF(CONSTITUENTS)CAVG(JS,JB,CN(1:NAC))=-99.0
  IF(DERIVED_CALC)CDAVG(JS,JB,CDN(1:NACD(JW),JW))=-99.0
END IF

! Inactive layers and total outflow

  IF (JS == NST) THEN
    WHERE (QOUT(:,JB) == 0.0) U(:,ID) = 0.0
  END IF
RETURN

!***********************************************************************************************************************************
!**                                             D O W N S T R E A M   W I T H D R A W A L  N E W                                 **
!***********************************************************************************************************************************
! 06/2023, Updated withdrawal calculation
! includes update of withdrawal limits, location of maximum velocity, overlap of withdrawal zones and single wet well option
ENTRY DOWNSTREAM_WITHDRAWAL_MULTI

! Water surface elevation

  ELR  = SINA(JB)*DLX(ID)*0.5
  WSELs = ELWS(ID)-ELR
  
   DO jsn=1,NSTR(JB)
  
  ! Structure layer

      DO K=KT,KB(ID)
        IF (EL(K,ID)-ELR < ESTR(jsn,JB)) EXIT
      END DO
      KSTRs(jsn) = MAX(K-1,KT)
      KSTRs(jsn) = MIN(KSTRs(jsn),KB(ID))

    ! Initial withdrawal limits

      KTOPs(jsn) = MAX(KTSW(jsn,JB),KT)
      IF (KSTRs(jsn) < KTOPs(jsn)) KTOPs(jsn) = KSTRs(jsn)
      KBOTs(jsn) = MIN(KBSW(jsn,JB),KB(ID))
      IF (KBOTs(jsn) <= KT .AND. KBOTs(jsn) /= KB(ID)) KBOTs(jsn) = KT+1
      IF (KBOTs(jsn) > KB(ID)) KBOTs(jsn) = KB(ID)
      ELSTRs(jsn) = ESTR(jsn,JB)
      IF (ESTR(jsn,JB) <= EL(KB(ID)+1,ID+1)-ELR) THEN
        KSTRs(jsn)  = KB(ID)
        ELSTRs(jsn) = EL(KB(ID),ID)-ELR
      END IF
      IF (ESTR(jsn,JB) > EL(KT,ID)-ELR) ELSTRs(jsn) = WSELs
      IF (KBSW(jsn,JB) < KSTRs(jsn)) THEN
        KSTRs(jsn)  = KT
        ELSTRs(jsn) = WSELs
      END IF
  end do

  if(wet_well(jb) .and. nstr(jb)>1)call single_wet_well   ! redistributes single well port flows using Howington(1990) research

  HSWTs = 0.0; HSWBs = 0.0; VNORMs = 0.0
  
  DO jsn=1,NSTR(JB)
    IF (QSTR(jsn,JB) /= 0.0) THEN
      
      tbi=.false.; bbi=.false.
      thetast=theta_STRT(jsn,jb)

    ! Boundary interference 

      COEF = 1.0
      IF (WSELs-(EL(KBOTs(jsn),ID)-ELR) /= 0.0) THEN   ! SR 11/2021
        RATIO = (ELSTRs(jsn)-(EL(KBOTs(jsn),ID)-ELR))/(WSELs-(EL(KBOTs(jsn),ID)-ELR))
        IF (RATIO < 0.1 .OR. RATIO > 0.9) COEF = 2.0
      END IF

    ! Withdrawal zone above structure

      DO K=KSTRs(jsn)-1,KTOPs(jsn),-1

    !** Density frequency

        HT    = (EL(K,ID)-ELR)-ELSTRs(jsn)
        RHOFT = MAX(SQRT((ABS(RHO(K,ID)-RHO(KSTRs(jsn),ID)))/(HT*RHO(KSTRs(jsn),ID)+NONZERO)*G),NONZERO)

    !** Thickness

        IF (POINT_SINK(jsn,JB)) THEN
          !HSWT = (COEF*QSTR(jsn,JB)/RHOFT)**0.333333
          HSWTs(jsn) = (COEF*QSTR(jsn,JB)/(RHOFT*thetast))**0.333333    
        ELSE
          !HSWT = SQRT(2.0*COEF*QSTR(jsn,JB)/(WSTR(jsn,JB)*RHOFT))
            HSWTs(jsn) = SQRT(2.0*COEF*QSTR(jsn,JB)/(WSTR(jsn,JB)*RHOFT*thetast))
        END IF
        IF (HT >= HSWTs(jsn)) THEN
          KTOPs(jsn) = K; EXIT
        END IF
      END DO

    ! Withdrawal zone below structure

      DO K=KSTRs(jsn)+1,KBOTs(jsn)

    !** Density frequency

        HB    = ELSTRs(jsn)-(EL(K,ID)-ELR)
        RHOFB = MAX(SQRT((ABS(RHO(K,ID)-RHO(KSTRs(jsn),ID)))/(HB*RHO(KSTRs(jsn),ID)+NONZERO)*G),NONZERO)

    !** Thickness

        IF (POINT_SINK(jsn,JB)) THEN
          !HSWB = (COEF*QSTR(jsn,JB)/RHOFB)**0.333333
            HSWBs(jsn) = (COEF*QSTR(jsn,JB)/(RHOFB*thetast))**0.333333  
        ELSE
          !HSWB = SQRT(2.0*COEF*QSTR(jsn,JB)/(WSTR(jsn,JB)*RHOFB))
            HSWBs(jsn) = SQRT(2.0*COEF*QSTR(jsn,JB)/(WSTR(jsn,JB)*RHOFB*thetast))  
        END IF
        IF (HB >= HSWBs(jsn)) THEN
          KBOTs(jsn) = K; EXIT
        END IF
      END DO  
  
    ! checking for top and bottom boundary interference
       dtop=WSELs-elr-ELSTRs(jsn)
       if(hswts(jsn) >  dtop)tbi=.true.
       dbot=ELSTRs(jsn)-(EL(KB(id)+1,ID)-ELR)
       if(hswbs(jsn) >  dbot)bbi=.true.
   
       ! recalculating of withdrawal zone if surface boundary interference occurs on one side
       if(tbi .and. .not. bbi)then  ! top boundary interference
           bwz=dtop
           qwflow=QSTR(jsn,JB)
          call WITHDRAWAL_ZONE
      
          KBOTs(jsn) = MIN(KBSW(jsn,JB),KB(ID))   ! need to recalculate KBOT for new HSWB
          IF (KBOTs(jsn) <= KT .AND. KBOTs(jsn) /= KB(ID)) KBOTs(jsn) = KT+1
          IF (KBOTs(jsn) > KB(ID)) KBOTs(jsn) = KB(ID)
           DO K=KSTRs(jsn)+1,KBOTs(jsn)
             HB    = ELSTRs(jsn)-(EL(K,ID)-ELR)
             IF (HB >= HSWBs(jsn)) THEN
               KBOTs(jsn) = K; EXIT
            END IF
          END DO 
       end if
       if(bbi .and. .not. tbi)then  ! bottom boundary interference
           bwz=dbot
           qwflow=QSTR(jsn,JB)
          call WITHDRAWAL_ZONE
          KTOPs(jsn)= MAX(KTSW(jsn,JB),KT)      ! need to recalculate KTOP for new HSWT
          IF (KSTRs(jsn) < KTOPs(jsn)) KTOPs(jsn) = KSTRs(jsn)
           DO K=KSTRs(jsn)-1,KTOPs(jsn),-1
             HT    = (EL(K,ID)-ELR)-ELSTRs(jsn)
             IF (HT >= HSWTs(jsn)) THEN
               KTOPs(jsn) = K; EXIT
             END IF
          END DO
       end if  
  
    ! Determining location of maximum velocity if surface and bottom interference do NOT occur
       elmvs(jsn)=ELSTRs(jsn)
       kmvs(jsn)=kstrs(jsn)
       if(.not. bbi .and. .not. tbi)then  ! eqn only applies if no boundary interference occurs
         z1estr=hswbs(jsn)
         hwz=hswbs(jsn)+hswts(jsn)
         y1estr= hwz * (sin(1.57*z1estr/hwz))**2
         elmvs(jsn)=y1estr+(ELSTRs(jsn)-hswbs(jsn))
         DO K=KSTRs(jsn)-1,KTOPs(jsn),-1
            if(elmvs(jsn) <= (EL(K+1,ID)-ELR))then
              kmvs(jsn)=k+1
              exit
            end if
         end do
        end if
  
    ! Reference density above maximum velocity

      IF ((ELSTRs(jsn)+HSWTs(jsn)) < WSELs) THEN
        !DLRHOT = ABS(RHO(KSTR,ID)-RHO(KTOP,ID))
        DLRHOT = ABS(RHO(kmvs(jsn),ID)-RHO(KTOPs(jsn),ID))
      ELSE IF (WSELs == ELSTRs(jsn)) THEN
        DLRHOT = NONZERO
       ELSE
        !DLRHOT = ABS(RHO(KSTR,ID)-RHO(KT,ID))*HSWT/(WSEL-ELSTR)
            DLRHOT = ABS(RHO(kmvs(jsn),ID)-RHO(KT,ID))*HSWTs(jsn)/(WSELs-elmvs(jsn))
      END IF
      DLRHOT = MAX(DLRHOT,NONZERO)
    
        ! Reference density below maximum velocity
  
      IF ((ELSTRs(jsn)-HSWBs(jsn)) > EL(KBOTs(jsn)+1,ID)) THEN
        !DLRHOB = ABS(RHO(KSTR,ID)-RHO(KBOT,ID))
        DLRHOB = ABS(RHO(kmvs(jsn),ID)-RHO(KBOTs(jsn),ID))
      ELSE IF ((EL(KBOTs(jsn)+1,ID)-ELR) == ELSTRs(jsn)) THEN                                                                          !SR 03/24/13
        DLRHOB = NONZERO                                                                                                   !SR 03/24/13
        ELSE
        !DLRHOB = ABS(RHO(KSTR,ID)-RHO(KBOT,ID))*HSWB/(ELSTR-(EL(KBOT+1,ID)-ELR))
        DLRHOB = ABS(RHO(kmvs(jsn),ID)-RHO(KBOTs(jsn),ID))*HSWBs(jsn)/(elmvs(jsn)-(EL(KBOTs(jsn)+1,ID)-ELR))
      END IF
      DLRHOB = MAX(DLRHOB,NONZERO)
    
    ! Velocity profile

      VSUMs(jsn)     = 0.0

      DO K=KTOPs(jsn),KBOTs(jsn)
 	       !IF(K.GT.KSTR)THEN
          IF(K.GT.kmvs(jsn))THEN
           DLRHOMAX = MAX(DLRHOB,1.0E-10)                          !GH 1/31/08
           ELSE
           DLRHOMAX = MAX(DLRHOT,1.0E-10)                          !GH 1/31/08
           ENDIF
         !VNORM(K) = 1.0-((RHO(K,ID)-RHO(KSTR,ID))/DLRHOMAX)**2
           VNORMs(K,jsn) = 1.0-((RHO(K,ID)-RHO(kmvs(jsn),ID))/DLRHOMAX)**2
 	     IF(VNORMs(K,jsn).GT.1.0) VNORMs(K,jsn)=1.0                         !GH 1/31/08
	     IF(VNORMs(K,jsn).LT.0.0) VNORMs(K,jsn)=0.0                         !GH 1/31/08
	     if(str_flow_prof_eqn(jsn,JB) == 1) VNORMs(K,jsn)=VNORMs(K,jsn)*BHR2(K,ID)
         VSUMs(jsn)     = VSUMs(jsn)+VNORMs(K,jsn)
      END DO
        
       veljs(:,jsn)=0
       DO K=KTOPs(jsn),KBOTs(jsn)
        Qjsk    = (VNORMs(K,jsn)/VSUMs(jsn))*QSTR(jsn,JB)
        veljs(k,jsn)=Qjsk/BHR2(K,ID)
      END DO
         
    end if
  end do
  
  ! finding pairs of outlets with withdrawal zones that overlap
  npair=0
  portup=.false.
  portdn=.false.
  if(nstr(jb) > 1)then
    do jsn=1,nstr(jb)
     IF (QSTR(jsn,JB) /= 0.0) THEN
      do jj=1,nstr(jb)
       if(qstr(jj,jb) /= 0.0 .and. kstrs(jj) < kstrs(jsn))then
        if(portup(jsn))then
          if( kstrs(jj) > kstrs(jsup(jsn)))then
            jsup(jsn)=jj 
          end if
        else
          portup(jsn)=.true.
          jsup(jsn)=jj
         end if
       end if
       if(qstr(jj,jb) /= 0.0 .and. kstrs(jj) > kstrs(jsn))then
        if(portdn(jsn))then
           if(kstrs(jj) < kstrs(jsdn(jsn)))then
            jsdn(jsn)=jj
           end if
        else
          jsdn(jsn)=jj
          portdn(jsn)=.true.
        end if
       end if
      end do
     end if
    end do
    
    do jsn=1,nstr(jb)    ! determining how may "port pairs" exist with overlapping zones
     IF (QSTR(jsn,JB) /= 0.0) THEN
      if(portup(jsn))then
        jspair=.false.
        do jj=1,npair        ! making sure port pair has not been counted yet
          if(jsup(jsn) == jspairup(jj) .and. jsn== jspairdn(jj))then
            jspair=.true.
          end if
        end do
        if(.not. jspair)then
          if(elstrs(jsup(jsn)) - hswbs(jsup(jsn)) < elstrs(jsn) + hswts(jsn))then  ! testing to see if the withdrawal zones overlap
            do k=kt,kb(id)  ! checking to see if overlap zone incluesd layer center with velocity prediction
              if(EL(K,ID)-ELR > elstrs(jsup(jsn)) - hswbs(jsup(jsn)) .and. EL(K,ID)-ELR < elstrs(jsn) + hswts(jsn))then
                npair=npair+1
                jspairup(npair)=jsup(jsn)
                jspairdn(npair)=jsn
                exit
              end if
            end do
          end if
        end if
      end if
      if(portdn(jsn))then
        jspair=.false.
        do jj=1,npair        ! making sure port pair has not been counted yet
          if(jsdn(jsn) == jspairdn(jj) .and. jsn== jspairup(jj))then
            jspair=.true.
          end if
        end do
        if(.not. jspair)then
          if(elstrs(jsn) - hswbs(jsn) < elstrs(jsdn(jsn)) + hswts(jsdn(jsn)))then  ! testing to see if the withdrawal zones overlap  
            do k=kt,kb(id)  ! checking to see if overlap zone include layer center with velocity prediction
              if(EL(K,ID)-ELR > elstrs(jsn) - hswbs(jsn) .and. EL(K,ID)-ELR < elstrs(jsdn(jsn)) + hswts(jsdn(jsn)))then
                npair=npair+1
                jspairup(npair)=jsn
                jspairdn(npair)=jsdn(jsn)
                exit
              end if
            end do
          end if
        end if
      end if
     end if   
    end do
  end if
  
  do ii=1,npair
      jup=jspairup(ii)
      jdn=jspairdn(ii)
   
      ! finding average velocities of each port in overlap zone
      kcnt=0
      vupsum=0.0
      vdnsum=0.0
      do k=kt,kb(id)  
        if(EL(K,ID)-ELR > elstrs(jup) - hswbs(jup) .and. EL(K,ID)-ELR < elstrs(jdn) + hswts(jdn))then
           vupsum=vupsum+veljs(k,jup)
           vdnsum=vdnsum+veljs(k,jdn)
           kcnt=kcnt+1
        end if
      end do
      vupavg=vupsum/real(kcnt)
      vdnavg=vdnsum/real(kcnt)
      hport=elstrs(jup)-elstrs(jdn)
      hovlap=(elstrs(jdn) + hswts(jdn)) - (elstrs(jup) - hswbs(jup))
      
       ! updating upper port's updated lower wihtdrawal limits and normalized velocities
      KSTRshft=kstrs(jup)  
      zmax=ELSTRs(jup)-(EL(KB(ID)+1,ID)-ELR)
      vozavg=vupavg
      kwzlim=kbots(jup)
      call WITHDRAWAL_ZONE_SHIFT
      if(wzsconv)hswbs(jup)=hswbs(jup)+zshift   ! shift solution will not converge if there are no vertical density differences
    
      KBOTs(jup) = MIN(KBSW(jup,JB),KB(ID))
      IF (KBOTs(jup) <= KT .AND. KBOTs(jup) /= KB(ID)) KBOTs(jup) = KT+1
      IF (KBOTs(jup) > KB(ID)) KBOTs(jup) = KB(ID)
      DO K=KSTRs(jup)+1,KBOTs(jup)
        HB    = ELSTRs(jup)-(EL(K,ID)-ELR)
        IF (HB >= HSWBs(jup)) THEN
          KBOTs(jup) = K; EXIT
        END IF
      END DO 
      
      ! Reference density above maximum velocity

      IF ((ELSTRs(jup)+HSWTs(jup)) < wsels) THEN
        DLRHOT = ABS(RHO(kmvs(jup),ID)-RHO(KTOPs(jup),ID))
      ELSE IF (wsels == ELSTRs(jup)) THEN
        DLRHOT = NONZERO
       ELSE
            DLRHOT = ABS(RHO(kmvs(jup),ID)-RHO(KT,ID))*HSWTs(jup)/(wsels-elmvs(jup))
      END IF
      DLRHOT = MAX(DLRHOT,NONZERO)
    
      ! Reference density below maximum velocity
  
      IF ((ELSTRs(jup)-HSWBs(jup)) > EL(KBOTs(jup)+1,ID)) THEN
        DLRHOB = ABS(RHO(kmvs(jup),ID)-RHO(KBOTs(jup),ID))
      ELSE IF ((EL(KBOTs(jup)+1,ID)-ELR) == ELSTRs(jup)) THEN                                                                          !SR 03/24/13
        DLRHOB = NONZERO                                                                                                   !SR 03/24/13
        ELSE
        DLRHOB = ABS(RHO(kmvs(jup),ID)-RHO(KBOTs(jup),ID))*HSWBs(jup)/(elmvs(jup)-(EL(KBOTs(jup)+1,ID)-ELR))
      END IF
      DLRHOB = MAX(DLRHOB,NONZERO)
      
      VSUMs(jup)     = 0.0
      DO K=KTOPs(jup),KBOTs(jup)
          IF(K.GT.kmvs(jup))THEN
           DLRHOMAX = MAX(DLRHOB,1.0E-10)  
           ELSE
           DLRHOMAX = MAX(DLRHOT,1.0E-10)                        
           ENDIF
           VNORMs(K,jup) = 1.0-((RHO(K,ID)-RHO(kmvs(jup),ID))/DLRHOMAX)**2
 	     IF(VNORMs(K,jup).GT.1.0) VNORMs(K,jup)=1.0                         
	     IF(VNORMs(K,jup).LT.0.0) VNORMs(K,jup)=0.0                        
	     if(str_flow_prof_eqn(jup,JB) == 1) VNORMs(K,jup)=VNORMs(K,jup)*BHR2(K,ID)
         VSUMs(jup)     = VSUMs(jup)+VNORMs(K,jup)
      END DO
      
      ! finding lower port's upper wihtdrawal limits and normalized velocities
      KSTRshft=kstrs(jdn)    
      zmax=elws(id)-ELSTRs(jdn)
      vozavg=vdnavg
      kwzlim=ktops(jdn)
      call WITHDRAWAL_ZONE_SHIFT
      if(wzsconv)hswts(jdn)=hswts(jdn)+zshift   ! shift solution will not converge if there are no vertical density differences
      
       KTOPs(jdn) = MAX(KTSW(jdn,JB),KT)
      IF (KSTRs(jdn) < KTOPs(jdn)) KTOPs(jdn) = KSTRs(jdn)
      DO K=KSTRs(jdn)-1,KTOPs(jdn),-1
        HT    = (EL(K,ID)-ELR)-ELSTRs(jdn)
        IF (HT >= HSWTs(jdn)) THEN
          KTOPs(jdn) = K; EXIT
        END IF
      END DO
      
      ! Reference density above maximum velocity

      IF ((ELSTRs(jdn)+HSWTs(jdn)) < wsels) THEN
        DLRHOT = ABS(RHO(kmvs(jdn),ID)-RHO(KTOPs(jdn),ID))
      ELSE IF (wsels == ELSTRs(jdn)) THEN
        DLRHOT = NONZERO
       ELSE
            DLRHOT = ABS(RHO(kmvs(jdn),ID)-RHO(KT,ID))*HSWTs(jdn)/(wsels-elmvs(jdn))
      END IF
      DLRHOT = MAX(DLRHOT,NONZERO)
    
      ! Reference density below maximum velocity
  
      IF ((ELSTRs(jdn)-HSWBs(jdn)) > EL(KBOTs(jdn)+1,ID)) THEN
        DLRHOB = ABS(RHO(kmvs(jdn),ID)-RHO(KBOTs(jdn),ID))
      ELSE IF ((EL(KBOTs(jdn)+1,ID)-ELR) == ELSTRs(jdn)) THEN                                                                          !SR 03/24/13
        DLRHOB = NONZERO                                                                                                   !SR 03/24/13
        ELSE
        DLRHOB = ABS(RHO(kmvs(jdn),ID)-RHO(KBOTs(jdn),ID))*HSWBs(jdn)/(elmvs(jdn)-(EL(KBOTs(jdn)+1,ID)-ELR))
      END IF
      DLRHOB = MAX(DLRHOB,NONZERO)
      
      VSUMs(jdn)     = 0.0
      DO K=KTOPs(jdn),KBOTs(jdn)
          IF(K.GT.kmvs(jdn))THEN
           DLRHOMAX = MAX(DLRHOB,1.0E-10)  
           ELSE
           DLRHOMAX = MAX(DLRHOT,1.0E-10)                        
           ENDIF
           VNORMs(K,jdn) = 1.0-((RHO(K,ID)-RHO(kmvs(jdn),ID))/DLRHOMAX)**2
 	     IF(VNORMs(K,jdn).GT.1.0) VNORMs(K,jdn)=1.0                         
	     IF(VNORMs(K,jdn).LT.0.0) VNORMs(K,jdn)=0.0                        
	     if(str_flow_prof_eqn(jdn,JB) == 1) VNORMs(K,jdn)=VNORMs(K,jdn)*BHR2(K,ID)
         VSUMs(jdn)     = VSUMs(jdn)+VNORMs(K,jdn)
      END DO
          
  end do

  !QNEW = 0.0
  DO jsn=1,NSTR(JB)
    IF (QSTR(jsn,JB) /= 0.0) THEN
    ! OUTFLOWS
      QSUMJS=0.0                                                  ! SW 7/30/09
      TAVG(jsn,JB)=0.0                                            ! CB 5/12/10
      QNEW = 0.0                                                  !ZZ 9/19/24
      IF(CONSTITUENTS)CAVG(jsn,JB,CN(1:NAC))=0.0
      IF(DERIVED_CALC)CDAVG(jsn,JB,CDN(1:NACD(JW),JW))=0.0
      DO K=KTOPs(jsn),KBOTs(jsn)
        QNEW(K)    = (VNORMs(K,jsn)/VSUMs(jsn))*QSTR(jsn,JB)
        QOUT(K,JB) =  QOUT(K,JB)+QNEW(K)
        QDSW(K,ID) =  QDSW(K,ID)+QNEW(K)                            ! For layer-specific output; ID redefined for SP/PI/PU/GT    !SR 12/19/2022
        TAVG(jsn,JB)=TAVG(jsn,JB)+QNEW(K)*T2(K,ID)                  ! SW 7/30/09
        IF(CONSTITUENTS)CAVG(jsn,JB,CN(1:NAC))=CAVG(jsn,JB,CN(1:NAC))+QNEW(K)*C2(K,ID,CN(1:NAC))  
        IF(DERIVED_CALC)CDAVG(jsn,JB,CDN(1:NACD(JW),JW))=CDAVG(jsn,JB,CDN(1:NACD(JW),JW))+QNEW(K)*CD(K,ID,CDN(1:NACD(JW),JW))
        QSUMJS=QSUMJS+QNEW(K)
      END DO
    IF(QSUMJS.GT.0.0)THEN
      TAVG(jsn,JB)=TAVG(jsn,JB)/QSUMJS
      IF(CONSTITUENTS)then                    ! cb 1/16/13
        CAVG(jsn,JB,CN(1:NAC))=CAVG(jsn,JB,CN(1:NAC))/QSUMJS
        if(tdgon)then
         IF (nnsg==1) THEN                                ! nnsg==1 is gate flow
               !
               ! systedg
               IF (SYSTDG) THEN                          
                   IF(GTNAME(jsg)) THEN                   
                      CALL  UPDATE_TDGC (0,palt(id),jsg,tavg(jsn,jb),cavg(jsn,jb,NDO))    
                      CALL  UPDATE_TDGC (1,palt(id),jsg,tavg(jsn,jb),cavg(jsn,jb,NN2))  
                      CALL  UPDATE_TDGC (2,palt(id),jsg,tavg(jsn,jb),cavg(jsn,jb,NDGP))
                   ELSE
                      call total_dissolved_gas (0,palt(id),nnsg,jsg,tavg(jsn,jb),cavg(jsn,jb,NDO))    
                      call total_dissolved_gas (1,palt(id),nnsg,jsg,tavg(jsn,jb),cavg(jsn,jb,NN2))     ! n2 GAS
                      call total_dissolved_gas (2,palt(id),nnsg,jsg,tavg(jsn,jb),cavg(jsn,jb,NDGP)) 
                   END IF
               ELSE
                   call total_dissolved_gas (0,palt(id),nnsg,jsg,tavg(jsn,jb),cavg(jsn,jb,NDO))    
                   call total_dissolved_gas (1,palt(id),nnsg,jsg,tavg(jsn,jb),cavg(jsn,jb,NN2))     ! n2 GAS
                   call total_dissolved_gas (2,palt(id),nnsg,jsg,tavg(jsn,jb),cavg(jsn,jb,NDGP))
               END IF
               !
          ELSE
          call total_dissolved_gas (0,palt(id),nnsg,jsg,tavg(jsn,jb),cavg(jsn,jb,NDO))    
          call total_dissolved_gas (1,palt(id),nnsg,jsg,tavg(jsn,jb),cavg(jsn,jb,NN2))     ! n2 GAS
          call total_dissolved_gas (2,palt(id),nnsg,jsg,tavg(jsn,jb),cavg(jsn,jb,NDGP))
          END IF
        end if
      end if
      IF(DERIVED_CALC)then                    
        CDAVG(jsn,JB,CDN(1:NACD(JW),JW))=CDAVG(jsn,JB,CDN(1:NACD(JW),JW))/QSUMJS
       
          dosat=exp(7.7117-1.31403*(log(tavg(jsn,jb)+45.93)))*palt(id)
          cdavg(jsn,jb,O2DG_DER)=(cavg(jsn,jb,ndo)/dosat)*100.0 
          
          If(CAC(NN2)== '      ON') THEN
              EA = DEXP(2.3026D0*(7.5D0*TDEW(JW)/(TDEW(JW)+237.3D0)+0.6609D0))*0.001316   ! in mm Hg   0.0098692atm=7.5006151mmHg  
              
              n2sat=1.5568D06*0.79*(PALT(ID)-EA)*(1.8816D-5 - 4.116D-7 * Tavg(jsn,jb) + 4.6D-9 * Tavg(jsn,jb)**2)
              cdavg(jsn,jb,TDG_DER)  = 100.*(0.79*(cavg(jsn,jb,NN2)/n2sat) + 0.21*(cavg(jsn,jb,ndo)/dosat))
          ELSE IF(CAC(NDGP)== '      ON') THEN
              cdavg(jsn,jb,TDG_DER)  = cavg(jsn,jb,NDGP)/palt(id)*100.0
          END IF
          
      end if
    ELSE
      TAVG(jsn,JB)=-99.0
      IF(CONSTITUENTS)CAVG(jsn,JB,CN(1:NAC))=-99.0
      IF(DERIVED_CALC)CDAVG(jsn,JB,CDN(1:NACD(JW),JW))=-99.0
    END IF

    ! Inactive layers and total outflow

      IF (jsn == NST) THEN
        WHERE (QOUT(:,JB) == 0.0) U(:,ID) = 0.0
      END IF
  end if
end do
   
RETURN


!***********************************************************************************************************************************
!**                                             D O W N S T R E A M   W I T H D R A W A L  ESTIMATE                               **
!***********************************************************************************************************************************

ENTRY DOWNSTREAM_WITHDRAWAL_ESTIMATE(JS,TEMPEST,ESTRTEST)

! VARIABLE INITIALIZATION

  HSWT = 0.0; HSWB = 0.0; VNORM = 0.0; QNEW = 0.0

! Water surface elevation

  ELR  = SINA(JB)*DLX(ID)*0.5
  WSEL = ELWS(ID)-ELR                   !EL(KT,ID)-Z(ID)*COSA(JB)                                                     !SR 12/19/2022

! Structure layer

  DO K=KT,KB(ID)
    IF (EL(K,ID)-ELR < estrtest) EXIT
  END DO
  KSTR = MAX(K-1,KT)
  KSTR = MIN(KSTR,KB(ID))

! Initial withdrawal limits

  KTOP = MAX(KTSW(JS,JB),KT)
  IF (KSTR < KTOP) KTOP = KSTR
  KBOT = MIN(KBSW(JS,JB),KB(ID))
  IF (KBOT <= KT .AND. KBOT /= KB(ID)) KBOT = KT+1
  IF (KBOT > KB(ID)) KBOT = KB(ID)                                                                                     !SW 06/03/02
  ELSTR = ESTRTEST
  IF (ESTRTEST <= EL(KB(ID)+1,ID+1)-ELR) THEN                                                                       !SW 10/17/01
    KSTR  = KB(ID)
    ELSTR = EL(KB(ID),ID)-ELR                                                                                          !SW 10/17/01
  END IF
  IF (ESTRTEST > EL(KT,ID)-ELR) ELSTR = WSEL
  IF (KBSW(JS,JB) < KSTR) THEN
    KSTR  = KTOP         !KT
    ELSTR = EL(KTOP,ID)  !WSEL                                                                                                       !SW 10/05/00
  END IF

! Boundary interference

  COEF = 1.0
  IF (WSEL-(EL(KBOT,ID)-ELR) /= 0.0) THEN    ! SR 11/2021
    RATIO = (ELSTR-(EL(KBOT,ID)-ELR))/(WSEL-(EL(KBOT,ID)-ELR))                                                         !SW 10/17/01
    IF (RATIO < 0.1 .OR. RATIO > 0.9) COEF = 2.0
  END IF

! Withdrawal zone above structure

  DO K=KSTR-1,KTOP,-1

!** Density frequency

    HT    = (EL(K,ID)-ELR)-ELSTR
    RHOFT = MAX(SQRT((ABS(RHO(K,ID)-RHO(KSTR,ID)))/(HT*RHO(KSTR,ID)+NONZERO)*G),NONZERO)

!** Thickness

    IF (POINT_SINK(JS,JB)) THEN
      HSWT = (COEF*QSTR(JS,JB)/RHOFT)**0.333333
    ELSE
      HSWT = SQRT(2.0*COEF*QSTR(JS,JB)/(WSTR(JS,JB)*RHOFT))
    END IF
    IF (HT >= HSWT) THEN
      KTOP = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELSTR+HSWT) < WSEL) THEN
    DLRHOT = ABS(RHO(KSTR,ID)-RHO(KTOP,ID))
  ELSE IF (WSEL == ELSTR) THEN
    DLRHOT = NONZERO
  ELSE
    DLRHOT = ABS(RHO(KSTR,ID)-RHO(KT,ID))*HSWT/(WSEL-ELSTR)
  END IF
  DLRHOT = MAX(DLRHOT,NONZERO)

! Withdrawal zone below structure

  DO K=KSTR+1,KBOT

!** Density frequency

    HB    = ELSTR-(EL(K,ID)-ELR)                                                                                       !SW 10/17/01
    RHOFB = MAX(SQRT((ABS(RHO(K,ID)-RHO(KSTR,ID)))/(HB*RHO(KSTR,ID)+NONZERO)*G),NONZERO)

!** Thickness

    IF (POINT_SINK(JS,JB)) THEN
      HSWB = (COEF*QSTR(JS,JB)/RHOFB)**0.333333
    ELSE
      HSWB = SQRT(2.0*COEF*QSTR(JS,JB)/(WSTR(JS,JB)*RHOFB))
    END IF
    IF (HB >= HSWB) THEN
      KBOT = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELSTR-HSWB) > EL(KBOT+1,ID)) THEN
    DLRHOB = ABS(RHO(KSTR,ID)-RHO(KBOT,ID))
  ELSE IF ((EL(KBOT+1,ID)-ELR) == ELSTR) THEN                                                                          !SR 03/24/13
    DLRHOB = NONZERO                                                                                                   !SR 03/24/13
  ELSE
    DLRHOB = ABS(RHO(KSTR,ID)-RHO(KBOT,ID))*HSWB/(ELSTR-(EL(KBOT+1,ID)-ELR))                                           !SW 10/17/01
  END IF
  DLRHOB = MAX(DLRHOB,NONZERO)

! Velocity profile

  VSUM     = 0.0
!  DLRHOMAX = MAX(DLRHOT,DLRHOB,1.0E-10)
  DO K=KTOP,KBOT
  IF(K.GT.KSTR)THEN
       DLRHOMAX = MAX(DLRHOB,1.0E-10)                          !GH 1/31/08
       ELSE
       DLRHOMAX = MAX(DLRHOT,1.0E-10)                          !GH 1/31/08
       ENDIF
       VNORM(K) = 1.0-((RHO(K,ID)-RHO(KSTR,ID))/DLRHOMAX)**2
     IF(VNORM(K).GT.1.0) VNORM(K)=1.0                          !GH 1/31/08
	 IF(VNORM(K).LT.0.0) VNORM(K)=0.0                          !GH 1/31/08
	 VNORM(K)=VNORM(K)*BHR2(K,ID)
     VSUM= VSUM+VNORM(K)
  END DO

! Outflows
    IF(VSUM==0.0)THEN
      WRITE(WRN,'(A,F12.3,A,I5,A,I5,A,I5,A,E12.4,A)')'DOWNSTREAM WITHDRAWAL ESTIMATE: VSUM=0.0 on JDAY:',JDAY,' KTOP:',KTOP,' KBOT:',KBOT,' KSTR:',KSTR,' DLRHOMAX:',DLRHOMAX,' SET TO EQUAL WITHDRAWALS WITH DEPTH'
      VSUM=1.0
      DO K=KTOP,KBOT
      VNORM(K)=1.0/(KTOP-KBOT+1)
      ENDDO
  ENDIF


  tempest=0.0
  DO K=KTOP,KBOT
    tempest=tempest+t2(k,id)*(VNORM(K)/VSUM)*QSTR(JS,JB)
  END DO

  if(qstr(js,jb).gt.0.0)tempest=tempest/qstr(js,jb)

RETURN
!***********************************************************************************************************************************
!**                                                L A T E R A L   W I T H D R A W A L                                            **
!***********************************************************************************************************************************

ENTRY LATERAL_WITHDRAWAL

! Variable initialization

  VNORM = 0.0; QSW(:,JWD) = 0.0; HWDT = 0.0; HWDB = 0.0

! Structure layer

  K = KT
  DO K=KT,KB(I)
    IF (EL(K,I) < EWD(JWD)) EXIT
  END DO
  KWD = MAX(K-1,KT)
  KWD = MIN(KWD,KB(I))

! Initial withdrawal limits

  KTOP = MAX(KTWD(JWD),KT)
  IF (KWD < KTOP) KTOP = KWD
  KBOT = MIN(KBWD(JWD),KB(I))
  IF (KBOT <= KT .AND. KB(I) /= KBOT) KBOT = KT+1
  IF (KBOT > KB(I)) KBOT = KB(I)
  ELWD = EWD(JWD)
  IF (EWD(JWD) <= EL(KB(I)+1,I)) THEN
    KWD  = KB(I)
    ELWD = EL(KB(I),I)
  END IF
  IF (EWD(JWD) > EL(KT,I)) ELWD = EL(KT,I)
  IF (KBWD(JWD) < KWD) THEN
    KWD  = KTOP        ! KT  !SW 8/13/2024
    ELWD = EL(KTOP,I)  ! EL(KT,I)
  END IF

! Boundary interference

  COEF = 1.0
  IF (KT /= KBOT) THEN
    RATIO = (ELWD-EL(KBOT,I))/(EL(KT,I)-EL(KBOT,I))
    IF (RATIO < 0.1 .OR. RATIO > 0.9) COEF = 2.0
  END IF

! Withdrawal zone above structure

  DO K=KWD-1,KTOP,-1

!** Density frequency

    HT    = EL(K,I)-ELWD
    RHOFT = MAX(SQRT((ABS(RHO(K,I)-RHO(KWD,I)))/(HT*RHO(KWD,I)+NONZERO)*G),NONZERO)

!** Thickness

    HWDT = (COEF*QWD(JWD)/RHOFT)**0.333333
    IF (HT >= HWDT) THEN
      KTOP = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELWD+HWDT) < EL(KT,I)) THEN
    DLRHOT = ABS(RHO(KWD,I)-RHO(KTOP,I))
  ELSE IF (EL(KT,I) == ELWD) THEN
    DLRHOT = NONZERO
  ELSE
    DLRHOT = ABS(RHO(KWD,I)-RHO(KT,I))*HWDT/(EL(KT,I)-ELWD)
  END IF
  DLRHOT = MAX(DLRHOT,NONZERO)

! Withdrawal zone below structure

  DO K=KWD+1,KBOT

!** Density frequency

    HB    = ELWD-EL(K,I)
    RHOFB = MAX(SQRT((ABS(RHO(K,I)-RHO(KWD,I)))/(HB*RHO(KWD,I)+NONZERO)*G),NONZERO)

!** Thickness

    HWDB = (COEF*QWD(JWD)/RHOFB)**0.333333
    IF (HB >= HWDB) THEN
      KBOT = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELWD-HWDB) > EL(KBOT+1,I)) THEN
    DLRHOB = ABS(RHO(KWD,I)-RHO(KBOT,I))
  ELSE IF (EL(KBOT+1,I) == ELWD) THEN                                                                                  !SR 03/24/13
    DLRHOB = NONZERO                                                                                                   !SR 03/24/13
  ELSE
    DLRHOB = ABS(RHO(KWD,I)-RHO(KBOT,I))*HWDB/(ELWD-EL(KBOT+1,I))
  END IF
  DLRHOB = MAX(DLRHOB,NONZERO)

! Velocity profile

  VSUM     = 0.0
!  DLRHOMAX = MAX(DLRHOT,DLRHOB,1.0E-10)                                                                             ! SW 1/24/05
  DO K=KTOP,KBOT
!    VNORM(K) = ABS(1.0-((RHO(K,I)-RHO(KWD,I))/DLRHOMAX)**2)*BHR2(K,I)
 	   IF(K.GT.KWD)THEN
       DLRHOMAX = MAX(DLRHOB,1.0E-10)                          !GH 1/31/08
       ELSE
       DLRHOMAX = MAX(DLRHOT,1.0E-10)                          !GH 1/31/08
       ENDIF
     VNORM(K) = 1.0-((RHO(K,I)-RHO(KWD,I))/DLRHOMAX)**2
 	 IF(VNORM(K).GT.1.0) VNORM(K)=1.0                         !GH 1/31/08
	 IF(VNORM(K).LT.0.0) VNORM(K)=0.0                         !GH 1/31/08
	 VNORM(K)=VNORM(K)*BHR2(K,I)
     VSUM     = VSUM+VNORM(K)
  END DO

! Outflows
  QSUMWD=0.0                                                  ! SW 7/30/09
  TAVGW(JWD)=0.0
  IF(CONSTITUENTS)CAVGW(JWD,CN(1:NAC))=0.0
  IF(DERIVED_CALC)CDAVGW(JWD,CDN(1:NACD(JW),JW))=0.0
  
    IF(VSUM==0.0)THEN
      WRITE(WRN,'(A,F12.3,A,I5,A,I5,A,I5,A,E12.4,A)')'LATERAL WITHDRAWAL: VSUM=0.0 on JDAY:',JDAY,' KTOP:',KTOP,' KBOT:',KBOT,' KWD:',KWD,' DLRHOMAX:',DLRHOMAX,' SET TO EQUAL WITHDRAWALS WITH DEPTH'
      VSUM=1.0
      DO K=KTOP,KBOT
      VNORM(K)=1.0/(KTOP-KBOT+1)
      ENDDO
  ENDIF

  DO K=KTOP,KBOT
    FRACV=(VNORM(K)/VSUM)
    QSW(K,JWD) = QSW(K,JWD)+FRACV*QWD(JWD)
    TAVGW(JWD)=TAVGW(JWD)+FRACV*QWD(JWD)*T2(K,I)                  ! SW 7/30/09
    IF(CONSTITUENTS)CAVGW(JWD,CN(1:NAC))=CAVGW(JWD,CN(1:NAC))+FRACV*QWD(JWD)*C2(K,I,CN(1:NAC))  
    IF(DERIVED_CALC)CDAVGW(JWD,CDN(1:NACD(JW),JW))=CDAVGW(JWD,CDN(1:NACD(JW),JW))+FRACV*QWD(JWD)*CD(K,I,CDN(1:NACD(JW),JW))
    QSUMWD=QSUMWD+FRACV*QWD(JWD)
  END DO
  ! Debug
  !if(qwd(jwd)>0.0 .and. qsumwd <= 0.0)then
  !    write(9575,'(A,f8.3,1x,i5,1x,i5,1x,i5,f8.4,1x,f8.4,1x,f10.2)')'JDAY, ktop, kbot, kwd, qwd, qsumwd, elwd:',JDAY, ktop, kbot, kwd, qwd, qsumwd, elwd
  !endif
  ! Debug
  IF(QSUMWD.GT.0.0)THEN
    TAVGW(JWD)=TAVGW(JWD)/QSUMWD               ! SW 7/30/09
    IF(CONSTITUENTS)then                       ! cb 1/16/13
      CAVGW(JWD,CN(1:NAC))=CAVGW(JWD,CN(1:NAC))/QSUMWD  
      if(tdgon)then
        IF (nnsg==1) THEN                                      ! systdg        nnsg==1 is gate flow
            IF (SYSTDG) THEN                                   ! systdg 
                IF (GTNAME(jsg)) THEN                          ! systdg 
                   CALL  UPDATE_TDGC (0,palt(i),jsg,tavgw(jwd),cavgw(jwd,NDO))         ! systdg 
                   CALL  UPDATE_TDGC(1,palt(i),jsg,tavgw(jwd),cavgw(jwd,NN2))          ! systdg 
                   CALL  UPDATE_TDGC(2,palt(i),jsg,tavgw(jwd),cavgw(jwd,NDGP))
                ELSE
                   call total_dissolved_gas (0,palt(i),nnsg,jsg,tavgw(jwd),cavgw(jwd,NDO))  
                   call total_dissolved_gas (1,palt(i),nnsg,jsg,tavgw(jwd),cavgw(jwd,NN2))     ! n2 GAS
                   call total_dissolved_gas (2,palt(i),nnsg,jsg,tavgw(jwd),cavgw(jwd,NDGP)) 
                END IF
            ELSE
        call total_dissolved_gas (0,palt(i),nnsg,jsg,tavgw(jwd),cavgw(jwd,NDO))  
        call total_dissolved_gas (1,palt(i),nnsg,jsg,tavgw(jwd),cavgw(jwd,NN2))     ! n2 GAS
        call total_dissolved_gas (2,palt(i),nnsg,jsg,tavgw(jwd),cavgw(jwd,NDGP))    
            END IF   
        ELSE
            call total_dissolved_gas (0,palt(i),nnsg,jsg,tavgw(jwd),cavgw(jwd,NDO))  
            call total_dissolved_gas (1,palt(i),nnsg,jsg,tavgw(jwd),cavgw(jwd,NN2))     ! n2 GAS
            call total_dissolved_gas (2,palt(i),nnsg,jsg,tavgw(jwd),cavgw(jwd,NDGP)) 
        END IF
      end if
    end if
    IF(DERIVED_CALC)then
      CDAVGW(JWD,CDN(1:NACD(JW),JW))=CDAVGW(JWD,CDN(1:NACD(JW),JW))/QSUMWD
      !if(tdgon)then                ! cb 11/6/17
        !cdavgw(jwd,O2DG_DER)  = (cavgw(jwd,ndo)/exp(7.7117-1.31403*(log(tavgw(jwd)+45.93)))*palt(i))*100.0       
        dosat=exp(7.7117-1.31403*(log(tavgw(jwd)+45.93)))*palt(i)
        cdavgw(jwd,O2DG_DER)=(cavgw(jwd,ndo)/dosat)*100.0 
        If(CAC(NN2)== '      ON') THEN
          EA = DEXP(2.3026D0*(7.5D0*TDEW(JW)/(TDEW(JW)+237.3D0)+0.6609D0))*0.001316   ! in mm Hg   0.0098692atm=7.5006151mmHg  
          !cdavgw(jwd,NDC)  = (cavgw(jwd,NGN2)/(1.5568D06*0.79*(PALT(I)-EA)*(1.8816D-5 - 4.116D-7 * Tavgw(jwd) + 4.6D-9 * Tavgw(jwd)**2)))*100.0    ! SW 10/27/15      
          n2sat=1.5568D06*0.79*(PALT(I)-EA)*(1.8816D-5 - 4.116D-7 * Tavgw(jwd) + 4.6D-9 * Tavgw(jwd)**2)
          cdavgw(jwd,TDG_DER)  = 100.*(0.79*(cavgw(jwd,NN2)/n2sat) + 0.21*(cavgw(jwd,ndo)/dosat))
        ELSE IF(CAC(NDGP)== '      ON') THEN
          cdavgw(jwd,TDG_DER)  =  cavgw(jwd,NDGP)/ palt(i) * 100.0
        END IF  
    end if
  ELSE
    TAVGW(JWD)=-99.0
    IF(CONSTITUENTS)CAVGW(JWD,CN(1:NAC))=-99.0 
    IF(DERIVED_CALC)CDAVGW(JWD,CDN(1:NACD(JW),JW))=-99.0
  ENDIF
  KTW(JWD) = KTOP
  KBW(JWD) = KBOT
  RETURN
!***********************************************************************************************************************************
!**                                                L A T E R A L   W I T H D R A W A L ESTIMATE                                   **
!***********************************************************************************************************************************

  ENTRY LATERAL_WITHDRAWAL_ESTIMATE (JJWD,TEMPEST,ESTRTEST)

! VARIABLE INITIALIZATION

  VNORM = 0.0; QSW(:,JJWD) = 0.0; HWDT = 0.0; HWDB = 0.0

! Structure layer

  K = KT
  DO K=KT,KB(I)
    IF (EL(K,I) < estrtest) EXIT
  END DO
  KWD = MAX(K-1,KT)
  KWD = MIN(KWD,KB(I))

! Initial withdrawal limits

  KTOP = MAX(KTWD(JJWD),KT)
  IF (KWD < KTOP) KTOP = KWD
  KBOT = MIN(KBWD(JJWD),KB(I))
  IF (KBOT <= KT .AND. KB(I) /= KBOT) KBOT = KT+1
  IF (KBOT > KB(I)) KBOT = KB(I)
  ELWD = ESTRTEST
  IF (ESTRTEST <= EL(KB(I)+1,I)) THEN
    KWD  = KB(I)
    ELWD = EL(KB(I),I)
  END IF
  IF (ESTRTEST > EL(KT,I)) ELWD = EL(KT,I)
  IF (KBWD(JJWD) < KWD) THEN
    KWD  = KTOP        ! KT
    ELWD = EL(KTOP,I)  ! EL(KT,I)
  END IF

! Boundary interference

  COEF = 1.0
  IF (KT /= KBOT) THEN
    RATIO = (ELWD-EL(KBOT,I))/(EL(KT,I)-EL(KBOT,I))
    IF (RATIO < 0.1 .OR. RATIO > 0.9) COEF = 2.0
  END IF

! Withdrawal zone above structure

  DO K=KWD-1,KTOP,-1

!** Density frequency

    HT    = EL(K,I)-ELWD
    RHOFT = MAX(SQRT((ABS(RHO(K,I)-RHO(KWD,I)))/(HT*RHO(KWD,I)+NONZERO)*G),NONZERO)

!** Thickness

    HWDT = (COEF*QWD(JJWD)/RHOFT)**0.333333
    IF (HT >= HWDT) THEN
      KTOP = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELWD+HWDT) < EL(KT,I)) THEN
    DLRHOT = ABS(RHO(KWD,I)-RHO(KTOP,I))
  ELSE IF (EL(KT,I) == ELWD) THEN
    DLRHOT = NONZERO
  ELSE
    DLRHOT = ABS(RHO(KWD,I)-RHO(KT,I))*HWDT/(EL(KT,I)-ELWD)
  END IF
  DLRHOT = MAX(DLRHOT,NONZERO)

! Withdrawal zone below structure

  DO K=KWD+1,KBOT

!** Density frequency

    HB    = ELWD-EL(K,I)
    RHOFB = MAX(SQRT((ABS(RHO(K,I)-RHO(KWD,I)))/(HB*RHO(KWD,I)+NONZERO)*G),NONZERO)

!** Thickness

    HWDB = (COEF*QWD(JJWD)/RHOFB)**0.333333
    IF (HB >= HWDB) THEN
      KBOT = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELWD-HWDB) > EL(KBOT+1,I)) THEN
    DLRHOB = ABS(RHO(KWD,I)-RHO(KBOT,I))
  ELSE IF (EL(KBOT+1,I) == ELWD) THEN                                                                                  !SR 03/24/13
    DLRHOB = NONZERO                                                                                                   !SR 03/24/13
  ELSE
    DLRHOB = ABS(RHO(KWD,I)-RHO(KBOT,I))*HWDB/(ELWD-EL(KBOT+1,I))
  END IF
  DLRHOB = MAX(DLRHOB,NONZERO)

! Velocity profile

  VSUM     = 0.0
!  DLRHOMAX = MAX(DLRHOT,DLRHOB,1.0E-10)                                                                             ! SW 1/24/05
  DO K=KTOP,KBOT
!    VNORM(K) = ABS(1.0-((RHO(K,I)-RHO(KWD,I))/DLRHOMAX)**2)*BHR2(K,I)
 	   IF(K.GT.KWD)THEN
       DLRHOMAX = MAX(DLRHOB,1.0E-10)                          !GH 1/31/08
       ELSE
       DLRHOMAX = MAX(DLRHOT,1.0E-10)                          !GH 1/31/08
       ENDIF
     VNORM(K) = 1.0-((RHO(K,I)-RHO(KWD,I))/DLRHOMAX)**2
 	 IF(VNORM(K).GT.1.0) VNORM(K)=1.0                          !GH 1/31/08
	 IF(VNORM(K).LT.0.0) VNORM(K)=0.0                          !GH 1/31/08
	 VNORM(K)=VNORM(K)*BHR2(K,I)
     VSUM = VSUM+VNORM(K)
  END DO

! Outflows
  
    IF(VSUM==0.0)THEN
      WRITE(WRN,'(A,F12.3,A,I5,A,I5,A,I5,A,E12.4,A)')'LATERAL WITHDRAWAL ESTIMATE: VSUM=0.0 on JDAY:',JDAY,' KTOP:',KTOP,' KBOT:',KBOT,' KWD:',KWD,' DLRHOMAX:',DLRHOMAX,' SET TO EQUAL WITHDRAWALS WITH DEPTH'
      VSUM=1.0
      DO K=KTOP,KBOT
      VNORM(K)=1.0/(KTOP-KBOT+1)
      ENDDO
  ENDIF


  TEMPEST = 0.0                                                                                                       !SR 12/19/2022
  DO K=KTOP,KBOT
    tempest=tempest+t2(k,i)*(VNORM(K)/VSUM)*QWD(JJWD)
  END DO
  if(qwd(Jjwd).gt.0.0)tempest=tempest/qwd(Jjwd)
  KTW(JJWD) = KTOP
  KBW(JJWD) = KBOT
  return

  END SUBROUTINE WITHDRAWAL

  
!***********************************************************************************************************************************
!**         W I T H D R A W A L    Z O N E                                                                                     **
!***********************************************************************************************************************************
subroutine WITHDRAWAL_ZONE
      
    USE GLOBAL; USE GEOMC; USE SELWC
    IMPLICIT NONE
    REAL :: elr
    INTEGER, PARAMETER ::JMAX=40    
    REAL :: X1, X2, FUNCVAL1, FUNCVAL2, XACC, FMID, FUNC1, RTBIS, DX, XMID
    INTEGER :: JJ,J  
                 
      wzconv=.false.
      ELR  = SINA(JB)*DLX(ID)*0.5
      
      ! FIRST, BRACKETING ROOT
      X1=0.001
      !X2=1.0
      if(bbi)then
        X2=elws(id)-ELSTRs(jsn)   ! finding zone above outlet
      else if(tbi)then
        X2=ELSTRs(jsn)-(EL(KB(ID)+1,ID)-ELR)
      end if
         
      CALL  WD_ZONE_EQN(x1,FUNCVAL1)     
      CALL  WD_ZONE_EQN(x2,FUNCVAL2)
                  
      DO JJ=1,JMAX
        IF(FUNCVAL1*FUNCVAL2 > 0.0)THEN
          IF(ABS(FUNCVAL1)  < ABS(FUNCVAL2))THEN
            X1=X1/2.0
            CALL  WD_ZONE_EQN(x1,FUNCVAL1)            
          ELSE
            X2=X2+1.5*(X2-X1)
            CALL  WD_ZONE_EQN(x2,FUNCVAL2)
          END IF                                          
        ELSE
          EXIT
        END IF      
      END DO                 
      
      ! FINDING ROOT BY BISECTION      
      XACC=0.01
      
      CALL  WD_ZONE_EQN(x2,FMID)  
      CALL  WD_ZONE_EQN(x1,FUNC1)      
  !    IF(FUNC1*FMID.GE.0.) PAUSE 'ROOT MUST BE BRACKETED IN RTBIS'
      IF(FUNC1.LT.0.)THEN
        RTBIS=X1
        DX=X2-X1
      ELSE
        RTBIS=X2
        DX=X1-X2
      ENDIF
      DO J=1,JMAX
        DX=DX*.5
        XMID=RTBIS+DX
        
        CALL  WD_ZONE_EQN(xmid,FMID)        
        IF(FMID.LE.0.)RTBIS=XMID
        IF(ABS(DX).LT.XACC .OR. FMID.EQ.0.)THEN          
            if(bbi)then
              hswts(jsn)=rtbis
            else if(tbi)then
              hswbs(jsn)=rtbis
            end if
          wzconv=.true.
          RETURN
        END IF
      END DO
      
      return
end subroutine WITHDRAWAL_ZONE

!***********************************************************************************************************************************
!**           W I T H D R A W A L    Z O N E     E Q U A T I O N                                                                  **
!***********************************************************************************************************************************
subroutine WD_ZONE_EQN(ddp,funcvalue2)    
    USE GLOBAL; USE GEOMC; USE SELWC
    IMPLICIT NONE
  
    REAL :: HT,HB,rhof,elr,ddp
    INTEGER :: K
    REAL    :: weqn_num, weqn_den
      
    REAL :: funcvalue2
    INTEGER :: JJ,J 
      
      ELR  = SINA(JB)*DLX(ID)*0.5
      HB=ddp-bwz
      if(bbi)then  ! bottom boundary interference
        DO K=KSTRs(jsn)-1,2,-1
          HT    = (EL(K,ID)-ELR)-ELSTRs(jsn)
          IF (HT >= HB) THEN
             EXIT
          END IF
        END DO
      end if
      
      if(tbi)then  ! top boundary interference
        DO K=KSTRs(jsn)+1,KB(ID)
          HT= ELSTRs(jsn)- (EL(K,ID)-ELR)
          IF (HB <= HT) THEN
            EXIT
          END IF
        END DO   
      end if
      
      RHOF = MAX(SQRT((ABS(RHO(K,ID)-RHO(KSTRs(jsn),ID)))/(HB*RHO(KSTRs(jsn),ID)+NONZERO)*G),NONZERO)  
      
      weqn_num= 1 +1/pi *sin((bwz/ddp)/(1-bwz/ddp)*pi)+(bwz/ddp)/(1-bwz/ddp)
      weqn_den= 1 + (bwz/ddp)/(1-bwz/ddp)
      funcvalue2= qwFLOW - ddp**3*rhof * 0.5 * weqn_num * thetast/pi  / (weqn_den)**3
      
      RETURN
end subroutine WD_ZONE_EQN

!***********************************************************************************************************************************
!**         W I T H D R A W A L    Z O N E    S H I F T                                                                           **
!***********************************************************************************************************************************
subroutine WITHDRAWAL_ZONE_SHIFT    
      USE SELWC, only: zshift, zmax, wzsconv
      IMPLICIT NONE
  
      INTEGER, PARAMETER ::JMAX=40    
      REAL :: X1, X2, FUNCVAL1, FUNCVAL2, XACC, FMID, FUNC1, RTBIS, DX, XMID
      INTEGER :: JJ,J  
               
      wzsconv=.false.
      
      ! FIRST, BRACKETING ROOT
      X2=zmax
      X1=0.001
          
      CALL  WD_ZONE_SHIFT_EQN(x1,FUNCVAL1)
    
      CALL  WD_ZONE_SHIFT_EQN(x2,FUNCVAL2)
                  
      DO JJ=1,JMAX
        IF(FUNCVAL1*FUNCVAL2 > 0.0)THEN
          IF(ABS(FUNCVAL1)  < ABS(FUNCVAL2))THEN
            X1=X1/2.0
            CALL  WD_ZONE_SHIFT_EQN(x1,FUNCVAL1)            
          ELSE
            X2=X2+1.5*(X2-X1)
            CALL  WD_ZONE_SHIFT_EQN(x2,FUNCVAL2)
          END IF                                          
        ELSE
          EXIT
        END IF      
      END DO                 
      
      ! FINDING ROOT BY BISECTION      
      XACC=0.01
      CALL  WD_ZONE_SHIFT_EQN(x2,FMID)
      CALL  WD_ZONE_SHIFT_EQN(x1,FUNC1)      
  !    IF(FUNC1*FMID.GE.0.) PAUSE 'ROOT MUST BE BRACKETED IN RTBIS'
      IF(FUNC1.LT.0.)THEN
        RTBIS=X1
        DX=X2-X1
      ELSE
        RTBIS=X2
        DX=X1-X2
      ENDIF
      DO J=1,JMAX
        DX=DX*.5
        XMID=RTBIS+DX
        CALL  WD_ZONE_SHIFT_EQN(xmid,FMID)        
        IF(FMID.LE.0.)RTBIS=XMID
        IF(ABS(DX).LT.XACC .OR. FMID.EQ.0.)THEN          
            zshift=rtbis
            wzsconv=.true.
          RETURN
        END IF
      END DO
      return

end subroutine WITHDRAWAL_ZONE_SHIFT
      
!***********************************************************************************************************************************
!**           W I T H D R A W A L    Z O N E    S H I F T    E Q U A T I O N                                                      **
!***********************************************************************************************************************************
subroutine WD_ZONE_SHIFT_EQN(zv,FUNCVALUE2)      
    USE GLOBAL, only: rho, g, nonzero, ID
    USE SELWC, only: vozavg, hovlap, hport, kstrshft, kwzlim
    IMPLICIT NONE
    
    REAL :: zv,funcvalue2
      
    FUNCVALUE2= vozavg - 0.7 * (hovlap/hport)**1.25 * SQRT((abs(RHO(KSTRshft,ID)-rho(kwzlim,id)))/(RHO(KSTRshft,ID)+NONZERO)*G*zv)  
      
    RETURN
end subroutine WD_ZONE_SHIFT_EQN

!***********************************************************************************************************************************
!**         S I N G L E    W E T    W E L L                                                                                       **
!***********************************************************************************************************************************
subroutine single_wet_well      
      USE GLOBAL; USE GEOMC; USE SELWC
      IMPLICIT NONE
      
      real    :: ejs
      integer :: jord
  
      INTEGER, PARAMETER ::JMAX=40    
      REAL :: X1, X2, FUNCVAL1, FUNCVAL2, XACC, FMID, FUNC1, RTBIS, DX, XMID
      INTEGER :: JJ,J, II
      
        QSUMport=0.0
        DO jsn=1,NSTR(JB)       ! summing specified port flows to get total flow through single well structure
          QSUMport=QSUMport+QSTR(jsn,JB)
        end do
        
        do jsn=1,nstr(jb)  ! sorting ports from highest elevation to lowest, adapted from Numerical Recipe's straight insertion subroutine
          jsorder(jsn)=jsn
          estrorder(jsn)=estr(jsn,jb)
        end do
        do jsn=2,nstr(jb)   ! pick out each element in turn
          ejs=estrorder(jsn)
          jord=jsorder(jsn)
          
          do ii=jsn-1,1,-1    ! look for a place to insert it
            if(estrorder(ii) >= ejs)go to 9112
            jsorder(ii+1)=jsorder(ii)
            estrorder(ii+1)=estrorder(ii)
          end do
          ii=0
9112      jsorder(ii+1)=jord      ! insert it
          estrorder(ii+1)=ejs
          
        end do
        
        IF (QSUMport > 0.0) THEN
          ! FIRST, BRACKETING ROOT
          !X2=wsel-estrorder(1)
          X2=wsels-el(kb(id)+1,id)
          X1=0.0
                
          CALL  Q_Single_Well_EQN(x1,FUNCVAL1)               
          CALL  Q_Single_Well_EQN(x2,FUNCVAL2)
                  
          DO JJ=1,JMAX
            IF(FUNCVAL1*FUNCVAL2 > 0.0)THEN
              IF(ABS(FUNCVAL1)  < ABS(FUNCVAL2))THEN
                X1=X1/2.0
                CALL  Q_Single_Well_EQN(x1,FUNCVAL1)            
              ELSE
                X2=X2+1.5*(X2-X1)
                CALL  Q_Single_Well_EQN(x2,FUNCVAL2)
              END IF                                          
            ELSE
              EXIT
            END IF      
          END DO                 
      
          ! FINDING ROOT BY BISECTION      
          XACC=0.00001
          CALL  Q_Single_Well_EQN(x2,FMID)
          CALL  Q_Single_Well_EQN(x1,FUNC1)      
      !    IF(FUNC1*FMID.GE.0.) PAUSE 'ROOT MUST BE BRACKETED IN RTBIS'
          IF(FUNC1*FMID.GE.0.)then
              continue
          end if
          IF(FUNC1.LT.0.)THEN
            RTBIS=X1
            DX=X2-X1
          ELSE
            RTBIS=X2
            DX=X1-X2
          ENDIF
          DO J=1,JMAX
            DX=DX*.5
            XMID=RTBIS+DX
            CALL  Q_Single_Well_EQN(xmid,FMID)        
            IF(FMID.LE.0.)RTBIS=XMID
            IF(ABS(DX).LT.XACC .OR. FMID.EQ.0.)THEN                         
              RETURN
            END IF
          END DO
        end if
              
      return
end subroutine single_wet_well
      
!***********************************************************************************************************************************
!**           S I N G L E    WET W E L L    F L O W    E Q U A T I O N                                                                **
!***********************************************************************************************************************************
subroutine Q_Single_Well_EQN(zxx,funcvalue2)     
    USE GLOBAL; USE GEOMC; USE SELWC
    use MAIN, only:aport, ploss,port_shape_type,aport_coef1,aport_coef2
    IMPLICIT NONE
    
    real :: rhodzsum, flowtot
    integer ::  jup, jord, k
    real :: zxx, funcvalue2
       
     hloss(1)=zxx
      
     do jsn=2,nstr(jb)
        jord=jsorder(jsn)
        jup=jsorder(jsn-1)
        rhodzsum=0.0
        do k=kstrs(jup),kstrs(jord)-1
          rhodzsum=rhodzsum+ ((RHO(k,ID)-RHO(KSTRs(jup),ID))*h(k,jw) + (RHO(k+1,ID)-RHO(KSTRs(jup),ID))*h(k+1,jw))/2.0
        end do
        
        hloss(jsn)=hloss(jsn-1)*RHO(KSTRs(jup),ID)/RHO(KSTRs(jord),ID)  + rhodzsum/RHO(KSTRs(jord),ID)
     end do
     
     flowtot=0.0
     do jsn=1,nstr(jb)
       jord=jsorder(jsn)
       QSTR(jsn,JB)=sqrt(2.0 * G * aport(jord,jb)**2 *hloss(jord)/ploss(jord,jb))
       flowtot=flowtot+QSTR(jsn,JB)
     end do
     funcvalue2=qsumport-flowtot
      
   RETURN
end subroutine Q_Single_Well_EQN

  
MODULE SELECTIVE1 
 REAL                                          :: NXTSTR, NXTTCD, NXTSPLIT,TCDFREQ,TFRQTMP
  CHARACTER(8)                                 :: TEMPC,TSPLTC
  CHARACTER(8), ALLOCATABLE, DIMENSION(:)      :: TCELEVCON,TCYEARLY,TCNTR,TSPLTCNTR,TSYEARLY,DYNSEL,ELCONTSPL,DYNSELSPLT
  INTEGER                                      :: NUMTEMPC,NUMTSPLT, TEMPN        
  INTEGER, ALLOCATABLE, DIMENSION(:)           :: TCNELEV,TCJB,TCJS,TCISEG,TSPLTJB,NOUTS,KSTRSPLT, JBMON, JSMON, NCOUNTCW,SELD
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: TCELEV, TEMPCRIT,QSTRFRAC
  REAL,          ALLOCATABLE, DIMENSION(:)     :: TCTEMP,TCTEND,TCTSRT,TCKLAY,TSPLTT,VOLM,QWDFRAC,TSTEND,TSTSRT,NXSEL,TEMP2,NXSELSPLT,TEMP3
  INTEGER, ALLOCATABLE, DIMENSION(:,:)         :: JSTSPLT, NCOUNTC, JSTSPLTT
  REAL,          ALLOCATABLE, DIMENSION(:,:)  :: VOLMC 
  LOGICAL, ALLOCATABLE, DIMENSION(:)          :: DYNSF,DYNSPF
  REAL, ALLOCATABLE, DIMENSION(:)             :: MINWL   ! Minimum water level above centerline of outlet if TCELEVCON is ON
END MODULE SELECTIVE1

SUBROUTINE SELECTIVEINIT

USE SELECTIVE1;   USE MAIN
  USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART

  IMPLICIT NONE
  
  INTEGER N, IFILE
  REAL DAYTEST
  CHARACTER(1) :: INFORMAT,CHAR1
  CHARACTER(8) :: CHAR8
  CHARACTER(30):: CHAR30

!**                                                   Task 2: Calculations                                                        **
!***********************************************************************************************************************************
      
      IFILE=1949
      TAVG=0.0
      TAVGW=0.0
      DO JB=1,NBR
        IF(NSTR(JB) > 0)THEN
        IFILE=IFILE+1
        WRITE (SEGNUM,'(I0)') JB
        SEGNUM = ADJUSTL(SEGNUM)
        L      = LEN_TRIM(SEGNUM)

        IF(RESTART_IN)THEN
                OPEN  (IFILE,FILE='str_br'//segnum(1:l)//'.csv',POSITION='APPEND')
                JDAY1=0.0
                REWIND (IFILE)
                READ   (IFILE,'(//)',END=13)
                DO WHILE (JDAY1 < JDAY)
                READ (IFILE,*,END=13) JDAY1            !READ (IFILE,'(F10.0)',END=13) JDAY1
                END DO
                BACKSPACE (IFILE)
                13     JDAY1=0.0    
        ELSE
                OPEN  (IFILE,FILE='str_br'//segnum(1:l)//'.csv',status='unknown')
                WRITE(IFILE,*)'Branch:,',jb,', # of structures:,',nstr(jb),', outlet temperatures'
                WRITE(IFILE,'("      JDAY,",<nstr(jb)>(6x,"T(C),"),<nstr(jb)>(3x,"Q(m3/s),"),<nstr(jb)>(4x,"ELEVCL,"))')
        ENDIF
        ENDIF
      END DO
 
      IF(NWD > 0)THEN
       IFILE=IFILE+1  
       IF(RESTART_IN)THEN
                OPEN  (IFILE,FILE='wd_out.opt',POSITION='APPEND')
                JDAY1=0.0
                REWIND (IFILE)
                READ   (IFILE,'(/)',END=14)
                DO WHILE (JDAY1 < JDAY)
                READ (IFILE,'(F10.0)',END=14) JDAY1
                END DO
                BACKSPACE (IFILE)
                14    JDAY1=0.0   
       ELSE
        OPEN  (IFILE,FILE='wd_out.opt',STATUS='unknown')
        WRITE(IFILE,*)'Withdrawals: # of withdrawals:',nwd,' outlet temperatures'
        WRITE(IFILE,'("      JDAY",<nwd>(6x,"T(C)"),<nwd>(3x,"Q(m3/s)"),<nwd>(4x,"ELEVCL"))')
       ENDIF
      end if
 
      
      OPEN(NUNIT,FILE='w2_selective.npt',STATUS='old')
      
      read(NUNIT,'(a)')CHAR1
      read(NUNIT,*)
      read(NUNIT,*)
      IF(CHAR1=='$')GO TO 7770
      ! check for commas another sign that it is comma delimited
      read(NUNIT,'(A)')CHAR30
     DO J=1,30
         IF(CHAR30(J:J)==',')THEN
             CHAR1='$'
             EXIT
         ENDIF
     ENDDO
     BACKSPACE(NUNIT)

7770  CONTINUE
      IF(CHAR1=='$')THEN
      READ(NUNIT,*)CHAR8,TFRQTMP
      READ(NUNIT,*)
      READ(NUNIT,*)
      READ(NUNIT,*)CHAR8,TEMPC,NUMTEMPC,TCDFREQ; TEMPC=ADJUSTR(TEMPC)
      ELSE
      READ(NUNIT,'(8X,F8.0)')TFRQTMP
      READ(NUNIT,'(//8X,A8,I8,F8.0)')TEMPC,NUMTEMPC,TCDFREQ
      ENDIF
      !NXTSTR=TMSTRT
      !NXTTCD=TMSTRT
      !NXTSPLIT=TMSTRT
      IF(JDAY>TMSTRT)THEN   ! SW 8/10/2023  During restart this allows continuing output
          NXTSTR=JDAY          
          NXTTCD=JDAY   
          NXTSPLIT=JDAY
      ELSE
          NXTSTR=TMSTRT
          NXTTCD=TMSTRT   
          NXTSPLIT=TMSTRT
      ENDIF
      
  ALLOCATE (TCNELEV(NUMTEMPC),TCJB(NUMTEMPC),TCJS(NUMTEMPC), TCELEV(NUMTEMPC,100),TCTEMP(NUMTEMPC),TCTEND(NUMTEMPC),TCTSRT(NUMTEMPC),NCOUNTC(NST,NBR),TCISEG(NUMTEMPC),TCKLAY(NUMTEMPC),TCELEVCON(NUMTEMPC)) 
  ALLOCATE (TCYEARLY(NUMTEMPC), JBMON(NUMTEMPC),JSMON(NUMTEMPC),TCNTR(NUMTEMPC)) 
  ALLOCATE (VOLM(NWB),NCOUNTCW(NWD),QWDFRAC(NWD),QSTRFRAC(NST,NBR),DYNSEL(NUMTEMPC),SELD(NUMTEMPC),NXSEL(NUMTEMPC),TEMP2(NUMTEMPC))       
  ALLOCATE (DYNSF(NUMTEMPC),MINWL(NUMTEMPC))    
  DYNSF=.FALSE.;MINWL=0.0
  
      DO J=1,2
      READ(NUNIT,*)
      END DO
      NCOUNTC=0
      DO J=1,NUMTEMPC
            IF(CHAR1=='$')THEN
              READ(NUNIT,*)CHAR8,TCNTR(J),TCJB(J),TCJS(J),TCYEARLY(J),TCTSRT(J),TCTEND(J),TCTEMP(J),TCNELEV(J),(TCELEV(J,N),N=1,TCNELEV(J))
              TCNTR(J)=ADJUSTR(TCNTR(J))
              TCYEARLY=ADJUSTR(TCYEARLY) 
            ELSE
              READ(NUNIT,'(8X,A8,I8,I8,A8,F8.0,F8.0,F8.0,I8,10(F8.0))')TCNTR(J),TCJB(J),TCJS(J),TCYEARLY(J),TCTSRT(J),TCTEND(J),TCTEMP(J),TCNELEV(J),(TCELEV(J,N),N=1,TCNELEV(J))
            ENDIF
            
        IF(TCNTR(J)=='      ST')THEN      
        TCELEV(J,TCNELEV(J)+1)=ESTR(TCJS(J),TCJB(J))   ! ALWAYS PUT THE ORIGINAL ELEVATION AS THE LAST ELEVATION
        ELSE
        TCELEV(J,TCNELEV(J)+1)=EWD(TCJS(J))   ! ALWAYS PUT THE ORIGINAL ELEVATION AS THE LAST ELEVATION
        ENDIF
      END DO
      DO J=1,2
      READ(NUNIT,*)
      END DO
      DO J=1,NUMTEMPC
                      IF(CHAR1=='$')THEN
                                READ(NUNIT,*)CHAR8,TCISEG(J),TCKLAY(J),DYNSEL(J);DYNSEL(J)=ADJUSTR(DYNSEL(J)) 
                      ELSE 
                                READ(NUNIT,'(8X,I8,F8.0,A8)')TCISEG(J),TCKLAY(J),DYNSEL(J) 
                      ENDIF
                      
      END DO
      DO J=1,2
      READ(NUNIT,*)
      END DO
      DO J=1,NUMTEMPC
                      IF(CHAR1=='$')THEN
                                      READ(NUNIT,*)CHAR8,TCELEVCON(J),MINWL(J);tcelevcon(J)=ADJUSTR(tcelevcon(J)) 
                      ELSE 
                                      READ(NUNIT,'(8X,A8,F8.0)')TCELEVCON(J),MINWL(J) 
                      ENDIF
      END DO
      DO J=1,2
      READ(NUNIT,*)
      END DO
                      IF(CHAR1=='$')THEN
                                    READ(NUNIT,*)CHAR8,TSPLTC,NUMTSPLT;TSPLTC=ADJUSTR(TSPLTC) 
                      ELSE 
                                    READ(NUNIT,'(8X,A8,I8)')TSPLTC,NUMTSPLT
                      ENDIF

      
      ALLOCATE(TSYEARLY(NUMTSPLT),TSTSRT(NUMTSPLT),TSTEND(NUMTSPLT),TSPLTJB(NUMTSPLT),TSPLTT(NUMTSPLT),NOUTS(NUMTSPLT),   &
               JSTSPLT(NUMTSPLT,10),KSTRSPLT(NUMTSPLT),TSPLTCNTR(NUMTSPLT))
      ALLOCATE(JSTSPLTT(NUMTSPLT,10),ELCONTSPL(NUMTSPLT),DYNSELSPLT(NUMTSPLT),DYNSPF(NUMTSPLT),NXSELSPLT(NUMTSPLT),TEMP3(NUMTSPLT))
      DYNSPF=.FALSE.
      
      DO J=1,2
      READ(NUNIT,*)
      END DO
      DO J=1,NUMTSPLT
              IF(CHAR1=='$')THEN
                    READ(NUNIT,*)CHAR8,TSPLTCNTR(J),TSPLTJB(J),TSYEARLY(J),TSTSRT(J),TSTEND(J),TSPLTT(J),NOUTS(J),(JSTSPLTT(J,N),N=1,2),ELCONTSPL(J),DYNSELSPLT(J)
                    TSPLTCNTR(J)=ADJUSTR(TSPLTCNTR(J))
                    TSYEARLY(J)=ADJUSTR(TSYEARLY(J))
                    ELCONTSPL(J)=ADJUSTR(ELCONTSPL(J))
                    DYNSELSPLT(J)=ADJUSTR(DYNSELSPLT(J))
              ELSE 
                    READ(NUNIT,'(8X,A8,I8,A8,F8.0,F8.0,F8.0,I8,2I8,A8,A8)')TSPLTCNTR(J),TSPLTJB(J),TSYEARLY(J),TSTSRT(J),TSTEND(J),TSPLTT(J),NOUTS(J),(JSTSPLTT(J,N),N=1,2),ELCONTSPL(J),DYNSELSPLT(J)
              ENDIF
      NOUTS(J)=2                ! NUMBER OF OUTLETS FOR EACH SPLIT FLOW PERIOD LIMITED TO 2
      !IF(NOUTS(J).GT.2)WRITE(*,*)'TCD NOUTS > 2 - ONLY FIRST 2 WILL BE USED'
      ENDDO
      JSTSPLT=JSTSPLTT                                                                                             ! CB 10/14/11 START
      DO J=1,NUMTSPLT  !REODERING OUTLETS SO THAT HIGHEST ELEVATION STRUCTURE ON TOP (ASSUMING 2 SPLIT OUTLETS) 
!        IF(TCNTR(J) == '      ST')THEN
        IF(TSPLTCNTR(J) == '      ST')THEN                                                                        ! cb 11/11/12
          IF(ESTR(JSTSPLTT(J,1),TSPLTJB(J)) < ESTR(JSTSPLTT(J,2),TSPLTJB(J)))THEN                               
            JSTSPLT(J,1)=JSTSPLTT(J,2)                                                                          
            JSTSPLT(J,2)=JSTSPLTT(J,1)                                                                          
          END IF                                                                                                
!        ELSE IF(TCNTR(J) == '      WD')THEN
        ELSE IF(TSPLTCNTR(J) == '      WD')THEN                                                                        ! cb 11/11/12
          IF(EWD(JSTSPLTT(J,1)) < EWD(JSTSPLTT(J,2)))THEN                                    
            JSTSPLT(J,1)=JSTSPLTT(J,2)                                                                          
            JSTSPLT(J,2)=JSTSPLTT(J,1)                                                                          
          END IF                                                                                                
        END IF
      END DO                                                                                                       ! CB 10/14/11 END
      DO J=1,2
      READ(NUNIT,*)
      END DO
               IF(CHAR1=='$')THEN
                    READ(NUNIT,*)CHAR8,TEMPN
              ELSE 
                    READ(NUNIT,'(8X,I8)')TEMPN
              ENDIF
      DO J=1,2
      READ(NUNIT,*)
      END DO
      ALLOCATE(TEMPCRIT(NWB,TEMPN),VOLMC(NWB,TEMPN))
      DO J=1,TEMPN
        IF(CHAR1=='$')THEN
        READ(NUNIT,*)CHAR8,(TEMPCRIT(JW,J),JW=1,NWB)   ! NOTE MAX OF 100 WATERBODIES   sw 4/20/15  
        ELSE
        READ(NUNIT,'(8X,100F8.0)')(TEMPCRIT(JW,J),JW=1,NWB)   ! NOTE MAX OF 100 WATERBODIES   sw 4/20/15
        ENDIF
      END DO
      CLOSE(NUNIT)

      
      DO JW=1,NWB
        IFILE=IFILE+1
        WRITE (SEGNUM,'(I0)') JW
        SEGNUM = ADJUSTL(SEGNUM)
        L      = LEN_TRIM(SEGNUM)
         IF(RESTART_IN)THEN
                OPEN  (IFILE,FILE='VOLUME_WB'//SEGNUM(1:L)//'.OPT',POSITION='APPEND')
                JDAY1=0.0
                REWIND (IFILE)
                READ   (IFILE,'(/)',END=15)
                DO WHILE (JDAY1 < JDAY)
                READ (IFILE,'(F10.0)',END=15) JDAY1
                END DO
                BACKSPACE (IFILE)
                15    JDAY1=0.0   
       ELSE
        OPEN  (IFILE,FILE='VOLUME_WB'//SEGNUM(1:L)//'.OPT',STATUS='UNKNOWN')
        WRITE(IFILE,4315)
       ENDIF
      ENDDO
      
4315  FORMAT("JDAY    VOLUME    ",<TEMPN>("VOLCRIT      "))


! INITIALIZING STRUCTURE ELEVATION IF STRUCTURE
IF(TEMPC=='      ON')THEN     
  DO JW=1,NWB
   DO JB=BS(JW),BE(JW)
    DO JS=1,NST
     DO J=1,NUMTEMPC        
       IF(TCJB(J) == JB .AND. TCJS(J) == JS .AND. TCNTR(J) == '      ST')THEN
           IF(TCYEARLY(J) == '     OFF')THEN
             DAYTEST=JDAY
           ELSE
             DAYTEST=REAL(JDAYG)+JDAY-INT(JDAY)
           END IF
           IF(DAYTEST >= TCTSRT(J) .AND. DAYTEST < TCTEND(J))THEN               
               ! MAKING SURE THAT STRUCTURE IS BELOW WATER SURFACE
             DO NN=1,TCNELEV(J)
               IF(TCELEV(J,NN) < ELWS(DS(JB)))THEN
                 NCOUNTC(JS,JB)=NN
                 ESTR(JS,JB)=TCELEV(J,NCOUNTC(JS,JB))
                 EXIT
               END IF                 
             END DO
		   END IF
	   END IF
	 END DO
	END DO
   END DO
  END DO
   
   
   ! INITIALIZING STRUCTURE ELEVATION IF WITHDRAWAL

  DO JWD=1,NWD
   
     DO J=1,NUMTEMPC        
       IF(TCJS(J) == JWD .AND. TCNTR(J) == '      WD')THEN
           IF(TCYEARLY(J) == '     OFF')THEN
             DAYTEST=JDAY
           ELSE
             DAYTEST=REAL(JDAYG)+JDAY-INT(JDAY)
           END IF
           IF(DAYTEST >= TCTSRT(J) .AND. DAYTEST < TCTEND(J))THEN               
               ! MAKING SURE THAT STRUCTURE IS BELOW WATER SURFACE
             DO NN=1,TCNELEV(J)
               IF(TCELEV(J,NN) < ELWS(IWD(JWD)))THEN
                 NCOUNTCW(JWD)=NN
                 EWD(JWD)=TCELEV(J,NCOUNTCW(JWD))
                 EXIT
               END IF                 
             END DO
		   END IF
	   END IF
	 END DO

  END DO

  
  ! OPEN DYNAMIC SELECTIVE WITHDRAWAL FILES
  
  DO J=1,numtempc
     if(DYNSEL(J) == '      ON')then
     WRITE (SEGNUM,'(I0)') J     
     SEGNUM = ADJUSTL(SEGNUM)
     L      = LEN_TRIM(SEGNUM)   
     SELD(J) = 1009+J  
     OPEN (SELD(J),FILE='dynselective'//SEGNUM(1:L)//'.npt',STATUS='OLD') 
     
      READ(SELD(J),'(A1)')INFORMAT    ! SW 8/28/2019
      IF(INFORMAT=='$')DYNSF(J)=.TRUE.
    
        IF(DYNSF(J))THEN
        READ (SELD(J),'(/)')
        READ (SELD(J),*) NXSEL(J),TEMP2(J)
        tctemp(J)=TEMP2(J)
        READ (SELD(J),*) NXSEL(J),TEMP2(J)       
        
        ELSE
        READ (SELD(J),'(//1000F8.0)') NXSEL(J),TEMP2(J)
        tctemp(J)=TEMP2(J)
        READ (SELD(J),'(1000F8.0)') NXSEL(J),TEMP2(J)
        ENDIF
     
      !READ (SELD(J),'(///1000F8.0)') NXSEL(J),TEMP2(J)
      !  tctemp(J)=TEMP2(J)
      !READ (SELD(J),'(1000F8.0)') NXSEL(J),TEMP2(J)
     END IF
  ENDDO 
END IF

IF(TSPLTC == '      ON')THEN
    DO J=1,NUMTSPLT
     IF(DYNSELSPLT(J) == '      ON')then
     WRITE (SEGNUM,'(I0)') J     
     SEGNUM = ADJUSTL(SEGNUM)
     L      = LEN_TRIM(SEGNUM)   
     SELD(J) = 1059+J  
     OPEN (SELD(J),FILE='dynselective'//SEGNUM(1:L)//'_splt.npt',STATUS='OLD') 
     
      READ(SELD(J),'(A1)')INFORMAT    ! SW 8/28/2019
      IF(INFORMAT=='$')DYNSPF(J)=.TRUE.
    
        IF(DYNSPF(J))THEN
        READ (SELD(J),'(/)')
        READ (SELD(J),*) NXSELSPLT(J),TEMP3(J)
        TSPLTT(J)=TEMP3(J)
        READ (SELD(J),*) NXSELSPLT(J),TEMP3(J)       
        
        ELSE
        READ (SELD(J),'(//1000F8.0)') NXSELSPLT(J),TEMP3(J)
        TSPLTT(J)=TEMP2(J)
        READ (SELD(J),'(1000F8.0)') NXSELSPLT(J),TEMP3(J)
        ENDIF
     
     END IF
    ENDDO 
    OPEN(2900,FILE='Split_Temp_Debug.csv',status='unknown')
        WRITE(2900,*)'JDAY,TSPLTT,TTOP,TBOT,QALL,QSTR1,QSTR2,ESTR1,ESTR2'
ENDIF

  
 RETURN

END SUBROUTINE SELECTIVEINIT

!***********************************************************************************************************************************
!**                                                   S E L E C T I V E                                                           **
!***********************************************************************************************************************************
SUBROUTINE SELECTIVE
 USE SELECTIVE1
  USE MAIN
  USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  
  IMPLICIT NONE
  !** Timestep violation entry point  210 CONTINUE                
  INTEGER JJ, JJW, KK, KS, IFILE, KSTR
  REAL DAYTEST, ELR, QALL, TCOMP, TEMPBOT, TEMPEST, TEMPTOP, TMOD, WSEL

  IF(TSPLTC=='      ON')THEN    
    DO J=1,NUMTSPLT
      IF(TSYEARLY(J) == '     OFF')THEN
        DAYTEST=JDAY
      ELSE
        DAYTEST=REAL(JDAYG)+JDAY-INT(JDAY)
      END IF
      IF(NXTSPLIT > TSTSRT(J) .AND. DAYTEST <= TSTSRT(J))THEN
        NXTSPLIT=TSTSRT(J)
      END IF


  IF(DYNSELSPLT(J) == '      ON')THEN
    SELD(J)=1059+J
     DO WHILE (JDAY >= NXSELSPLT(J))
        TSPLTT(J)=TEMP3(J)
                IF(DYNSPF(J))THEN
                    READ (SELD(J),*) NXSELSPLT(J),TEMP3(J)
                ELSE
                    READ (SELD(J),'(1000F8.0)') NXSELSPLT(J),TEMP3(J)
                ENDIF   
    END DO
   ENDIF
   ENDDO
END IF

 IF(TSPLTC=='      ON'.AND.JDAY.GE.NXTSPLIT)THEN  
 
  DO J=1,NUMTSPLT
        IF(TSYEARLY(J) == '     OFF')THEN
            DAYTEST=JDAY
          ELSE
        DAYTEST=REAL(JDAYG)+JDAY-INT(JDAY)
        END IF
   IF(DAYTEST >= TSTSRT(J) .AND. DAYTEST < TSTEND(J))THEN 
    ! DO STRUCTURES FIRST
    DO JW=1,NWB
        DO JB=BS(JW),BE(JW)
            IF(TSPLTJB(J) == JB .AND. TSPLTCNTR(J) == '      ST')THEN
                QALL=0.0
                DO JJ=1,NOUTS(J)
                QALL=QALL+QSTR(JSTSPLT(J,JJ),TSPLTJB(J))   ! SUM UP ALL THE FLOWS
                ELR  = SINA(JB)*DLX(DS(JB))*0.5
                    DO K=KTWB(JW),KB(DS(JB))
                    IF (EL(K,DS(JB))-ELR < ESTR(JSTSPLT(J,JJ),TSPLTJB(J))) EXIT                                                                               !SW 10/17/01
                    END DO
                KSTR = K-1
                KSTRSPLT(JJ) = MIN(KSTR,KB(DS(JB)))
                ENDDO               
              DO JJ=1,NOUTS(J)               ! cb 11/11/12 dividing total flow between outlets for temperature test - if no flow there is no temperature test
                  QSTR(JSTSPLT(J,JJ),TSPLTJB(J)) = qall/real(nouts(j))
              ENDDO               
              ID=DS(JB)
              ELR  = SINA(JB)*DLX(ID)*0.5          ! CB 10/14/11
              WSEL = ELWS(ID)-ELR                  ! CB 10/14/11
              kt=ktwb(jw)      ! cb 07/24/19
              CALL DOWNSTREAM_WITHDRAWAL_ESTIMATE(JSTSPLT(J,1),TEMPTOP,ESTR(JSTSPLT(J,1),TSPLTJB(J)))
              CALL DOWNSTREAM_WITHDRAWAL_ESTIMATE(JSTSPLT(J,2),TEMPBOT,ESTR(JSTSPLT(J,2),TSPLTJB(J)))
             IF(ESTR(JSTSPLT(J,1),TSPLTJB(J)) > WSEL .AND. ELCONTSPL(J) =='     OFF') THEN   ! NO FLOWS THROUG THIS OUTLET IF WSEL BELOW LEVEL OF OUTLET  ! CB 10/14/11
               QSTR(JSTSPLT(J,1),TSPLTJB(J))=0.0
               QSTRFRAC(JSTSPLT(J,1),TSPLTJB(J))=0.0
              
            ELSE IF(TEMPTOP > TSPLTT(J)  .AND.  TEMPBOT > TSPLTT(J) ) THEN   ! NO FLOWS THROUG THIS OUTLET IF T1 AND T2 > TCRITERIA
               QSTR(JSTSPLT(J,1),TSPLTJB(J))=0.0
               QSTRFRAC(JSTSPLT(J,1),TSPLTJB(J))=0.0

              !ELSEIF(T2(KSTRSPLT(1),DS(JB)) < TSPLTT(J)) THEN   ! ALL FLOWS FROM TOP IF TCRITERIA < TOUTLET
              ELSEIF(TEMPTOP < TSPLTT(J)) THEN   ! ALL FLOWS FROM TOP IF TCRITERIA < TOUTLET
               QSTR(JSTSPLT(J,1),TSPLTJB(J))=QALL
               QSTRFRAC(JSTSPLT(J,1),TSPLTJB(J))=1.0

              ELSE
                !QSTR(JSTSPLT(J,1),TSPLTJB(J))=QALL*(TSPLTT(J)-T2(KSTRSPLT(2),DS(JB)))/(T2(KSTRSPLT(1),DS(JB))-T2(KSTRSPLT(2),DS(JB)))
                IF(ABS(TEMPTOP-TEMPBOT) < 0.0001)THEN
                  QSTR(JSTSPLT(J,1),TSPLTJB(J))=QALL
                  QSTRFRAC(JSTSPLT(J,1),TSPLTJB(J))=1.0
                ELSE
                  QSTR(JSTSPLT(J,1),TSPLTJB(J))=QALL*(TSPLTT(J)-TEMPBOT)/(TEMPTOP-TEMPBOT)
                  QSTRFRAC(JSTSPLT(J,1),TSPLTJB(J))=QSTR(JSTSPLT(J,1),TSPLTJB(J))/QALL
                END IF
              ENDIF

              QSTR(JSTSPLT(J,2),TSPLTJB(J))=QALL-QSTR(JSTSPLT(J,1),TSPLTJB(J))
              QSTRFRAC(JSTSPLT(J,2),TSPLTJB(J))=QSTR(JSTSPLT(J,2),TSPLTJB(J))/QALL
              write(2900,'(9(f12.4,","))')jday,tspltt(j),TEMPTOP,TEMPBOT,QALL,QSTR(JSTSPLT(J,1),TSPLTJB(J)),QSTR(JSTSPLT(J,2),TSPLTJB(J)),ESTR(JSTSPLT(J,1),TSPLTJB(J)),ESTR(JSTSPLT(J,2),TSPLTJB(J))
             EXIT
             END IF
         END DO
       END DO
    ! DO WITHDRAWALS NEXT
      DO JWD=1,NWD
            IF(TSPLTCNTR(J) == '      WD')THEN
                QALL=0.0
               DO JJB=1,NBR
                 IF (IWD(JWD) >= US(JJB) .AND. IWD(JWD) <= DS(JJB)) EXIT
               END DO
               DO JJW=1,NWB
                 IF (JJB >= BS(JJW) .AND. JJB <= BE(JJW)) EXIT
               END DO
               DO JJ=1,NOUTS(J)
                QALL=QALL+QWD(JSTSPLT(J,JJ))   ! SUM UP ALL THE FLOWS
                ELR  = SINA(JJB)*DLX(IWD(JWD))*0.5
                    DO K=KTWB(JJW),KB(IWD(JWD))
                    IF (EL(K,IWD(JWD))-ELR < EWD(JSTSPLT(J,JJ))) EXIT                                                                               !SW 10/17/01
                    END DO
                KSTR = K-1
                KSTRSPLT(JJ) = MIN(KSTR,KB(IWD(JWD)))
               ENDDO
               JJ=1               ! ASSIGN FLOW TO FIRST OUTLET               
               WSEL = ELWS(IWD(JWD))-ELR                  ! CB 10/14/11
               I=IWD(JWD)
               kt=ktwb(jjw)      ! cb 07/24/19
               CALL LATERAL_WITHDRAWAL_ESTIMATE(JSTSPLT(J,1),TEMPTOP,EWD(JSTSPLT(J,1)))
               CALL LATERAL_WITHDRAWAL_ESTIMATE(JSTSPLT(J,2),TEMPBOT,EWD(JSTSPLT(J,2)))              
              IF(EWD(JSTSPLT(J,1)) > WSEL .AND. TCELEVCON(J) =='     OFF') THEN
                QWD(JSTSPLT(J,1))=0.0
                QWDFRAC(JSTSPLT(J,1))=0.0             
             ELSE IF(TEMPTOP > TSPLTT(J)  .AND.  TEMPBOT > TSPLTT(J) ) THEN   ! NO FLOWS THROUG THIS OUTLET IF T1 AND T2 > TCRITERIA
               QWD(JSTSPLT(J,1))=0.0
               QWDFRAC(JSTSPLT(J,1))=0.0

              ELSEIF(TEMPTOP < TSPLTT(J)) THEN   ! ALL FLOWS FROM TOP IF TCRITERIA < TOUTLET
               QWD(JSTSPLT(J,1))=QALL
               QWDFRAC(JSTSPLT(J,1))=1.0

              ELSE
                !QWD(JSTSPLT(J,1))=QALL*(TSPLTT(J)-T2(KSTRSPLT(2),IWD(JWD)))/(T2(KSTRSPLT(1),IWD(JWD))-T2(KSTRSPLT(2),IWD(JWD)))
                IF(ABS(TEMPTOP-TEMPBOT) < 0.0001)THEN
                  QWD(JSTSPLT(J,1))=QALL
                  QWDFRAC(JSTSPLT(J,1))=1.0
                ELSE
                  QWD(JSTSPLT(J,1))=QALL*(TSPLTT(J)-TEMPBOT)/(TEMPTOP-TEMPBOT)
                  QWDFRAC(JSTSPLT(J,1))=QWD(JSTSPLT(J,1))/QALL
                END IF
              ENDIF

              QWD(JSTSPLT(J,2))=QALL-QWD(JSTSPLT(J,1))
              QWDFRAC(JSTSPLT(J,2))=QWD(JSTSPLT(J,2))/QALL
             EXIT
             END IF
       END DO
     ENDIF
     ENDDO
  
   NXTSPLIT=NXTSPLIT+TCDFREQ
  END IF
  IF(TSPLTC=='      ON')THEN

    DO J=1,NUMTSPLT
    IF(TSYEARLY(J) == '     OFF')THEN
            DAYTEST=JDAY
          ELSE
        DAYTEST=REAL(JDAYG)+JDAY-INT(JDAY)
        END IF
    IF(DAYTEST >= TSTSRT(J) .AND. DAYTEST < TSTEND(J))THEN 
    ! DO STRUCTURES FIRST
      DO JW=1,NWB
        DO JB=BS(JW),BE(JW)
            IF(TSPLTJB(J) == JB .AND. TSPLTCNTR(J) == '      ST')THEN
                QALL=0.0
                DO JJ=1,NOUTS(J)
                QALL=QALL+QSTR(JSTSPLT(J,JJ),TSPLTJB(J))   ! SUM UP ALL THE FLOWS
                ELR  = SINA(JB)*DLX(DS(JB))*0.5
                    DO K=KTWB(JW),KB(DS(JB))
                    IF (EL(K,DS(JB))-ELR < ESTR(JSTSPLT(J,JJ),TSPLTJB(J))) EXIT                                                                               !SW 10/17/01
                    END DO
                KSTR = K-1
                KSTRSPLT(JJ) = MIN(KSTR,KB(DS(JB)))
                ENDDO               
              QSTR(JSTSPLT(J,1),TSPLTJB(J))=QSTRFRAC(JSTSPLT(J,1),TSPLTJB(J))*QALL
              QSTR(JSTSPLT(J,2),TSPLTJB(J))=QSTRFRAC(JSTSPLT(J,2),TSPLTJB(J))*QALL
             EXIT
             END IF
        END DO
      END DO
    ! DO WITHDRAWALS NEXT
      DO JWD=1,NWD
            IF(TSPLTCNTR(J) == '      WD')THEN
                QALL=0.0
                DO JJB=1,NBR
                  IF (IWD(JWD) >= US(JJB) .AND. IWD(JWD) <= DS(JJB)) EXIT
                END DO
                DO JJW=1,NWB
                  IF (JJB >= BS(JJW) .AND. JJB <= BE(JJW)) EXIT
                END DO
               DO JJ=1,NOUTS(J)
                QALL=QALL+QWD(JSTSPLT(J,JJ))   ! SUM UP ALL THE FLOWS
                ELR  = SINA(JJB)*DLX(IWD(JWD))*0.5
                    DO K=KTWB(JJW),KB(IWD(JWD))
                    IF (EL(K,IWD(JWD))-ELR < EWD(JSTSPLT(J,JJ))) EXIT                                                                               !SW 10/17/01
                    END DO
                KSTR = K-1
                KSTRSPLT(JJ) = MIN(KSTR,KB(IWD(JWD)))
               ENDDO               
               QWD(JSTSPLT(J,1))=  QWDFRAC(JSTSPLT(J,1))*QALL
               QWD(JSTSPLT(J,2))=  QWDFRAC(JSTSPLT(J,2))*QALL
             EXIT
             END IF
      END DO
      ENDIF
    ENDDO
   
  ENDIF   
 

      IF (JDAY.GE.NXTSTR) THEN
        NXTSTR = NXTSTR+TFRQTMP   
        IFILE=1949
        DO JB=1,NBR
            IF(NSTR(JB) > 0)THEN
            IFILE=IFILE+1
            WRITE (IFILE,'(F10.4,",",<NSTR(JB)>(F10.2,","),<NSTR(JB)>(F10.2,","),<NSTR(JB)>(F10.2,","))') JDAY,(TAVG(I,JB),I=1,NSTR(JB)),(QSTR(I,JB),I=1,NSTR(JB)),(ESTR(I,JB),I=1,NSTR(JB))
            END IF
         ENDDO          
          IF(NWD > 0)THEN
            IFILE=IFILE+1
            WRITE (IFILE,'(F10.4,<NWD>F10.2,<NWD>F10.2,<NWD>F10.2)') JDAY,(TAVGW(I),I=1,NWD),(QWD(I),I=1,NWD),(EWD(I),I=1,NWD)
          END IF
         ! TEMPERATURE CONTROL LOGIC 

         ! COMPUTING RESERVOIR VOLUME AND VOLUME BELOW 'TEMPCRIT'        
        VOLMC=0.0
        VOLM=0.0
        DO JW=1,NWB
         KT = KTWB(JW)
           DO JB=BS(JW),BE(JW)           
             DO I=CUS(JB),DS(JB)
               VOLM(JW) = VOLM(JW) +BH2(KT,I)*DLX(I)               
               DO K=KT+1,KB(I)
                 VOLM(JW) = VOLM(JW)+BH(K,I)*DLX(I)               
               END DO
               DO KK=1,TEMPN                                         
                 IF(T2(KT,I).LE.TEMPCRIT(JW,KK))VOLMC(JW,KK) = VOLMC(JW,KK)+BH2(KT,I)*DLX(I)                                                 
                 DO K=KT+1,KB(I)                 
                   IF(T2(K,I).LE.TEMPCRIT(JW,KK))VOLMC(JW,KK) = VOLMC(JW,KK)+BH(K,I)*DLX(I)
                 END DO
               END DO               
             END DO         
           END DO
     
         IFILE=IFILE+1
         WRITE(IFILE,5315)JDAY,VOLM(JW),(VOLMC(JW,KK), KK=1,TEMPN)
5315     FORMAT(F8.2,100(G12.4,G12.4))
       ENDDO
         
      ENDIF



IF(TEMPC=='      ON'.AND.JDAY.GE.NXTTCD)THEN  

! IF DYNAMIC SELECTIVE CHANGE TEMPERATURE 

 DO J=1,NUMTEMPC
  IF(DYNSEL(J) == '      ON')THEN
    SELD(J)=1009+J
     DO WHILE (JDAY >= NXSEL(J))
        TCTEMP(J)=TEMP2(J)
                IF(DYNSF(J))THEN
                    READ (SELD(J),*) NXSEL(J),TEMP2(J)
                ELSE
                    READ (SELD(J),'(1000F8.0)') NXSEL(J),TEMP2(J)
                ENDIF   
    END DO
   ENDIF
  ENDDO

  
! STRUCTURES  
  
  DO JW=1,NWB
   DO JB=BS(JW),BE(JW)
    DO JS=1,NST
     DO J=1,NUMTEMPC      

          
      IF(TCJB(J) == JB .AND. TCJS(J) == JS .AND.  TCNTR(J) == '      ST')THEN
          IF(TCISEG(J).EQ.0)THEN
            TCOMP=TAVG(TCJS(J),TCJB(J))   !CB 9/8/06   TAVG(JSMON(J),JBMON(J))
          ELSEIF(TCISEG(J) < 0)THEN
            TCOMP=TWDO(ABS(TCISEG(J)))      ! SW 11/26/10       
          ELSE

! CHECKING TO SEE IF THE MONITORING SEGMENT TCISEG IS IN THE SAME BRANCH AND WATER BODY AS THE STRUCTURE
            DO JJB=1,NBR
              IF (TCISEG(J) >= US(JJB) .AND. TCISEG(J) <= DS(JJB)) EXIT
            END DO
            DO JJW=1,NWB
              IF (JJB >= BS(JJW) .AND. JJB <= BE(JJW)) EXIT
            END DO

            IF (TCKLAY(J)< 0) THEN                                                                                       
              K = INT(ABS(TCKLAY(J)))
            ELSE
              DO K=KTWB(JJW),KB(TCISEG(J))
                IF (DEPTHB(K,TCISEG(J)) > TCKLAY(J)) EXIT                                                                      
              END DO
              K = MIN(K,KB(TCISEG(J)))                                                                                         
            END IF
            TCOMP=T2(K,TCISEG(J))
          ENDIF
          IF(TCYEARLY(J) == '     OFF')THEN
            DAYTEST=JDAY
          ELSE
            DAYTEST=REAL(JDAYG)+JDAY-INT(JDAY)
          END IF          
          IF(DAYTEST >= TCTSRT(J) .AND. DAYTEST < TCTEND(J))THEN
               ACTIVE_RULE_W2SELECTIVE(JS,JB)=.TRUE.
               IF(TCOMP > TCTEMP(J) .AND. TCNELEV(J) > NCOUNTC(JS,JB))THEN
               ! MAKING SURE THAT THE NEXT LOWER STRUCTURE FOR A PARTICULAR 'J' IS FOUND
                 DO NN=NCOUNTC(JS,JB)+1,TCNELEV(J)                 
                   IF(TCELEV(J,NN) < ESTR(JS,JB))THEN
                      NCOUNTC(JS,JB)=NN
                      ESTR(JS,JB)=TCELEV(J,NCOUNTC(JS,JB))
                      EXIT
                   END IF                 
                 END DO                                               
               ELSEIF(TCOMP < TCTEMP(J) .AND.  NCOUNTC(JS,JB).GT. 1)THEN
                 ! TO PREVENT THIS HAPPENING AT EACH TIME IT CHECKS IT AND HENCE OSCIALLTING BACK AND FORTH - CHECK THE TEMP AT THE UPPER OUTLET ALSO
                 IF(TCISEG(J) > 0)THEN  
                    IF(JB.EQ.JJB)THEN
                      DO KS=KTWB(JW),KB(DS(JB))
                        IF (DEPTHB(KS,TCISEG(J)) > TCELEV(J,NCOUNTC(JS,JB)-1)) EXIT                                                                          !TC 01/03/02
                      END DO
                      KS = MIN(KS,KB(TCISEG(J)))
                      TMOD= T2(KS,DS(JB))
                    ELSE
                      TMOD=T2(K,TCISEG(J))                      
                    END IF
                    IF(TMOD < TCTEMP(J) .AND. TCELEV(J,NCOUNTC(JS,JB)-1) < ELWS(DS(JB)))THEN                      
                      ! MAKING SURE THAT THE NEXT UPPER STRUCTURE FOR A PARTICULAR 'J' IS FOUND
                      DO NN=NCOUNTC(JS,JB)-1,1,-1
                        IF(TCELEV(J,NN) > ESTR(JS,JB))THEN
                          NCOUNTC(JS,JB)=NN
                          ESTR(JS,JB)=TCELEV(J,NCOUNTC(JS,JB))
                          EXIT
                        END IF                 
                      END DO
                    ENDIF
                 END IF  ! CB 9/8/06
                 IF(TCISEG(J).EQ.0)THEN
! CALCULATE THE ESTIMATED OUTFLOW TEMPERATURE AT HIGHER PORTS WHEN TCOMP<TCTEMP(J), AND MOVE UP IF HIGHER PORT STILL MEETS TO CRITERIA - THIS DOESN'T HAPPEN WHEN TCISEG < 0
                   DO NN=1,NCOUNTC(JS,JB)-1
                     ID=DS(JB)
                     KT=KTWB(JW)                     
                     CALL DOWNSTREAM_WITHDRAWAL_ESTIMATE(JS,TEMPEST,TCELEV(J,NN))
                     IF(TEMPEST < TCTEMP(J) .AND. TCELEV(J,NN) < ELWS(DS(JB)))THEN
                       NCOUNTC(JS,JB)=NN
                       ESTR(JS,JB)=TCELEV(J,NCOUNTC(JS,JB))
                       EXIT
                     END IF
                   END DO
                 END IF
               ENDIF
               IF(TCELEVCON(J) =='      ON' .AND. TCNELEV(J) > NCOUNTC(JS,JB).AND. ESTR(JS,JB) > (ELWS(DS(JB))-MINWL(J)))THEN  
                 NCOUNTC(JS,JB)=NCOUNTC(JS,JB)+1
                 ESTR(JS,JB)=TCELEV(J,NCOUNTC(JS,JB))
               END IF
          ELSE
              ACTIVE_RULE_W2SELECTIVE(JS,JB)=.FALSE.
          ENDIF          
      ENDIF            
     END DO
    END DO
   END DO
  ENDDO
  
 ! Withdrawals 
  

  DO JWD=1,NWD
     DO J=1,NUMTEMPC
      IF(TCJS(J) == JWD .AND.  TCNTR(J) == '      WD')THEN
          IF(TCISEG(J).EQ.0)THEN
!           TCOMP=TOUT(JB)
            TCOMP=TAVGW(TCJS(J))   !CB 9/8/06   TAVGW(JSMON(J))
          ELSEIF(TCISEG(J) < 0)THEN  
            TCOMP=TWDO(ABS(TCISEG(J)))
          ELSE

! CHECKING TO SEE IF THE MONITORING SEGMENT TCISEG IS IN THE SAME BRANCH AND WATER BODY AS THE WITHDRAWAL
            DO JJB=1,NBR
              IF (TCISEG(J) >= US(JJB) .AND. TCISEG(J) <= DS(JJB)) EXIT
            END DO
            DO JJW=1,NWB
              IF (JJB >= BS(JJW) .AND. JJB <= BE(JJW)) EXIT
            END DO

            IF (TCKLAY(J)< 0) THEN                                                                                       
              K = INT(ABS(TCKLAY(J)))
            ELSE
              DO K=KTWB(JJW),KB(TCISEG(J))
                IF (DEPTHB(K,TCISEG(J)) > TCKLAY(J)) EXIT                                                                      
              END DO
              K = MIN(K,KB(TCISEG(J)))                                                                                         
            END IF
            TCOMP=T2(K,TCISEG(J))
          ENDIF
          IF(TCYEARLY(J) == '     OFF')THEN
            DAYTEST=JDAY
          ELSE
            DAYTEST=REAL(JDAYG)+JDAY-INT(JDAY)
          END IF
          IF(DAYTEST >= TCTSRT(J) .AND. DAYTEST < TCTEND(J))THEN
               IF(TCOMP > TCTEMP(J) .AND. TCNELEV(J) > NCOUNTCW(JWD))THEN
               ! MAKING SURE THAT THE NEXT LOWER STRUCTURE FOR A PARTICULAR 'J' IS FOUND
                 DO NN=NCOUNTCW(JWD)+1,TCNELEV(J)                 
                   IF(TCELEV(J,NN) < EWD(JWD))THEN
                      NCOUNTCW(JWD)=NN
                      EWD(JWD)=TCELEV(J,NCOUNTCW(JWD))
                      EXIT
                   END IF                 
                 END DO                                               
               ELSEIF(TCOMP < TCTEMP(J) .AND.  NCOUNTCW(JWD).GT. 1)THEN
                 ! TO PREVENT THIS HAPPENING AT EACH TIME IT CHECKS IT AND HENCE OSCIALLTING BACK AND FORTH - CHECK THE TEMP AT THE UPPER OUTLET ALSO
                 IF(TCISEG(J) >  0)THEN  
                      TMOD=T2(K,TCISEG(J))                      
                    IF(TMOD < TCTEMP(J) .AND. TCELEV(J,NCOUNTCW(JWD)-1) < ELWS(IWD(JWD)))THEN                      
                      ! MAKING SURE THAT THE NEXT UPPER STRUCTURE FOR A PARTICULAR 'J' IS FOUND
                      DO NN=NCOUNTCW(JWD)-1,1,-1
                        IF(TCELEV(J,NN) > EWD(JWD))THEN
                          NCOUNTCW(JWD)=NN
                          EWD(JWD)=TCELEV(J,NCOUNTCW(JWD))
                          EXIT
                        END IF                 
                      END DO
                    ENDIF
                 END IF  ! CB 9/8/06
                 IF(TCISEG(J) == 0)THEN
! CALCULATE THE ESTIMATED OUTFLOW TEMPERATURE AT HIGHER PORTS WHEN TCOMP<TCTEMP(J), AND MOVE UP IF HIGHER PORT STILL MEETS TO CRITERIA
                   I         = MAX(CUS(JBWD(JWD)),IWD(JWD))
                   DO NN=1,NCOUNTCW(JWD)-1                     
                     CALL LATERAL_WITHDRAWAL_ESTIMATE(JWD,TEMPEST,TCELEV(J,NN))
                     IF(TEMPEST < TCTEMP(J) .AND. TCELEV(J,NN) < ELWS(IWD(JWD)))THEN
                       NCOUNTCW(JWD)=NN
                       EWD(JWD)=TCELEV(J,NCOUNTCW(JWD))
                       EXIT
                     END IF
                   END DO
                 END IF
               ENDIF
               IF(TCELEVCON(J) =='      ON' .AND. TCNELEV(J) > NCOUNTCW(JWD).AND. EWD(JWD) > ELWS(IWD(JWD)))THEN  
                 NCOUNTCW(JWD)=NCOUNTCW(JWD)+1
                 EWD(JWD)=TCELEV(J,NCOUNTCW(JWD))
               END IF
          ENDIF          
      ENDIF            
     END DO
  ENDDO
  
  NXTTCD = NXTTCD+TCDFREQ    
ENDIF  
  
RETURN
ENTRY DEALLOCATE_SELECTIVE
CLOSE(2900)
  DEALLOCATE (TCNELEV,TCJB,TCJS, TCELEV,TCTEMP,TCTEND,TCTSRT,NCOUNTC,TCISEG,TCKLAY,TCELEVCON,ELCONTSPL) 
  DEALLOCATE (TSPLTJB,TSPLTT,NOUTS,JSTSPLT,KSTRSPLT,TCYEARLY, JBMON,JSMON,TCNTR,TSPLTCNTR,JSTSPLTT) 
  DEALLOCATE (VOLM,NCOUNTCW,QWDFRAC,QSTRFRAC,MINWL)    
  DEALLOCATE(TEMPCRIT,VOLMC,DYNSEL,SELD,NXSEL,TEMP2,TSYEARLY,TSTEND,TSTSRT,DYNSF,DYNSPF,DYNSELSPLT,TEMP3,NXSELSPLT)
RETURN   

END SUBROUTINE SELECTIVE
!***********************************************************************************************************************************
!**                                            S E L E C T I V E   I N I T   U S G S                                              **
!***********************************************************************************************************************************

Module Selective1USGS
  INTEGER                                      :: numtempc, numtsplt, tempn, ng1, ng2
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: tcnelev, tcjb, tcjs, tciseg, kstrsplt, ncountcw, seld
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: tspltjb, tsseld, nouts, nout0, nout1, nout2
  INTEGER,       ALLOCATABLE, DIMENSION(:,:)   :: jstsplt, ncountc, tsprior
  REAL                                         :: nxtstr, nxttcd, nxtsplit, tcdfreq, tfrqtmp, tspltfreq
  REAL                                         :: sum_minfrac1, sum_minfrac2, qfrac1, qfrac2, tsconv
  REAL,          ALLOCATABLE, DIMENSION(:)     :: tctemp, tctend, tctsrt, tcklay, tspltt, volm, qwdfrac, tstend, tstsrt, nxsel,temp2
  REAL,          ALLOCATABLE, DIMENSION(:)     :: nxtssel, tstemp2, ewdsav, minfrac1, minfrac2, maxfrac1, maxfrac2, splt2t, splt2e
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: tcelev, tempcrit, qstrfrac, volmc
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: tsdepth, tsminfrac, tsminhead, tsmaxhead, tsmaxflow, estrsav
  CHARACTER(8)                                 :: tempc, tspltc
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: tcelevcon, tcyearly, tcntr, dynsel
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: tspltcntr, tsyearly, elcontspl, tsdynsel
  CHARACTER(5),  ALLOCATABLE, DIMENSION(:,:)   :: tstype
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: wd_active, share_flow
  LOGICAL,       ALLOCATABLE, DIMENSION(:,:)   :: no_flow, str_active
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: DYNSF, DYNTF
  CHARACTER(30)                                :: CHAR30
  REAL, ALLOCATABLE, DIMENSION(:)              :: MINWL   ! Minimum water level above centerline of outlet if TCELEVCON is ON

End Module Selective1USGS


Subroutine SelectiveInitUSGS

  Use Selective1USGS; USE MAIN
  USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  IMPLICIT NONE
  
  integer      :: ifile, nj, N, JJ
  REAL         :: DAYTEST
  character(8) :: tsshare, AID1
  CHARACTER(1) :: INFORMAT, CHAR1        
  LOGICAL      :: CSVFORMAT
  INTEGER      :: NoOutlets=30   !ZZ 4/2023 maximum no of blending outlets
  
  ifile=1949
  tavg=0.0
  tavgw=0.0
  do jb=1,nbr
    if(nstr(jb) > 0)then
      ifile=ifile+1
      write (segnum,'(i0)') jb
      segnum = adjustl(segnum)
      l      = len_trim(segnum)

      IF(RESTART_IN)THEN
        OPEN  (ifile,file='str_br'//segnum(1:l)//'.csv',POSITION='APPEND')
        JDAY1=0.0
        REWIND (IFILE)
        READ   (IFILE,'(//)',END=13)
        DO WHILE (JDAY1 < JDAY)
          READ (IFILE,*,END=13) JDAY1                !'(F10.0)'
        END DO
        BACKSPACE (IFILE)
13      JDAY1=0.0
      ELSE
        open  (ifile,file='str_br'//segnum(1:l)//'.csv',status='unknown')
        write (ifile,*)'Branch:,',jb,', # of structures:,',nstr(jb),', outlet temperatures'
        write (ifile,'("      JDAY,",<nstr(jb)>(6x,"T(C),"),<nstr(jb)>(3x,"Q(m3/s),"),<nstr(jb)>(4x,"ELEVCL,"))')
      ENDIF
    endif
  end do

  if(nwd > 0)then
    ifile=ifile+1
    IF(RESTART_IN)THEN
      OPEN  (ifile,file='wd_out.opt',POSITION='APPEND')
      JDAY1=0.0
      REWIND (IFILE)
      READ   (IFILE,'(/)',END=14)
      DO WHILE (JDAY1 < JDAY)
        READ (IFILE,'(F10.0)',END=14) JDAY1
      END DO
      BACKSPACE (IFILE)
14    JDAY1=0.0
    ELSE
      open  (ifile,file='wd_out.opt',status='unknown')
      write (ifile,*)'Withdrawals: # of withdrawals:',nwd,' outlet temperatures'
      write (ifile,'("      JDAY",<nwd>(6x,"T(C)"),<nwd>(3x,"Q(m3/s)"),<nwd>(4x,"ELEVCL"))')
    ENDIF
  end if

  !ZZ 4/2023, add a CSV format input 
  open (NUNIT,file='w2_selective.npt',status='old')   
  READ (NUNIT,'(A)') AID1 
  IF (INDEX(AID1, "$") == 0) THEN
    CSVFORMAT=.FALSE.
        ! Check for commas -- another sign that it is comma-delimited
    READ (NUNIT,'(A)') CHAR30
    DO J=1,30
      IF (CHAR30(J:J)==',') THEN
        CSVFORMAT=.TRUE.
        EXIT
      END IF
    END DO
    BACKSPACE(NUNIT)     
  ELSE
    CSVFORMAT=.TRUE.
  END IF
  IF (CSVFORMAT) THEN
    READ (NUNIT,*)
    READ (NUNIT,*)
    READ (NUNIT,*) AID1, tfrqtmp
    READ (NUNIT,*)
    READ (NUNIT,*)    
    READ (NUNIT,*) AID1, tempc, numtempc, tcdfreq
    tempc=ADJUSTR(tempc)
  ELSE  
    read (NUNIT,'(//8x,f8.0)') tfrqtmp
    read (NUNIT,'(//8x,a8,i8,f8.0)') tempc, numtempc, tcdfreq
  END IF

  IF(JDAY>TMSTRT)THEN   ! SW 8/10/2023  During restart this allows continuing output
      NXTSTR=JDAY          
      NXTTCD=JDAY   
      NXTSPLIT=JDAY
  ELSE
      NXTSTR=TMSTRT
      NXTTCD=TMSTRT   
      NXTSPLIT=TMSTRT
  ENDIF

  allocate (tcnelev(numtempc),tcjb(numtempc),tcjs(numtempc),tcelev(numtempc,100),tctemp(numtempc),tctend(numtempc),tctsrt(numtempc))
  allocate (ncountc(nst,nbr),tciseg(numtempc),tcklay(numtempc),tcelevcon(numtempc))
  Allocate (tcyearly(numtempc), tcntr(numtempc))
  allocate (volm(nwb),ncountcw(nwd),qwdfrac(nwd),qstrfrac(nst,nbr),DYNSEL(numtempc),SELD(numtempc),NXSEL(numtempc),TEMP2(numtempc))
  ALLOCATE(DYNSF(NUMTEMPC),MINWL(NUMTEMPC))
  DYNSF=.FALSE.;MINWL=0.0
  
  ncountc=0
  IF(CSVFORMAT) THEN
    READ (NUNIT,*)
    READ (NUNIT,*)
  ELSE
  read (NUNIT,'(/)')
  END IF
  do j=1,numtempc
    IF(CSVFORMAT) THEN
      READ (NUNIT,*) AID1, tcntr(j),tcjb(j),tcjs(j),tcyearly(j),tctsrt(j),tctend(j),tctemp(j),tcnelev(j),(tcelev(j,n),n=1,tcnelev(j))
      tcntr(j)=ADJUSTR(tcntr(j))
    ELSE
      read(NUNIT,'(8x,a8,i8,i8,a8,f8.0,f8.0,f8.0,i8,10(f8.0))')tcntr(j),tcjb(j),tcjs(j),tcyearly(j),tctsrt(j),tctend(j),tctemp(j),tcnelev(j),(tcelev(j,n),n=1,tcnelev(j))
    END IF
    if(tcntr(j)=='      ST')then
      tcelev(j,tcnelev(j)+1)=ESTR(tcjs(j),tcjb(j))   ! always put the original elevation as the last elevation
    else
      tcelev(j,tcnelev(j)+1)=EWD(tcjs(j))   ! always put the original elevation as the last elevation
    endif
  end do

  IF(CSVFORMAT) THEN
    READ (NUNIT,*)
    READ (NUNIT,*)
  ELSE
  read (NUNIT,'(/)')
  END IF
  do j=1,numtempc
    IF(CSVFORMAT) THEN
      READ (NUNIT,*) AID1, tciseg(j), tcklay(j), DYNSEL(J)
      DYNSEL(J)=ADJUSTR(DYNSEL(J))
    ELSE
      read (NUNIT,'(8x,i8,f8.0,A8)') tciseg(j), tcklay(j), DYNSEL(J)
    END IF
  end do
  IF(CSVFORMAT) THEN
    READ (NUNIT,*)
    READ (NUNIT,*)
  ELSE
  read (NUNIT,'(/)')
  END IF
  do j=1,numtempc
    IF(CSVFORMAT) THEN
      READ (NUNIT,*) AID1, tcelevcon(j),MINWL(j)
      tcelevcon(j)=ADJUSTR(tcelevcon(j))
    ELSE
    read (NUNIT,'(8x,a8,F8.0)') tcelevcon(j),MINWL(j)
    END IF
  end do
  IF(CSVFORMAT) THEN
    READ (NUNIT,*)
    READ (NUNIT,*)
    READ (NUNIT,*) AID1, tspltc, numtsplt, tspltfreq, tsconv
    tspltc=ADJUSTR(tspltc)
  ELSE  
  read (NUNIT,'(//8x,a8,i8,2f8.0)') tspltc, numtsplt, tspltfreq, tsconv
  END IF
  
  allocate (tsyearly(numtsplt), tstsrt(numtsplt), tstend(numtsplt), tspltjb(numtsplt), tspltt(numtsplt), nouts(numtsplt))
  allocate (ewdsav(nwd), wd_active(nwd), estrsav(nst,nbr), str_active(nst,nbr))
  !
  allocate (jstsplt(numtsplt,NoOutlets), kstrsplt(NoOutlets), tspltcntr(numtsplt), elcontspl(numtsplt))
  allocate (tsdepth(numtsplt,NoOutlets), tstype(numtsplt,NoOutlets), tsminfrac(numtsplt,NoOutlets), tsprior(numtsplt,NoOutlets))
  allocate (tsminhead(numtsplt,NoOutlets), tsmaxhead(numtsplt,NoOutlets), tsmaxflow(numtsplt,NoOutlets), no_flow(numtsplt,NoOutlets), share_flow(numtsplt))
  allocate (tsdynsel(numtsplt), tsseld(numtsplt), nxtssel(numtsplt), tstemp2(numtsplt))
  allocate (nout0(NoOutlets), nout1(NoOutlets), nout2(NoOutlets), minfrac1(NoOutlets), maxfrac1(NoOutlets), minfrac2(NoOutlets), maxfrac2(NoOutlets), splt2t(NoOutlets), splt2e(NoOutlets))  
  ALLOCATE (DYNTF(numtsplt))
  DYNTF=.FALSE.
  
   IF(CSVFORMAT) THEN
    READ (NUNIT,*)
    READ (NUNIT,*)
  ELSE
  read (NUNIT,'(/)')
  END IF
  do j=1,numtsplt
    IF(CSVFORMAT) THEN
      READ (NUNIT,*) AID1, tspltcntr(j), tspltjb(j), tsyearly(j), tstsrt(j), tstend(j), tspltt(j), tsdynsel(j), elcontspl(j), nouts(j), tsshare
      tspltcntr(j)=ADJUSTR(tspltcntr(j)); tsyearly(j)=ADJUSTR(tsyearly(j)); tsdynsel(j)=ADJUSTR(tsdynsel(j)); elcontspl(j)=ADJUSTR(elcontspl(j)); tsshare=ADJUSTR(tsshare)
    ELSE
      read (NUNIT,'(8x,a8,i8,a8,3f8.0,2a8,i8,a8)') tspltcntr(j), tspltjb(j), tsyearly(j), tstsrt(j), tstend(j),                       &
                                                tspltt(j), tsdynsel(j), elcontspl(j), nouts(j), tsshare
    END IF
    if (tspltc == '      ON') then
      if (nouts(j) < 2) then
        write (w2err, '(A,I0)') 'ERROR-- Less than two outlets specified for blending group ',j
        ERROR_OPEN = .TRUE.     ! will trigger the program to end when this subroutine is completed
        return
      else if (nouts(j) > NoOutlets) then
        write (w2err, '(A,I0)') 'ERROR-- More than ",NoOutlets," outlets specified for blending group ',j
        ERROR_OPEN = .TRUE.
        return
      end if
    end if
    share_flow(j) = tsshare == '      ON'
  end do

  IF(CSVFORMAT) THEN
    READ (NUNIT,*)
    READ (NUNIT,*)
  ELSE
  read (NUNIT,'(/)')
  END IF
  do j=1,numtsplt
    IF(CSVFORMAT) THEN
      READ (NUNIT,*) AID1, (jstsplt(j,n),n=1,nouts(j))
    ELSE
      read (NUNIT,'(8x,30i8)') (jstsplt(j,n),n=1,nouts(j))
    END IF
  end do
  IF(CSVFORMAT) THEN
    READ (NUNIT,*)
    READ (NUNIT,*)
  ELSE
  read (NUNIT,'(/)')
  END IF
  do j=1,numtsplt
    IF(CSVFORMAT) THEN
      READ (NUNIT,*) AID1, (tsdepth(j,n),n=1,nouts(j))
    ELSE
      read (NUNIT,'(8x,20f8.0)') (tsdepth(j,n),n=1,nouts(j))
    END IF
  end do
  IF(CSVFORMAT) THEN
    READ (NUNIT,*)
    READ (NUNIT,*)
  ELSE
  read (NUNIT,'(/)')
  END IF
  do j=1,numtsplt
    IF(CSVFORMAT) THEN
      READ (NUNIT,*) AID1, (tsminfrac(j,n),n=1,nouts(j))
    ELSE
      read (NUNIT,'(8x,20f8.0)') (tsminfrac(j,n),n=1,nouts(j))
    END IF
  end do
  IF(CSVFORMAT) THEN
    READ (NUNIT,*)
    READ (NUNIT,*)
  ELSE
  read (NUNIT,'(/)')
  END IF
  do j=1,numtsplt
    IF(CSVFORMAT) THEN
      READ (NUNIT,*) AID1, (tsprior(j,n),n=1,nouts(j))
    ELSE
      read (NUNIT,'(8x,20i8)') (tsprior(j,n),n=1,nouts(j))
    END IF
  end do
  IF(CSVFORMAT) THEN
    READ (NUNIT,*)
    READ (NUNIT,*)
  ELSE
  read (NUNIT,'(/)')
  END IF
  do j=1,numtsplt
    IF(CSVFORMAT) THEN
      READ (NUNIT,*) AID1, (tsminhead(j,n),n=1,nouts(j))
    ELSE
      read (NUNIT,'(8x,20f8.0)') (tsminhead(j,n),n=1,nouts(j))
    END IF
  end do
  IF(CSVFORMAT) THEN
    READ (NUNIT,*)
    READ (NUNIT,*)
  ELSE
  read (NUNIT,'(/)')
  END IF
  do j=1,numtsplt
    IF(CSVFORMAT) THEN
      READ (NUNIT,*) AID1, (tsmaxhead(j,n),n=1,nouts(j))
    ELSE
      read (NUNIT,'(8x,20f8.0)') (tsmaxhead(j,n),n=1,nouts(j))
    END IF
  end do
  IF(CSVFORMAT) THEN
    READ (NUNIT,*)
    READ (NUNIT,*)
  ELSE
  read (NUNIT,'(/)')
  END IF
  do j=1,numtsplt
    IF(CSVFORMAT) THEN
      READ (NUNIT,*) AID1, (tsmaxflow(j,n),n=1,nouts(j))
    ELSE
      read (NUNIT,'(8x,20f8.0)') (tsmaxflow(j,n),n=1,nouts(j))
    END IF
  end do

  estrsav = estr    ! Save the original structure elevations
  ewdsav  = ewd     ! Save the original withdrawal elevations
  do j=1,numtsplt
    do n=1,nouts(j)
      tstype(j,n) = "FIXED"
      if (tsdepth(j,n)   > 0.0) tstype(j,n) = "FLOAT"
      if (tsminfrac(j,n) > 1.0) tsminfrac(j,n) = 1.0    ! remove unrealistic input value
      if (tsminhead(j,n) < 0.0) tsminhead(j,n) = 0.0    ! remove unrealistic input value
      if (tsmaxhead(j,n) < 0.0) tsmaxhead(j,n) = 0.0    ! remove unrealistic input value
      if (tsmaxflow(j,n) < 0.0) tsmaxflow(j,n) = 0.0    ! remove unrealistic input value
    end do
  end do
  if (tsconv <= 0.0) tsconv = 0.005   ! constrain the convergence criterion to be > 0.0 and <= 0.1
  if (tsconv >  0.1) tsconv = 0.1

  IF(CSVFORMAT) THEN
    READ (NUNIT,*)
    READ (NUNIT,*)
    READ (NUNIT,*) AID1, tempn
  ELSE
  read (NUNIT,'(//8x,i8)') tempn
  END IF
  allocate (tempcrit(nwb,tempn),volmc(nwb,tempn))
  IF(CSVFORMAT) THEN
    READ (NUNIT,*)
    READ (NUNIT,*)
  ELSE
  read (NUNIT,'(/)')
  END IF
  do j=1,tempn
    IF(CSVFORMAT) THEN
      READ (NUNIT,*) AID1, (tempcrit(jw,j), jw=1,nwb) 
    ELSE
      read (NUNIT,'(8x,10f8.0)') (tempcrit(jw,j), jw=1,nwb)   ! Note max of 10 waterbodies
    END IF
  end do
  close (NUNIT)

  do jw=1,nwb
    ifile=ifile+1
    write (segnum,'(i0)') jw
    segnum = adjustl(segnum)
    l      = len_trim(segnum)
    IF(RESTART_IN)THEN
      OPEN  (ifile,file='Volume_wb'//segnum(1:l)//'.opt',POSITION='APPEND')
      JDAY1=0.0
      REWIND (IFILE)
      READ   (IFILE,'(/)',END=15)
      DO WHILE (JDAY1 < JDAY)
        READ (IFILE,'(F10.0)',END=15) JDAY1
      END DO
      BACKSPACE (IFILE)
15    JDAY1=0.0
    ELSE
      open (ifile,file='Volume_wb'//segnum(1:l)//'.opt',status='unknown')
      write(ifile,4315)
    END IF
  end do

4315  format("jday    Volume    ",<tempn>("Volcrit      "))

  if (tempc == '      ON') then
    do j=1,numtempc
      if (tcyearly(j) == '     OFF') then
        daytest=jday
      else
        daytest=real(jdayg)+jday-int(jday)
      end if
      if (daytest >= tctsrt(j) .and. daytest < tctend(j)) then

      ! initializing structure elevation
        if (tcntr(j) == '      ST') then
          jb = tcjb(j)                                         ! set branch index
          js = tcjs(j)                                         ! set structure index
          do nj=1,tcnelev(j)                                   ! making sure that structure is below water surface
            if (tcelev(j,nj) < elws(ds(jb))) then
              ncountc(js,jb)=nj
            ! ESTR(js,jb)=tcelev(j,ncountc(js,jb))    ! don't alter the elevation at this point.  Set it later.
              exit
            end if
          end do

      ! initializing withdrawal elevation
        else if (tcntr(j) == '      WD') then
          jwd = tcjs(j)                                        ! set withdrawal index
          do nj=1,tcnelev(j)                                   ! making sure that structure is below water surface
            if (tcelev(j,nj) < elws(iwd(jwd))) then
              ncountcw(jwd)=nj
            ! EWD(jwd)=tcelev(j,ncountcw(jwd))        ! don't alter the elevation at this point.  Set it later.
              exit
            end if
          end do
        end if
      end if

    ! Open dynamic selective withdrawal files
      if (DYNSEL(J) == '      ON') then
        WRITE (SEGNUM,'(I0)') J
        SEGNUM  = ADJUSTL(SEGNUM)
        L       = LEN_TRIM(SEGNUM)
        SELD(J) = 1009+J
        OPEN (SELD(J),FILE='dynselective'//SEGNUM(1:L)//'.npt',STATUS='OLD')
        
        READ(SELD(J),'(A1)')INFORMAT    ! SW 8/28/2019
        IF(INFORMAT=='$')DYNSF(J)=.TRUE.
    
        IF(DYNSF(J))THEN
        READ (SELD(J),'(/)')
        READ (SELD(J),*) NXSEL(J),TEMP2(J)
        tctemp(J)=TEMP2(J)
        READ (SELD(J),*) NXSEL(J),TEMP2(J)       
        
        ELSE
        READ (SELD(J),'(//1000F8.0)') NXSEL(J),TEMP2(J)
        tctemp(J)=TEMP2(J)
        READ (SELD(J),'(1000F8.0)') NXSEL(J),TEMP2(J)
        ENDIF       
        
        !READ (SELD(J),'(///1000F8.0)') NXSEL(J),TEMP2(J)
        !tctemp(J)=TEMP2(J)
        !READ (SELD(J),'(1000F8.0)') NXSEL(J),TEMP2(J)
      end if
    end do
  end if

! Open dynamic temperature target files for blending outlets
  if (tspltc == '      ON') then
    do j=1,numtsplt
      if (tsdynsel(j) == '      ON') then
        write (segnum,'(i0)') j
        segnum    = adjustl(segnum)
        L         = len_trim(segnum)
        tsseld(j) = 1009+numtempc+j
        open (tsseld(j),file='dynsplit_selective'//segnum(1:L)//'.npt',status='old')
        !
        READ(tsseld(j),'(A1)') CHAR1
        IF(CHAR1=='$') DYNTF(J)=.TRUE.
        IF(DYNTF(J))THEN
          READ (tsseld(j),'(/)')
          READ (tsseld(j),*) nxtssel(j), tstemp2(j)
          tctemp(J)=TEMP2(J)
          READ (tsseld(j),*) nxtssel(j), tstemp2(j)       
        ELSE
          READ (tsseld(j),'(//2F8.0)') nxtssel(j), tstemp2(j)
          tspltt(j) = tstemp2(j)
          READ (tsseld(j),'(2F8.0)') nxtssel(j), tstemp2(j)
        ENDIF 
      end if
    end do
  end if

! Test to see if the user specified inconsistent inputs. If so, stop with an error message.
  if (tspltc == '      ON') then
    if (tspltfreq <= 0.0) then                                                                                        !SR 12/19/2022
      write (w2err, '(A)') 'ERROR-- Update frequency for temperature blending must be greater than zero.'             !SR 12/19/2022
      ERROR_OPEN = .TRUE.                                                                                             !SR 12/19/2022
    end if                                                                                                            !SR 12/19/2022
    do j=1,numtsplt
      do n=1,nouts(j)-1
        do nj=n+1,nouts(j)
          if (jstsplt(j,n) == jstsplt(j,nj)) then
            write (w2err, '(A,I0)') 'w2_selective.npt USGS ERROR-- Duplicate split outlet numbers in group ', j
            ERROR_OPEN = .TRUE.     ! will trigger the program to end when this subroutine is completed
          end if
        end do
      end do
    end do
    do j=1,numtsplt-1
      do jj=j+1,numtsplt
        if ((tstsrt(jj) >= tstsrt(j) .and. tstsrt(jj) <  tstend(j)) .or.                                                           &
            (tstend(jj) >  tstsrt(j) .and. tstend(jj) <= tstend(j))) then
          if (tspltcntr(j) .eq. tspltcntr(jj) .and. (tspltcntr(j) .eq. '      WD' .or. tspltjb(j) == tspltjb(jj))) then
            do n=1,nouts(j)
              do nj=1,nouts(jj)
                if (jstsplt(j,n) == jstsplt(jj,nj)) then
                  write (w2err, '(A,I0,A)') 'w2_selective.npt USGS ERROR-- Split outlet number ', jstsplt(j,n), ' used in more than one group at a time.'
                  ERROR_OPEN = .TRUE.       ! will trigger the program to end when this subroutine is completed
                end if
              end do
            end do
          end if
        end if
      end do
    end do
    do j=1,numtsplt
      do n=1,nouts(j)
        if (tsprior(j,n) < -1) then
          tsprior(j,n) = -1                                                           ! reassign, rather than error   !SR 12/19/2022
          write (wrn, '(A,I0,A,I0,A)') 'w2_selective USGS: WARNING-- Priority input for outlet ', jstsplt(j,n),                                       &
                                       ' in group ', j, ' reassigned to -1.'                                          !SR 12/19/2022
          WARNING_OPEN = .TRUE.                                                                                       !SR 12/19/2022
        end if
        if (tsminhead(j,n) > 0.0 .and. tsmaxhead(j,n) > 0.0 .and. tsminhead(j,n) > tsmaxhead(j,n)) then
          write (wrn, '(A,I0,A,I0,A)') 'w2_selective.npt USGS WARNING-- Minimum and maximum head constraints for outlet ', jstsplt(j,n), ' in group ',   &
                                        j, ' are such that the outlet cannot ever be used.'
          WARNING_OPEN = .TRUE.
        end if
        if (tsdepth(j,n) > 0.0 .and. tsminhead(j,n) > 0.0 .and. tsdepth(j,n) < tsminhead(j,n)) then
          write (wrn, '(A,I0,A,I0,A)') 'w2_selective.npt USGS WARNING-- Depth of floating outlet ', jstsplt(j,n), ' in group ', j,                       &
              ' is shallower than the minimum head constraint.  To honor the head constraint, no flow is possible for that outlet.'
          WARNING_OPEN = .TRUE.
        end if
      end do
    end do

    if (tempc == '      ON') then
      do j=1,numtsplt
        do jj=1,numtempc
          if ((tctsrt(jj) >= tstsrt(j) .and. tctsrt(jj) <  tstend(j)) .or.                                                         &
              (tctend(jj) >  tstsrt(j) .and. tctend(jj) <= tstend(j))) then
            if (tspltcntr(j) .eq. tcntr(jj) .and. (tspltcntr(j) .eq. '      WD' .or. tspltjb(j) == tcjb(jj))) then
              do n=1,nouts(j)
                if (jstsplt(j,n) == tcjs(jj)) then
                  write (w2err, '(A,I0,A)') 'w2_selective.npt USGS ERROR-- Outlet number ',tcjs(jj),' used in tower and blending group at same time.'
                  ERROR_OPEN = .TRUE.       ! will trigger the program to end when this subroutine is completed
                end if
              end do
            end if
          end if
        end do
      end do
    end if
  end if

  if (tempc == '      ON') then
    do j=1,numtempc-1
      do jj=j+1,numtempc
        if ((tctsrt(jj) >= tctsrt(j) .and. tctsrt(jj) <  tctend(j)) .or.                                                           &
            (tctend(jj) >  tctsrt(j) .and. tctend(jj) <= tctend(j))) then
          if (tcntr(j) .eq. tcntr(jj) .and. (tcntr(j) .eq. '      WD' .or. tcjb(j) == tcjb(jj))) then
            if (tcjs(j) == tcjs(jj)) then
              write (w2err, '(A,I0,A)') 'w2_selective.npt USGS ERROR-- Tower outlet number ', tcjs(j), ' used more than once for overlapping dates.'
              ERROR_OPEN = .TRUE.       ! will trigger the program to end when this subroutine is completed
            end if
          end if
        end if
      end do
    end do
  end if

return

End subroutine SelectiveInitUSGS


!***********************************************************************************************************************************
!**                                                   S E L E C T I V E                                                           **
!***********************************************************************************************************************************

Subroutine SelectiveUSGS
  Use Selective1USGS; USE MAIN
  USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  IMPLICIT NONE

  integer :: ifile, j2lo, j2hi, j2max, j2min, j2pref, jj, jjw, jst, kk, ks
  integer :: n, ng0, ng1max, ng1min, nj, num_left, num_noflow, prior1, prior2
  real    :: addfrac, blendfrac, daytest, elev1, elev2, elr, etemp, etemp1, etemp2, excess_frac
  real    :: lastfrac, lastfrac2, maxelev, maxtemp, minelev, mintemp, q_notblended, qall
  real    :: sum_minfrac0, sum_maxfrac1, sum_maxfrac2, sumelev, sumfrac, sumtemp, tcomp, tempest, tmod, ttarg, wsel  
  
  str_active = .FALSE.
  wd_active  = .FALSE.
  j2pref     = 1

! Update some date variables and determine which outlets are being actively blended or adjusted.
  if (tspltc == '      ON') then
    do j=1,numtsplt
      if (tsyearly(j) == '     OFF') then
        daytest = jday
      else
        daytest = real(jdayg) + jday - int(jday)
      end if
      if (nxtsplit > tstsrt(j) .and. daytest <= tstsrt(j)) then
        nxtsplit = tstsrt(j)
      end if
      if (daytest >= tstsrt(j) .and. daytest < tstend(j)) then
        do jj=1,nouts(j)
          if (tspltcntr(j) == '      ST') then
            str_active(jstsplt(j,jj),tspltjb(j)) = .TRUE.
          else if (tspltcntr(j) == '      WD') then
            wd_active(jstsplt(j,jj)) = .TRUE.
          end if
        end do
      end if
    end do
  end if
  if (tempc == '      ON') then
    do j=1,numtempc
      if (tcyearly(j) == '     OFF') then
        daytest = jday
      else
        daytest = real(jdayg) + jday - int(jday)
      end if
      if (daytest >= tctsrt(j) .and. daytest < tctend(j)) then
        if (tcntr(j) == '      ST') then
          str_active(tcjs(j),tcjb(j)) = .TRUE.
        else if (tcntr(j) == '      WD') then
          wd_active(tcjs(j)) = .TRUE.
        end if
      end if
    end do
  end if

! Reset outlet elevations back to original values outside of control periods
  do jb=1,nbr                                                               ! swapped loops: jb, then jst             !SR 12/19/2022
    do jst=1,nstr(jb)                                                       ! changed from nst to nstr(jb)            !SR 12/19/2022
      if (.not. str_active(jst,jb)) estr(jst,jb) = estrsav(jst,jb)
    end do
  end do
  do jwd=1,nwd
    if (.not. wd_active(jwd)) ewd(jwd) = ewdsav(jwd)
  end do

! Check to see if it's time to update temperature targets and flow fractions for blended groups.
  if (tspltc=='      ON' .and. jday .ge. nxtsplit) then

  ! Update the temperature targets
    do j=1,numtsplt
      if (tsdynsel(j) == '      ON') then
        do while (jday >= nxtssel(j))
          tspltt(j) = tstemp2(j)
          IF(DYNTF(J))THEN
            read (tsseld(j),*) nxtssel(j), tstemp2(j)
          ELSE
            read (tsseld(j),'(2F8.0)') nxtssel(j), tstemp2(j)
          END IF
        end do
      end if
    end do

  ! Loop over blending groups
    do j=1,numtsplt
      qall    = 0.0                                                          ! sum up all the flows
      sumfrac = 0.0                                                          ! sum of flow fraction multipliers
      do jj=1,nouts(j)
        if (tspltcntr(j) == '      ST') then
          qall    = qall    + qstrsav(jstsplt(j,jj),tspltjb(j))                                                       !SR 06/29/2021
          sumfrac = sumfrac + qstrfrac(jstsplt(j,jj),tspltjb(j))
        else if (tspltcntr(j) == '      WD') then
          qall    = qall    + qwdsav(jstsplt(j,jj))                                                                   !SR 06/29/2021
          sumfrac = sumfrac + qwdfrac(jstsplt(j,jj))
        end if
      end do
      if (tsyearly(j) == '     OFF') then
        daytest = jday
      else
        daytest = real(jdayg) + jday - int(jday)
      end if

    ! Do blending calculations if date is in correct window or the window was just entered.
    ! If flows are zero and flow fractions are already initialized, then leave everything alone.
      if (daytest >= tstsrt(j) .and. daytest < tstend(j) .and.                                                                     &
          (qall > 0.0 .or. sumfrac < 0.001 .or. daytest < tstsrt(j) + tspltfreq)) then

      ! Do structures first
        if (tspltcntr(j) == '      ST') then
          jb = tspltjb(j)                                                    ! set branch index
          do jw=1,nwb                                                        ! set waterbody index
            if (jb >= bs(jw) .and. jb <= be(jw)) exit
          end do
          no_flow(j,:) = .FALSE.
          num_noflow   = 0
          do jj=1,nouts(j)
            jst          = jstsplt(j,jj)
            elr          = sina(jb) * dlx(ds(jb)) * 0.5
            wsel         = elws(ds(jb)) - elr                                ! compute water-surface elevation
            estr(jst,jb) = estrsav(jst,jb)                                   ! reset outlet elevation to original
            if (tstype(j,jj) .eq. "FLOAT") then
              estr(jst,jb) = wsel - tsdepth(j,jj)
            else if (estr(jst,jb) > wsel) then
              if (elcontspl(j) == '     OFF') then
                no_flow(j,jj) = .TRUE.                                       ! no flow-- high and dry
                num_noflow    = num_noflow + 1
              else
                estr(jst,jb) = wsel                                          ! poor man's floating outlet
              end if
            end if
            if (.not. no_flow(j,jj) .and. tsminhead(j,jj) > 0.0 .and. wsel - estr(jst,jb) < tsminhead(j,jj)) then
              no_flow(j,jj) = .TRUE.                                         ! minimum head criterion not met -- no flow
              num_noflow    = num_noflow + 1
            end if
            if (.not. no_flow(j,jj) .and. tsmaxhead(j,jj) > 0.0 .and. wsel - estr(jst,jb) > tsmaxhead(j,jj)) then
              no_flow(j,jj) = .TRUE.                                         ! maximum head criterion exceeded -- no flow
              num_noflow    = num_noflow + 1
            end if
            do k=ktwb(jw),kb(ds(jb))
              if (el(k,ds(jb))-elr < estr(jst,jb)) exit
            end do
            kstrsplt(jj)     = min(k-1,kb(ds(jb)))
            qstrfrac(jst,jb) = 0.0                                           ! initialize flow fractions
          end do

        ! Use priority inputs to determine which outlets to use
          prior1 = -999
          do jj=1,nouts(j)
            if (.not. no_flow(j,jj) .and. tsprior(j,jj) >= 0) then
              if (prior1 == -999 .or. tsprior(j,jj) < prior1) prior1 = tsprior(j,jj)
            end if
          end do
          prior2 = -999
          do jj=1,nouts(j)
            if (.not. no_flow(j,jj) .and. tsprior(j,jj) >= 0 .and. tsprior(j,jj) > prior1) then
              if (prior2 == -999 .or. tsprior(j,jj) < prior2) prior2 = tsprior(j,jj)
            end if
          end do

        ! Outlets with a priority of -1 get used, but are not blended
          ng0 = 0
          q_notblended = 0.0
          do jj=1,nouts(j)
            jst = jstsplt(j,jj)
            if (.not. no_flow(j,jj) .and. tsprior(j,jj) == -1) then
              ng0 = ng0 + 1
              nout0(ng0) = jj
              if (qstrsav(jst,jb) > tsmaxflow(j,jj) .and. tsmaxflow(j,jj) > 0.0) then                                 !SR 06/29/2021
                q_notblended = q_notblended + tsmaxflow(j,jj)
                qstrfrac(jst,jb) = tsmaxflow(j,jj) / qall
              else if (qall > 0.0) then
                q_notblended = q_notblended + qstrsav(jst,jb)                                                         !SR 06/29/2021
                qstrfrac(jst,jb) = qstrsav(jst,jb) / qall                                                             !SR 06/29/2021
              end if
            end if
          end do
          sum_minfrac0 = 0.0
          if (qall > 0.0) sum_minfrac0 = q_notblended / qall

        ! Outlets with priority 1 and 2 may be used and blended.
          ng1 = 0
          ng2 = 0
          sum_minfrac1 = 0.0
          sum_minfrac2 = 0.0
          sum_maxfrac1 = 0.0
          sum_maxfrac2 = 0.0
          do jj=1,nouts(j)
            if (.not. no_flow(j,jj)) then
              if (tsprior(j,jj) == prior1) then
                ng1 = ng1 + 1
                nout1(ng1) = jj
                maxfrac1(ng1) = 1.0
                if (qall > 0.0 .and. tsmaxflow(j,jj) > 0.0) maxfrac1(ng1) = min(1.0, tsmaxflow(j,jj)/qall)
                minfrac1(ng1) = tsminfrac(j,jj)
                if (tsminfrac(j,jj) < 0.0) then
                  minfrac1(ng1) = 0.0
                  if (qall > 0.0) minfrac1(ng1) = min(1.0, abs(tsminfrac(j,jj))/qall)
                end if
                if (minfrac1(ng1) > maxfrac1(ng1)) minfrac1(ng1) = maxfrac1(ng1)
                sum_minfrac1 = sum_minfrac1 + minfrac1(ng1)
                sum_maxfrac1 = sum_maxfrac1 + maxfrac1(ng1)

              else if (tsprior(j,jj) == prior2) then
                ng2 = ng2 + 1
                nout2(ng2) = jj
                maxfrac2(ng2) = 1.0
                if (qall > 0.0 .and. tsmaxflow(j,jj) > 0.0) maxfrac2(ng2) = min(1.0, tsmaxflow(j,jj)/qall)
                minfrac2(ng2) = tsminfrac(j,jj)
                if (tsminfrac(j,jj) < 0.0) then
                  minfrac2(ng2) = 0.0
                  if (qall > 0.0) minfrac2(ng2) = min(1.0, abs(tsminfrac(j,jj))/qall)
                end if
                if (minfrac2(ng2) > maxfrac2(ng2)) minfrac2(ng2) = maxfrac2(ng2)
                sum_minfrac2 = sum_minfrac2 + minfrac2(ng2)
                sum_maxfrac2 = sum_maxfrac2 + maxfrac2(ng2)
              end if
            end if
          end do

        ! If minimum flows are overspecified, then the priority 2 minimum flow fractions are decreased.
          if (ng2 > 0 .and. sum_minfrac0 + sum_minfrac1 + sum_minfrac2 > 1.0) then
            if (sum_minfrac0 + sum_minfrac1 >= 1.0) then
              ng2 = 0
              sum_minfrac2 = 0.0
            else
              do n=1,ng2
                minfrac2(n) = minfrac2(n) * (1.0 - sum_minfrac0 - sum_minfrac1) / sum_minfrac2
              end do
              sum_minfrac2 = 1.0 - sum_minfrac0 - sum_minfrac1
            end if
          end if

        ! If minimum flows are still overspecified, then the priority 1 minimum flow fractions are decreased.
          if (ng1 > 0 .and. sum_minfrac0 + sum_minfrac1 > 1.0) then
            if (sum_minfrac0 >= 1.0) then
              ng1 = 0
              sum_minfrac1 = 0.0
            else
              do n=1,ng1
                minfrac1(n) = minfrac1(n) * (1.0 - sum_minfrac0) / sum_minfrac1
              end do
              sum_minfrac1 = 1.0 - sum_minfrac0
            end if
          end if

        ! If group 1 has 3 or more outlets and group 2 has no outlets, then redistribute priorities based on elevation.
        ! Keep the highest and lowest elevation outlets in group 1, and put other active outlets into nonblended category
        ! with their minimum flows.  If ties in elevation exist, go with the first in the list.
          if (ng1 > 2 .and. ng2 == 0) then
            ng1max  = 1
            ng1min  = 1
            jst     = jstsplt(j,nout1(1))
            maxelev = estr(jst,jb)
            minelev = estr(jst,jb)
            do n=2,ng1
              jst = jstsplt(j,nout1(n))
              if (estr(jst,jb) > maxelev) then
                maxelev = estr(jst,jb)
                ng1max  = n
              else if (estr(jst,jb) < minelev) then
                minelev = estr(jst,jb)
                ng1min  = n
              end if
            end do
            blendfrac = 1.0 - sum_minfrac0 - sum_minfrac1 + minfrac1(ng1max) + minfrac1(ng1min)   ! flow frac to ng1max and ng1min
            if (maxfrac1(ng1max) + maxfrac1(ng1min) < blendfrac) then                             ! cannot handle all intended flow
              if (sum_maxfrac1 < 1.0 - sum_minfrac0) then                                         ! max flows exceeded
                write (wrn,'(A,I0,A,F0.3)') 'WARNING-- Maximum flows for outlets exceeded for group ', j, ' at day ', jday
                WARNING_OPEN = .TRUE.
                do n=1,ng1
                  minfrac1(n) = maxfrac1(n)                        ! all group 1 outlets at maximum flow              !SR 06/29/2021
                end do
              else                                                                                ! push excess to other outlets
                excess_frac = blendfrac - maxfrac1(ng1max) - maxfrac1(ng1min)                     ! ng1max and ng1min will be at max
                num_left = ng1 - 2
                do nj=1,ng1                                        ! iterative process to redistribute excess flows
                  if (num_left > 0 .and. excess_frac > 0.0) then
                    addfrac = excess_frac / num_left
                    do n=1,ng1
                      if (n .ne. ng1max .and. n .ne. ng1min .and. maxfrac1(n) - minfrac1(n) > 0.00001) then
                        if (minfrac1(n) + addfrac > maxfrac1(n)) then
                          num_left    = num_left - 1
                          excess_frac = excess_frac - (maxfrac1(n) - minfrac1(n))
                          minfrac1(n) = maxfrac1(n)
                        else
                          excess_frac = excess_frac - addfrac
                          minfrac1(n) = minfrac1(n) + addfrac
                        end if
                      end if
                    end do
                  end if
                end do
                minfrac1(ng1max) = maxfrac1(ng1max)                                               ! set to max flow   !SR 06/29/2021
                minfrac1(ng1min) = maxfrac1(ng1min)                                               ! set to max flow   !SR 06/29/2021
              end if
            end if
            do n=1,ng1                                       ! assign the other priority 1 outlets to nonblended status
              if (n .ne. ng1max .and. n .ne. ng1min) then
                ng0              = ng0 + 1
                nout0(ng0)       = nout1(n)
                jst              = jstsplt(j,nout1(n))
                sum_minfrac0     = sum_minfrac0 + minfrac1(n)
                q_notblended     = q_notblended + qall * minfrac1(n)
                qstrfrac(jst,jb) = minfrac1(n)
              end if
            end do
            ng1 = 1                                ! rearrange outlets-- one in each priority group, but same priority
            ng2 = 1
            nout2(1)     = nout1(ng1min)
            minfrac2(1)  = minfrac1(ng1min)
            maxfrac2(1)  = maxfrac1(ng1min)
            sum_minfrac2 = minfrac1(ng1min)
            sum_maxfrac2 = maxfrac1(ng1min)
            nout1(1)     = nout1(ng1max)
            minfrac1(1)  = minfrac1(ng1max)
            maxfrac1(1)  = maxfrac1(ng1max)
            sum_minfrac1 = minfrac1(ng1max)
            sum_maxfrac1 = maxfrac1(ng1max)
            prior2       = prior1
          end if

        ! If only two blended outlets, ensure that they are in separate groups.
          if (ng1 == 2 .and. ng2 == 0) then
            ng1 = 1
            ng2 = 1
            nout2(1)     = nout1(2)
            minfrac2(1)  = minfrac1(2)
            maxfrac2(1)  = maxfrac1(2)
            sum_minfrac2 = minfrac1(2)
            sum_maxfrac2 = maxfrac1(2)
            sum_minfrac1 = minfrac1(1)
            sum_maxfrac1 = maxfrac1(1)
            prior2       = prior1
          end if

        ! Begin the blending decisions.
        ! No usable outlets.  All flow fractions remain at zero.
          if (nouts(j) == num_noflow) then
            write (wrn,'(A,I0,A,F0.3)') 'Warning-- All outlets dry or unusable for group ', j, ' at day ', jday
            WARNING_OPEN = .TRUE.

        ! Only nonblended outlets.
          else if (nouts(j) == ng0) then
            write (wrn,'(A,I0,A,F0.3)') 'Warning-- Only nonblended outlets present in group ', j, ' at day ', jday
            WARNING_OPEN = .TRUE.

        ! Only one blended outlet.  It gets all of the blended flow, but must not exceed its maximum flow criterion.
          else if (ng1 + ng2 == 1) then
            jst = jstsplt(j,nout1(1))
            qstrfrac(jst,jb) = 1.0 - sum_minfrac0
            if (qall - q_notblended > tsmaxflow(j,nout1(1)) .and. tsmaxflow(j,nout1(1)) > 0.0) then
              qstrfrac(jst,jb) = tsmaxflow(j,nout1(1)) / qall
              write (wrn,'(A,A,I0,A,I0,A,F0.3)') 'Warning-- Total release flow rate decreased to comply with maximum flow ',       &
                                                 'criterion for structure ', jst, ' in group ', j, ' at day ', jday
              WARNING_OPEN = .TRUE.
            end if

        ! Minimum flows comprise entire release.  No blending calculations required.
          else if (abs(1.0 - sum_minfrac0 - sum_minfrac1 - sum_minfrac2) <= 0.000001) then
            do n=1,ng1
              jst = jstsplt(j,nout1(n))
              qstrfrac(jst,jb) = minfrac1(n)
            end do
            do n=1,ng2
              jst = jstsplt(j,nout2(n))
              qstrfrac(jst,jb) = minfrac2(n)
            end do

        ! Nonblended outlets all at maximum flows.  No blending calculations required.                                !SR 06/29/2021
          else if (ng1 > 0 .and. sum_minfrac1 == sum_maxfrac1 .and. (ng2 == 0 .or. sum_minfrac2 == sum_maxfrac2)) then !SR 06/29/2021
            do n=1,ng1                                                                                                !SR 06/29/2021
              jst = jstsplt(j,nout1(n))                                                                               !SR 06/29/2021
              qstrfrac(jst,jb) = minfrac1(n)                                                                          !SR 06/29/2021
            end do                                                                                                    !SR 06/29/2021
            do n=1,ng2                                                                                                !SR 06/29/2021
              jst = jstsplt(j,nout2(n))                                                                               !SR 06/29/2021
              qstrfrac(jst,jb) = minfrac2(n)                                                                          !SR 06/29/2021
            end do                                                                                                    !SR 06/29/2021

        ! More than one usable outlet, and blending among priority 1 outlet(s) and priority 2 outlet(s) required.
          else
            id = ds(jb)                                                      ! needed for downstream_withdrawal_estimate
            kt = ktwb(jw)                                                    ! needed for downstream_withdrawal_estimate

          ! Warn the user if maximum flow criteria are likely to decrease the specified outflows.  Use 0.999999 for round-off.
            if (sum_minfrac0 + sum_maxfrac1 + sum_maxfrac2 < 0.999999) then                                           !SR 06/29/2021
              write (wrn,'(A,A,I0,A,F0.3)') 'WARNING-- Total release flow rate may be decreased to comply with maximum flow ',     &
                                            'criteria for structures in group ', j, ' at day ', jday
              WARNING_OPEN = .TRUE.
            end if

          ! Set the initial release fractions.  Ensure that maximum flows are not exceeded.
            qfrac1 = sum_minfrac1 + 0.5 * (1.0 - sum_minfrac0 - sum_minfrac1 - sum_minfrac2)
            qfrac1 = min(sum_maxfrac1, qfrac1)
            qfrac2 = 1.0 - sum_minfrac0 - qfrac1
            if (qfrac2 > sum_maxfrac2) then
              excess_frac = qfrac2 - sum_maxfrac2
              if (qfrac1 + excess_frac <= sum_maxfrac1) then
                qfrac1 = qfrac1 + excess_frac
              else
                qfrac1 = sum_maxfrac1
              end if
              qfrac2 = sum_maxfrac2
            end if
            j2pref = 1
            call Set_Flow_Fracs(j, jb, j2pref)                               ! set flow fractions; redistribute if maxfrac exceeded

          ! If priority 2 outlets are not sharing flows, identify the two with the highest and lowest elevations
            if (.not. share_flow(j) .and. ng2 > 1) then
              j2hi    = 1
              j2lo    = 1
              jst     = jstsplt(j,nout2(1))
              maxelev = estr(jst,jb)
              minelev = estr(jst,jb)
              do n=2,ng2
                jst = jstsplt(j,nout2(n))
                if (estr(jst,jb) > maxelev) then
                  maxelev = estr(jst,jb)
                  j2hi    = n
                else if (estr(jst,jb) < minelev) then
                  minelev = estr(jst,jb)
                  j2lo    = n
                end if
              end do
            end if

          ! Get weighted blend of all nonblended release temperatures
            ttarg = tspltt(j)
            if (sum_minfrac0 > 0.0) then
              sumtemp = 0.0
              do n=1,ng0
                jst = jstsplt(j,nout0(n))
                qstr(jst,jb) = qall * qstrfrac(jst,jb)
                if (qstr(jst,jb) > 0.0) then
                  call downstream_withdrawal_estimate(jst,etemp,estr(jst,jb))
                else
                  etemp = t2(kstrsplt(nout0(n)),ds(jb))                      ! Use temperature at outlet elevation if no flow
                end if
                sumtemp = sumtemp + qstrfrac(jst,jb) * etemp
              end do
              etemp = sumtemp / sum_minfrac0
              ttarg = (ttarg - sum_minfrac0 * etemp) / (1.0 - sum_minfrac0)  ! New temperature target for blended releases
            end if

          ! Need an iterative approach because released T depends on Q
            lastfrac = qfrac1
            do jj=1,8                                                        ! Maximum of eight iterations
              lastfrac2 = lastfrac
              lastfrac  = qfrac1

              sumtemp = 0.0
              sumelev = 0.0
              do n=1,ng1                                                     ! Get weighted temp and elevation for group 1
                jst = jstsplt(j,nout1(n))
                qstr(jst,jb) = qall * qstrfrac(jst,jb)
                if (qstr(jst,jb) > 0.0) then
                  call downstream_withdrawal_estimate(jst,etemp,estr(jst,jb))
                else
                  etemp = t2(kstrsplt(nout1(n)),ds(jb))                      ! Use temperature at outlet elevation if no flow
                end if
                if (qfrac1 > 0.0) then
                  sumtemp = sumtemp + qstrfrac(jst,jb) * etemp
                  sumelev = sumelev + qstrfrac(jst,jb) * estr(jst,jb)
                else
                  sumtemp = sumtemp + etemp
                  sumelev = sumelev + estr(jst,jb)
                end if
              end do
              if (qfrac1 > 0.0) then
                etemp1 = sumtemp / qfrac1                                    ! Weighted temperature from group 1 outlets
                elev1  = sumelev / qfrac1                                    ! Weighted elevation of group 1 outlets
              else
                etemp1 = sumtemp / ng1
                elev1  = sumelev / ng1
              end if

              if (share_flow(j) .or. ng2 < 2) then                           ! Get weighted temp and elevation for group 2
                sumtemp = 0.0                                                ! ...when flows are shared among outlets
                sumelev = 0.0
                do n=1,ng2
                  jst = jstsplt(j,nout2(n))
                  qstr(jst,jb) = qall * qstrfrac(jst,jb)
                  if (qstr(jst,jb) > 0.0) then
                    call downstream_withdrawal_estimate(jst,etemp,estr(jst,jb))
                  else
                    etemp = t2(kstrsplt(nout2(n)),ds(jb))                    ! Use temperature at outlet elevation if no flow
                  end if
                  if (qfrac2 > 0.0) then
                    sumtemp = sumtemp + qstrfrac(jst,jb) * etemp
                    sumelev = sumelev + qstrfrac(jst,jb) * estr(jst,jb)
                  else
                    sumtemp = sumtemp + etemp
                    sumelev = sumelev + estr(jst,jb)
                  end if
                end do
                if (qfrac2 > 0.0) then
                  etemp2 = sumtemp / qfrac2                                  ! Weighted temperature from group 2 outlets
                  elev2  = sumelev / qfrac2                                  ! Weighted elevation of group 2 outlets
                else
                  etemp2 = sumtemp / ng2
                  elev2  = sumelev / ng2
                end if

              else                                                           ! ...and when flows are not shared
                if (qfrac2 == 0.0) then
                  do n=1,ng2
                    jst       = jstsplt(j,nout2(n))
                    splt2t(n) = t2(kstrsplt(nout2(n)),ds(jb))
                    splt2e(n) = estr(jst,jb)
                  end do
                else
                  do nj=1,ng2                                                ! Find the temperatures produced in group 2
                    sumtemp = 0.0                                            ! by testing when each outlet is preferred
                    sumelev = 0.0
                    call Set_Flow_Fracs2(j, jb, nj)
                    do n=1,ng2
                      jst = jstsplt(j,nout2(n))
                      qstr(jst,jb) = qall * qstrfrac(jst,jb)
                      if (qstr(jst,jb) > 0.0) then
                        call downstream_withdrawal_estimate(jst,etemp,estr(jst,jb))
                      else
                        etemp = t2(kstrsplt(nout2(n)),ds(jb))
                      end if
                      sumtemp = sumtemp + qstrfrac(jst,jb) * etemp
                      sumelev = sumelev + qstrfrac(jst,jb) * estr(jst,jb)
                    end do
                    splt2t(nj) = sumtemp / qfrac2
                    splt2e(nj) = sumelev / qfrac2
                  end do
                end if
                j2max   = 1
                j2min   = 1
                maxtemp = splt2t(1)
                mintemp = splt2t(1)
                do n=2,ng2
                  if (splt2t(n) > maxtemp) then
                    maxtemp = splt2t(n)
                    j2max   = n
                  else if (splt2t(n) < mintemp) then
                    mintemp = splt2t(n)
                    j2min   = n
                  end if
                end do
                if (ttarg < etemp1 - 0.001) then                             ! need a colder temp from group 2
                  if (maxtemp - mintemp > 0.001) then
                    etemp2 = splt2t(j2min)                                   ! preferred outlet is the coldest one
                    elev2  = splt2e(j2min)
                    j2pref = j2min
                  else
                    etemp2 = splt2t(j2lo)                                    ! preferred outlet is the lowest one
                    elev2  = splt2e(j2lo)
                    j2pref = j2lo
                  end if
                else if (ttarg > etemp1 + 0.001) then                        ! need a warmer temp from group 2
                  if (maxtemp - mintemp > 0.001) then
                    etemp2 = splt2t(j2max)                                   ! preferred outlet is the warmest one
                    elev2  = splt2e(j2max)
                    j2pref = j2max
                  else
                    etemp2 = splt2t(j2hi)                                    ! preferred outlet is the highest one
                    elev2  = splt2e(j2hi)
                    j2pref = j2hi
                  end if
                else
                  etemp2 = splt2t(1)                                         ! if temp is close to target, choose first outlet
                  elev2  = splt2e(1)
                  j2pref = 1
                end if
              end if

            ! Target temperature is less than either outlet temperature.
              if (ttarg < etemp1 .and. ttarg < etemp2) then
                qfrac1 = sum_minfrac1                                        ! default for if/then cases
                if (abs(etemp1 - etemp2) < 0.001) then
                  if (prior1 == prior2) then                                 ! If each outlet has the same priority level, then...
                    if (elev1 <= elev2) then                                 ! Choose lower outlet if both have same temperature.
                      qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2
                    end if
                  else if (prior1 < prior2) then
                    qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2               ! Choose higher priority outlet if temps are the same.
                  end if
                else if (etemp1 < etemp2) then                               ! If temps are different, choose the one closer
                  qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2                 ! to target temperature
                end if

            ! Target temperature is greater than either outlet temperature.
              else if (ttarg > etemp1 .and. ttarg > etemp2) then
                qfrac1 = sum_minfrac1                                        ! default for if/then cases
                if (abs(etemp1 - etemp2) < 0.001) then
                  if (prior1 == prior2) then                                 ! If each outlet has the same priority level, then...
                    if (elev1 >= elev2) then                                 ! Choose upper outlet if both have same temperature.
                      qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2
                    end if
                  else if (prior1 < prior2) then
                    qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2               ! Choose higher priority outlet if temps are the same.
                  end if
                else if (etemp1 > etemp2) then                               ! If temps are different, choose the one closer
                  qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2                 ! to target temperature
                end if

            ! Target temperature is essentially the same as the two outlet temperatures.
              else if (abs(etemp1 - etemp2) < 0.001) then
                qfrac1 = sum_minfrac1                                        ! default for if/then cases
                if (prior1 == prior2) then                                   ! If each outlet has the same priority level, then...
                  qfrac1 = sum_minfrac1 + 0.5 * (1.0 - sum_minfrac0 - sum_minfrac1 - sum_minfrac2)    ! Split the flow equally.
                else if (prior1 < prior2) then
                  qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2                 ! Choose higher priority outlet if temps are the same.
                end if

            ! Target temperature is between the two outlet temperatures.
              else
                qfrac1 = (1.0 - sum_minfrac0) * abs((ttarg-etemp2)/(etemp1-etemp2+NONZERO))
                qfrac1 = max(sum_minfrac1, qfrac1)
                qfrac1 = min(1.0 - sum_minfrac0 - sum_minfrac2, qfrac1)
              end if
              qfrac1 = min(sum_maxfrac1, qfrac1)
              qfrac2 = 1.0 - sum_minfrac0 - qfrac1
              if (qfrac2 > sum_maxfrac2) then
                excess_frac = qfrac2 - sum_maxfrac2
                if (qfrac1 + excess_frac <= sum_maxfrac1) then
                  qfrac1 = qfrac1 + excess_frac
                else
                  qfrac1 = sum_maxfrac1
                end if
                qfrac2 = sum_maxfrac2
              end if

            ! Set flow fractions for individual outlets and redistribute flows if maximum flow fractions exceeded.
              call Set_Flow_Fracs(j, jb, j2pref)

            ! Exit the loop if the latest flow fraction calculation agrees with the previous one.
            ! Exit if no flow, because no iteration requried in that case.
              if (abs(lastfrac - qfrac1) < tsconv .or. qall == 0.0) exit
            end do

          ! Check to see if iterative solution did not converge.
            if (abs(lastfrac - qfrac1) >= tsconv .and. qall > 0.0) then
              write (wrn,'(A,F0.3,3(A,F0.4))') 'Flow fraction calculations not converging at day ', jday,                          &
                                               '  Current: ', qfrac1, ' Last: ', lastfrac, ' Next-to-last: ', lastfrac2
              WARNING_OPEN = .TRUE.

            ! Check to see if the iterative solution is unstable.  If so, use priorities to assign releases.
            ! Criteria:  change is at least 0.1 and most recent change is in opposite direction from previous change.
              if (abs(lastfrac - qfrac1) >= 0.1 .and. (qfrac1-lastfrac)*(lastfrac-lastfrac2) < 0.0) then
                qfrac1 = sum_minfrac1                                        ! default for if/then cases
                if (prior1 < prior2) then                                    ! group 1 is higher priority
                  qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2
                else                                                         ! else, fulfill minima and split the rest
                  qfrac1 = sum_minfrac1 + 0.5 * (1.0 - sum_minfrac0 - sum_minfrac1 - sum_minfrac2)
                end if
                qfrac1 = min(sum_maxfrac1, qfrac1)
                qfrac2 = 1.0 - sum_minfrac0 - qfrac1
                if (qfrac2 > sum_maxfrac2) then
                  excess_frac = qfrac2 - sum_maxfrac2
                  if (qfrac1 + excess_frac <= sum_maxfrac1) then
                    qfrac1 = qfrac1 + excess_frac
                  else
                    qfrac1 = sum_maxfrac1
                  end if
                  qfrac2 = sum_maxfrac2
                end if
                call Set_Flow_Fracs(j, jb, j2pref)                           ! set flow fractions; redistribute if maxfrac exceeded
              end if
            end if
          end if

        ! Set final flows to go with the flow fractions.  May not be necessary, but do it anyway.
          do jj=1,nouts(j)
            qstr(jstsplt(j,jj),jb) = qall * qstrfrac(jstsplt(j,jj),jb)
          end do


      ! Do Withdrawals next
        else if (tspltcntr(j) == '      WD') then
          jwd = jstsplt(j,1)                                                 ! assume withdrawals are from same branch and waterbody
          do jb=1,nbr
            if (iwd(jwd) >= us(jb) .and. iwd(jwd) <= ds(jb)) exit
          end do
          do jw=1,nwb
            if (jb >= bs(jw) .and. jb <= be(jw)) exit
          end do
          no_flow(j,:) = .FALSE.
          num_noflow   = 0
          do jj=1,nouts(j)
            jwd      = jstsplt(j,jj)
            elr      = sina(jb) * dlx(iwd(jwd)) * 0.5
            wsel     = elws(iwd(jwd)) - elr                                  ! compute water-surface elevation
            ewd(jwd) = ewdsav(jwd)                                           ! reset outlet elevation to original
            if (tstype(j,jj) .eq. "FLOAT") then
              ewd(jwd) = wsel - tsdepth(j,jj)
            else if (ewd(jwd) > wsel) then
              if (elcontspl(j) == '     OFF') then
                no_flow(j,jj) = .TRUE.                                       ! no flow-- high and dry
                num_noflow    = num_noflow + 1
              else
                ewd(jwd) = wsel                                              ! poor man's floating outlet
              end if
            end if
            if (.not. no_flow(j,jj) .and. tsminhead(j,jj) > 0.0 .and. wsel - ewd(jwd) < tsminhead(j,jj)) then
              no_flow(j,jj) = .TRUE.                                         ! minimum head criterion not met -- no flow
              num_noflow    = num_noflow + 1
            end if
            if (.not. no_flow(j,jj) .and. tsmaxhead(j,jj) > 0.0 .and. wsel - ewd(jwd) > tsmaxhead(j,jj)) then
              no_flow(j,jj) = .TRUE.                                         ! maximum head criterion exceeded -- no flow
              num_noflow    = num_noflow + 1
            end if
            do k=ktwb(jw),kb(iwd(jwd))
              if (el(k,iwd(jwd))-elr < ewd(jwd)) exit
            end do
            kstrsplt(jj) = min(k-1,kb(iwd(jwd)))
            qwdfrac(jwd) = 0.0                                               ! initialize flow fractions
          end do

        ! Use priority inputs to determine which outlets to use
          prior1 = -999
          do jj=1,nouts(j)
            if (.not. no_flow(j,jj) .and. tsprior(j,jj) >= 0) then
              if (prior1 == -999 .or. tsprior(j,jj) < prior1) prior1 = tsprior(j,jj)
            end if
          end do
          prior2 = -999
          do jj=1,nouts(j)
            if (.not. no_flow(j,jj) .and. tsprior(j,jj) >= 0 .and. tsprior(j,jj) > prior1) then
              if (prior2 == -999 .or. tsprior(j,jj) < prior2) prior2 = tsprior(j,jj)
            end if
          end do

        ! Outlets with a priority of -1 get used, but are not blended
          ng0 = 0
          q_notblended = 0.0
          do jj=1,nouts(j)
            jwd = jstsplt(j,jj)
            if (.not. no_flow(j,jj) .and. tsprior(j,jj) == -1) then
              ng0 = ng0 + 1
              nout0(ng0) = jj
              if (qwdsav(jwd) > tsmaxflow(j,jj) .and. tsmaxflow(j,jj) > 0.0) then                                     !SR 06/29/2021
                q_notblended = q_notblended + tsmaxflow(j,jj)
                qwdfrac(jwd) = tsmaxflow(j,jj) / qall
              else if (qall > 0.0) then
                q_notblended = q_notblended + qwdsav(jwd)                                                             !SR 06/29/2021
                qwdfrac(jwd) = qwdsav(jwd) / qall                                                                     !SR 06/29/2021
              end if
            end if
          end do
          sum_minfrac0 = 0.0
          if (qall > 0.0) sum_minfrac0 = q_notblended / qall

        ! Outlets with priority 1 and 2 may be used and blended.
          ng1 = 0
          ng2 = 0
          sum_minfrac1 = 0.0
          sum_minfrac2 = 0.0
          sum_maxfrac1 = 0.0
          sum_maxfrac2 = 0.0
          do jj=1,nouts(j)
            if (.not. no_flow(j,jj)) then
              if (tsprior(j,jj) == prior1) then
                ng1 = ng1 + 1
                nout1(ng1) = jj
                maxfrac1(ng1) = 1.0
                if (qall > 0.0 .and. tsmaxflow(j,jj) > 0.0) maxfrac1(ng1) = min(1.0, tsmaxflow(j,jj)/qall)
                minfrac1(ng1) = tsminfrac(j,jj)
                if (tsminfrac(j,jj) < 0.0) then
                  minfrac1(ng1) = 0.0
                  if (qall > 0.0) minfrac1(ng1) = min(1.0, abs(tsminfrac(j,jj))/qall)
                end if
                if (minfrac1(ng1) > maxfrac1(ng1)) minfrac1(ng1) = maxfrac1(ng1)
                sum_minfrac1 = sum_minfrac1 + minfrac1(ng1)
                sum_maxfrac1 = sum_maxfrac1 + maxfrac1(ng1)

              else if (tsprior(j,jj) == prior2) then
                ng2 = ng2 + 1
                nout2(ng2) = jj
                maxfrac2(ng2) = 1.0
                if (qall > 0.0 .and. tsmaxflow(j,jj) > 0.0) maxfrac2(ng2) = min(1.0, tsmaxflow(j,jj)/qall)
                minfrac2(ng2) = tsminfrac(j,jj)
                if (tsminfrac(j,jj) < 0.0) then
                  minfrac2(ng2) = 0.0
                  if (qall > 0.0) minfrac2(ng2) = min(1.0, abs(tsminfrac(j,jj))/qall)
                end if
                if (minfrac2(ng2) > maxfrac2(ng2)) minfrac2(ng2) = maxfrac2(ng2)
                sum_minfrac2 = sum_minfrac2 + minfrac2(ng2)
                sum_maxfrac2 = sum_maxfrac2 + maxfrac2(ng2)
              end if
            end if
          end do

        ! If minimum flows are overspecified, then the priority 2 minimum flow fractions are decreased.
          if (ng2 > 0 .and. sum_minfrac0 + sum_minfrac1 + sum_minfrac2 > 1.0) then
            if (sum_minfrac0 + sum_minfrac1 >= 1.0) then
              ng2 = 0
              sum_minfrac2 = 0.0
            else
              do n=1,ng2
                minfrac2(n) = minfrac2(n) * (1.0 - sum_minfrac0 - sum_minfrac1) / sum_minfrac2
              end do
              sum_minfrac2 = 1.0 - sum_minfrac0 - sum_minfrac1
            end if
          end if

        ! If minimum flows are still overspecified, then the priority 1 minimum flow fractions are decreased.
          if (ng1 > 0 .and. sum_minfrac0 + sum_minfrac1 > 1.0) then
            if (sum_minfrac0 >= 1.0) then
              ng1 = 0
              sum_minfrac1 = 0.0
            else
              do n=1,ng1
                minfrac1(n) = minfrac1(n) * (1.0 - sum_minfrac0) / sum_minfrac1
              end do
              sum_minfrac1 = 1.0 - sum_minfrac0
            end if
          end if

        ! If group 1 has 3 or more outlets and group 2 has no outlets, then redistribute priorities based on elevation.
        ! Keep the highest and lowest elevation outlets in group 1, and put other active outlets into nonblended category
        ! with their minimum flows.  If ties in elevation exist, go with the first in the list.
          if (ng1 > 2 .and. ng2 == 0) then
            ng1max  = 1
            ng1min  = 1
            jwd     = jstsplt(j,nout1(1))
            maxelev = ewd(jwd)
            minelev = ewd(jwd)
            do n=2,ng1
              jwd = jstsplt(j,nout1(n))
              if (ewd(jwd) > maxelev) then
                maxelev = ewd(jwd)
                ng1max  = n
              else if (ewd(jwd) < minelev) then
                minelev = ewd(jwd)
                ng1min  = n
              end if
            end do
            blendfrac = 1.0 - sum_minfrac0 - sum_minfrac1 + minfrac1(ng1max) + minfrac1(ng1min)
            if (maxfrac1(ng1max) + maxfrac1(ng1min) < blendfrac) then
              if (sum_maxfrac1 < 1.0 - sum_minfrac0) then
                write (wrn,'(A,I0,A,F0.3)') 'Warning-- Maximum flows for outlets exceeded for group ', j, ' at day ', jday
                WARNING_OPEN = .TRUE.
                do n=1,ng1
                  minfrac1(n) = maxfrac1(n)                        ! all group 1 outlets at maximum flow              !SR 06/29/2021
                end do
              else                                                                                ! push excess to other outlets
                excess_frac = blendfrac - maxfrac1(ng1max) - maxfrac1(ng1min)                     ! ng1max and ng1min will be at max
                num_left = ng1 - 2
                do nj=1,ng1                                        ! iterative process to redistribute excess flows
                  if (num_left > 0 .and. excess_frac > 0.0) then
                    addfrac = excess_frac / num_left
                    do n=1,ng1
                      if (n .ne. ng1max .and. n .ne. ng1min .and. maxfrac1(n) - minfrac1(n) > 0.00001) then
                        if (minfrac1(n) + addfrac > maxfrac1(n)) then
                          num_left    = num_left - 1
                          excess_frac = excess_frac - (maxfrac1(n) - minfrac1(n))
                          minfrac1(n) = maxfrac1(n)
                        else
                          excess_frac = excess_frac - addfrac
                          minfrac1(n) = minfrac1(n) + addfrac
                        end if
                      end if
                    end do
                  end if
                end do
                minfrac1(ng1max) = maxfrac1(ng1max)                                               ! set to max flow   !SR 06/29/2021
                minfrac1(ng1min) = maxfrac1(ng1min)                                               ! set to max flow   !SR 06/29/2021
              end if
            end if
            do n=1,ng1                                       ! assign the other priority 1 outlets to nonblended status
              if (n .ne. ng1max .and. n .ne. ng1min) then
                ng0          = ng0 + 1
                nout0(ng0)   = nout1(n)
                jwd          = jstsplt(j,nout1(n))
                sum_minfrac0 = sum_minfrac0 + minfrac1(n)
                q_notblended = q_notblended + qall * minfrac1(n)
                qwdfrac(jwd) = minfrac1(n)
              end if
            end do
            ng1 = 1                                ! rearrange outlets-- one in each priority group, but same priority
            ng2 = 1
            nout2(1)     = nout1(ng1min)
            minfrac2(1)  = minfrac1(ng1min)
            maxfrac2(1)  = maxfrac1(ng1min)
            sum_minfrac2 = minfrac1(ng1min)
            sum_maxfrac2 = maxfrac1(ng1min)
            nout1(1)     = nout1(ng1max)
            minfrac1(1)  = minfrac1(ng1max)
            maxfrac1(1)  = maxfrac1(ng1max)
            sum_minfrac1 = minfrac1(ng1max)
            sum_maxfrac1 = maxfrac1(ng1max)
            prior2       = prior1
          end if

        ! If only two blended outlets, ensure that they are in separate groups.
          if (ng1 == 2 .and. ng2 == 0) then
            ng1 = 1
            ng2 = 1
            nout2(1)     = nout1(2)
            minfrac2(1)  = minfrac1(2)
            maxfrac2(1)  = maxfrac1(2)
            sum_minfrac2 = minfrac1(2)
            sum_maxfrac2 = maxfrac1(2)
            sum_minfrac1 = minfrac1(1)
            sum_maxfrac1 = maxfrac1(1)
            prior2       = prior1
          end if


        ! Begin the blending decisions.
        ! No usable outlets.  All flow fractions remain at zero.
          if (nouts(j) == num_noflow) then
            write (wrn,'(A,I0,A,F0.3)') 'Warning-- All outlets dry or unusable for group ', j, ' at day ', jday
            WARNING_OPEN = .TRUE.

        ! Only nonblended outlets.
          else if (nouts(j) == ng0) then
            write (wrn,'(A,I0,A,F0.3)') 'Warning-- Only nonblended outlets present in group ', j, ' at day ', jday
            WARNING_OPEN = .TRUE.

        ! Only one blended outlet.  It gets all of the blended flow, but must not exceed its maximum flow criterion.
          else if (ng1 + ng2 == 1) then
            jwd = jstsplt(j,nout1(1))
            qwdfrac(jwd) = 1.0 - sum_minfrac0
            if (qall - q_notblended > tsmaxflow(j,nout1(1)) .and. tsmaxflow(j,nout1(1)) > 0.0) then
              qwdfrac(jwd) = tsmaxflow(j,nout1(1)) / qall
              write (wrn,'(A,A,I0,A,I0,A,F0.3)') 'Warning-- Total release flow rate decreased to comply with maximum flow ',       &
                                                 'criterion for withdrawal ', jwd, ' in group ', j, ' at day ', jday
              WARNING_OPEN = .TRUE.
            end if

        ! Minimum flows comprise entire release.  No blending calculations required.
          else if (abs(1.0 - sum_minfrac0 - sum_minfrac1 - sum_minfrac2) <= 0.000001) then
            do n=1,ng1
              jwd = jstsplt(j,nout1(n))
              qwdfrac(jwd) = minfrac1(n)
            end do
            do n=1,ng2
              jwd = jstsplt(j,nout2(n))
              qwdfrac(jwd) = minfrac2(n)
            end do

        ! Nonblended outlets all at maximum flows.  No blending calculations required.                                !SR 06/29/2021
          else if (ng1 > 0 .and. sum_minfrac1 == sum_maxfrac1 .and. (ng2 == 0 .or. sum_minfrac2 == sum_maxfrac2)) then !SR 06/29/2021
            do n=1,ng1                                                                                                !SR 06/29/2021
              jwd = jstsplt(j,nout1(n))                                                                               !SR 06/29/2021
              qwdfrac(jwd) = minfrac1(n)                                                                              !SR 06/29/2021
            end do                                                                                                    !SR 06/29/2021
            do n=1,ng2                                                                                                !SR 06/29/2021
              jwd = jstsplt(j,nout2(n))                                                                               !SR 06/29/2021
              qwdfrac(jwd) = minfrac2(n)                                                                              !SR 06/29/2021
            end do                                                                                                    !SR 06/29/2021

        ! More than one usable outlet, and blending among priority 1 outlet(s) and priority 2 outlet(s) required.
          else
            I  = iwd(jstsplt(j,nout1(1)))                                    ! needed for lateral_withdrawal_estimate
            kt = ktwb(jw)                                                    ! needed for lateral_withdrawal_estimate

          ! Warn the user if maximum flow criteria are likely to decrease the specified outflows.  Use 0.999999 for round-off.
            if (sum_minfrac0 + sum_maxfrac1 + sum_maxfrac2 < 0.999999) then                                           !SR 06/29/2021
              write (wrn,'(A,A,I0,A,F0.3)') 'WARNING-- Total release flow rate may be decreased to comply with maximum flow ',     &
                                            'criteria for withdrawals in group ', j, ' at day ', jday
              WARNING_OPEN = .TRUE.
            end if

          ! Set the initial release fractions.  Ensure that maximum flows are not exceeded.
            qfrac1 = sum_minfrac1 + 0.5 * (1.0 - sum_minfrac0 - sum_minfrac1 - sum_minfrac2)
            qfrac1 = min(sum_maxfrac1, qfrac1)
            qfrac2 = 1.0 - sum_minfrac0 - qfrac1
            if (qfrac2 > sum_maxfrac2) then
              excess_frac = qfrac2 - sum_maxfrac2
              if (qfrac1 + excess_frac <= sum_maxfrac1) then
                qfrac1 = qfrac1 + excess_frac
              else
                qfrac1 = sum_maxfrac1
              end if
              qfrac2 = sum_maxfrac2
            end if
            j2pref = 1
            call Set_Flow_Fracs(j, jb, j2pref)                               ! set flow fractions; redistribute if maxfrac exceeded

          ! If priority 2 outlets are not sharing flows, identify the two with the highest and lowest elevations
            if (.not. share_flow(j) .and. ng2 > 1) then
              j2hi    = 1
              j2lo    = 1
              jwd     = jstsplt(j,nout2(1))
              maxelev = ewd(jwd)
              minelev = ewd(jwd)
              do n=2,ng2
                jwd = jstsplt(j,nout2(n))
                if (ewd(jwd) > maxelev) then
                  maxelev = ewd(jwd)
                  j2hi    = n
                else if (ewd(jwd) < minelev) then
                  minelev = ewd(jwd)
                  j2lo    = n
                end if
              end do
            end if

          ! Get weighted blend of all nonblended release temperatures
            ttarg = tspltt(j)
            if (sum_minfrac0 > 0.0) then
              sumtemp = 0.0
              do n=1,ng0
                jwd = jstsplt(j,nout0(n))
                qwd(jwd) = qall * qwdfrac(jwd)
                if (qwd(jwd) > 0.0) then
                  call lateral_withdrawal_estimate(jwd,etemp,ewd(jwd))       ! Get an estimate of the temperature of the outflow
                else
                  etemp = t2(kstrsplt(nout0(n)),iwd(jwd))                    ! Use temperature at outlet elevation if no flow
                end if
                sumtemp = sumtemp + qwdfrac(jwd) * etemp
              end do
              etemp = sumtemp / sum_minfrac0
              ttarg = (ttarg - sum_minfrac0 * etemp) / (1.0 - sum_minfrac0)  ! New temperature target for blended releases
            end if

          ! Need an iterative approach because released T depends on Q
            lastfrac = qfrac1
            do jj=1,8                                                        ! Maximum of eight iterations
              lastfrac2 = lastfrac
              lastfrac  = qfrac1

              sumtemp = 0.0
              sumelev = 0.0
              do n=1,ng1                                                     ! Get weighted temp and elevation for group 1
                jwd = jstsplt(j,nout1(n))
                qwd(jwd) = qall * qwdfrac(jwd)
                if (qwd(jwd) > 0.0) then
                  call lateral_withdrawal_estimate(jwd,etemp,ewd(jwd))
                else
                  etemp = t2(kstrsplt(nout1(n)),iwd(jwd))                    ! Use temperature at outlet elevation if no flow
                end if
                if (qfrac1 > 0.0) then
                  sumtemp = sumtemp + qwdfrac(jwd) * etemp
                  sumelev = sumelev + qwdfrac(jwd) * ewd(jwd)
                else
                  sumtemp = sumtemp + etemp
                  sumelev = sumelev + ewd(jwd)
                end if
              end do
              if (qfrac1 > 0.0) then
                etemp1 = sumtemp / qfrac1                                    ! Weighted temperature from group 1 outlets
                elev1  = sumelev / qfrac1                                    ! Weighted elevation of group 1 outlets
              else
                etemp1 = sumtemp / ng1
                elev1  = sumelev / ng1
              end if

              if (share_flow(j) .or. ng2 < 2) then                           ! Get weighted temp and elevation for group 2
                sumtemp = 0.0                                                ! ...when flows are shared among outlets
                sumelev = 0.0
                do n=1,ng2
                  jwd = jstsplt(j,nout2(n))
                  qwd(jwd) = qall * qwdfrac(jwd)
                  if (qwd(jwd) > 0.0) then
                    call lateral_withdrawal_estimate(jwd,etemp,ewd(jwd))
                  else
                    etemp = t2(kstrsplt(nout2(n)),iwd(jwd))                  ! Use temperature at outlet elevation if no flow
                  end if
                  if (qfrac2 > 0.0) then
                    sumtemp = sumtemp + qwdfrac(jwd) * etemp
                    sumelev = sumelev + qwdfrac(jwd) * ewd(jwd)
                  else
                    sumtemp = sumtemp + etemp
                    sumelev = sumelev + ewd(jwd)
                  end if
                end do
                if (qfrac2 > 0.0) then
                  etemp2 = sumtemp / qfrac2                                  ! Weighted temperature from group 2 outlets
                  elev2  = sumelev / qfrac2                                  ! Weighted elevation of group 2 outlets
                else
                  etemp2 = sumtemp / ng2
                  elev2  = sumelev / ng2
                end if

              else                                                           ! ...and when flows are not shared
                if (qfrac2 == 0.0) then
                  do n=1,ng2
                    jwd       = jstsplt(j,nout2(n))
                    splt2t(n) = t2(kstrsplt(nout2(n)),iwd(jwd))
                    splt2e(n) = ewd(jwd)
                  end do
                else
                  do nj=1,ng2                                                ! Find the temperatures produced in group 2
                    sumtemp = 0.0                                            ! by testing when each outlet is preferred
                    sumelev = 0.0
                    call Set_Flow_Fracs2(j, jb, nj)
                    do n=1,ng2
                      jwd = jstsplt(j,nout2(n))
                      qwd(jwd) = qall * qwdfrac(jwd)
                      if (qwd(jwd) > 0.0) then
                        call lateral_withdrawal_estimate(jwd,etemp,ewd(jwd))
                      else
                        etemp = t2(kstrsplt(nout2(n)),iwd(jwd))
                      end if
                      sumtemp = sumtemp + qwdfrac(jwd) * etemp
                      sumelev = sumelev + qwdfrac(jwd) * ewd(jwd)
                    end do
                    splt2t(nj) = sumtemp / qfrac2
                    splt2e(nj) = sumelev / qfrac2
                  end do
                end if
                j2max   = 1
                j2min   = 1
                maxtemp = splt2t(1)
                mintemp = splt2t(1)
                do n=2,ng2
                  if (splt2t(n) > maxtemp) then
                    maxtemp = splt2t(n)
                    j2max   = n
                  else if (splt2t(n) < mintemp) then
                    mintemp = splt2t(n)
                    j2min   = n
                  end if
                end do
                if (ttarg < etemp1 - 0.001) then                             ! need a colder temp from group 2
                  if (maxtemp - mintemp > 0.001) then
                    etemp2 = splt2t(j2min)                                   ! preferred outlet is the coldest one
                    elev2  = splt2e(j2min)
                    j2pref = j2min
                  else
                    etemp2 = splt2t(j2lo)                                    ! preferred outlet is the lowest one
                    elev2  = splt2e(j2lo)
                    j2pref = j2lo
                  end if
                else if (ttarg > etemp1 + 0.001) then                        ! need a warmer temp from group 2
                  if (maxtemp - mintemp > 0.001) then
                    etemp2 = splt2t(j2max)                                   ! preferred outlet is the warmest one
                    elev2  = splt2e(j2max)
                    j2pref = j2max
                  else
                    etemp2 = splt2t(j2hi)                                    ! preferred outlet is the highest one
                    elev2  = splt2e(j2hi)
                    j2pref = j2hi
                  end if
                else
                  etemp2 = splt2t(1)                                         ! if temp is close to target, choose first outlet
                  elev2  = splt2e(1)
                  j2pref = 1
                end if
              end if

            ! Target temperature is less than either outlet temperature.
              if (ttarg < etemp1 .and. ttarg < etemp2) then
                qfrac1 = sum_minfrac1                                        ! default for if/then cases
                if (abs(etemp1 - etemp2) < 0.001) then
                  if (prior1 == prior2) then                                 ! If each outlet has the same priority level, then...
                    if (elev1 <= elev2) then                                 ! Choose lower outlet if both have same temperature.
                      qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2
                    end if
                  else if (prior1 < prior2) then
                    qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2               ! Choose higher priority outlet if temps are the same.
                  end if
                else if (etemp1 < etemp2) then                               ! If temps are different, choose the one closer
                  qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2                 ! to target temperature
                end if

            ! Target temperature is greater than either outlet temperature.
              else if (ttarg > etemp1 .and. ttarg > etemp2) then
                qfrac1 = sum_minfrac1                                        ! default for if/then cases
                if (abs(etemp1 - etemp2) < 0.001) then
                  if (prior1 == prior2) then                                 ! If each outlet has the same priority level, then...
                    if (elev1 >= elev2) then                                 ! Choose upper outlet if both have same temperature.
                      qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2
                    end if
                  else if (prior1 < prior2) then
                    qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2               ! Choose higher priority outlet if temps are the same.
                  end if
                else if (etemp1 > etemp2) then                               ! If temps are different, choose the one closer
                  qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2                 ! to target temperature
                end if

            ! Target temperature is essentially the same as the two outlet temperatures.
              else if (abs(etemp1 - etemp2) < 0.001) then
                qfrac1 = sum_minfrac1                                        ! default for if/then cases
                if (prior1 == prior2) then                                   ! If each outlet has the same priority level, then...
                  qfrac1 = sum_minfrac1 + 0.5 * (1.0 - sum_minfrac0 - sum_minfrac1 - sum_minfrac2)    ! Split the flow equally.
                else if (prior1 < prior2) then
                  qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2                 ! Choose higher priority outlet if temps are the same.
                end if

            ! Target temperature is between the two outlet temperatures.
              else
                qfrac1 = (1.0 - sum_minfrac0) * abs((ttarg-etemp2)/(etemp1-etemp2+NONZERO))
                qfrac1 = max(sum_minfrac1, qfrac1)
                qfrac1 = min(1.0 - sum_minfrac0 - sum_minfrac2, qfrac1)
              end if
              qfrac1 = min(sum_maxfrac1, qfrac1)
              qfrac2 = 1.0 - sum_minfrac0 - qfrac1
              if (qfrac2 > sum_maxfrac2) then
                excess_frac = qfrac2 - sum_maxfrac2
                if (qfrac1 + excess_frac <= sum_maxfrac1) then
                  qfrac1 = qfrac1 + excess_frac
                else
                  qfrac1 = sum_maxfrac1
                end if
                qfrac2 = sum_maxfrac2
              end if

            ! Set flow fractions for individual outlets and redistribute flows if maximum flow fractions exceeded.
              call Set_Flow_Fracs(j, jb, j2pref)

            ! Exit the loop if the latest flow fraction calculation agrees with the previous one.
            ! Exit if no flow, because no iteration requried in that case.
              if (abs(lastfrac - qfrac1) < tsconv .or. qall == 0.0) exit
            end do

          ! Check to see if iterative solution did not converge.
            if (abs(lastfrac - qfrac1) >= tsconv .and. qall > 0.0) then
              write (wrn,'(A,F0.3,3(A,F0.4))') 'Flow fraction calculations not converging at day ', jday,                          &
                                               '  Current: ', qfrac1, ' Last: ', lastfrac, ' Next-to-last: ', lastfrac2
              WARNING_OPEN = .TRUE.

            ! Check to see if the iterative solution is unstable.  If so, use priorities to assign releases.
            ! Criteria:  change is at least 0.1 and most recent change is in opposite direction from previous change.
              if (abs(lastfrac - qfrac1) >= 0.1 .and. (qfrac1-lastfrac)*(lastfrac-lastfrac2) < 0.0) then
                qfrac1 = sum_minfrac1                                        ! default for if/then cases
                if (prior1 < prior2) then                                    ! group 1 is higher priority
                  qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2
                else                                                         ! else, fulfill minima and split the rest
                  qfrac1 = sum_minfrac1 + 0.5 * (1.0 - sum_minfrac0 - sum_minfrac1 - sum_minfrac2)
                end if
                qfrac1 = min(sum_maxfrac1, qfrac1)
                qfrac2 = 1.0 - sum_minfrac0 - qfrac1
                if (qfrac2 > sum_maxfrac2) then
                  excess_frac = qfrac2 - sum_maxfrac2
                  if (qfrac1 + excess_frac <= sum_maxfrac1) then
                    qfrac1 = qfrac1 + excess_frac
                  else
                    qfrac1 = sum_maxfrac1
                  end if
                  qfrac2 = sum_maxfrac2
                end if
                call Set_Flow_Fracs(j, jb, j2pref)                           ! set flow fractions; redistribute if maxfrac exceeded
              end if
            end if
          end if

        ! Set final flows to go with the flow fractions.  May not be necessary, but do it anyway.
          do jj=1,nouts(j)
            qwd(jstsplt(j,jj)) = qall * qwdfrac(jstsplt(j,jj))
          end do
        end if
      end if
    end do

    nxtsplit = nxtsplit + tspltfreq
  end if

! Use the flow fractions to set flows in blended groups.
  if (tspltc=='      ON') then
    do j=1,numtsplt
      if (tsyearly(j) == '     OFF') then
        daytest = jday
      else
        daytest = real(jdayg) + jday - int(jday)
      end if
      if (daytest >= tstsrt(j) .and. daytest < tstend(j)) then
        qall = 0.0

      ! Do structures first
        if (tspltcntr(j) == '      ST') then
          do jj=1,nouts(j)
            qall = qall + qstrsav(jstsplt(j,jj),tspltjb(j))                         ! sum up all the flows            !SR 06/29/2021
          end do
          do jj=1,nouts(j)                                                          ! set the flows and honor the maximum flow
            jst = jstsplt(j,jj)
            qstr(jst,tspltjb(j)) = qstrfrac(jst,tspltjb(j)) * qall
            if (tsmaxflow(j,jj) > 0.0 .and. qstr(jst,tspltjb(j)) > tsmaxflow(j,jj)) qstr(jst,tspltjb(j)) = tsmaxflow(j,jj)
          end do

      ! Do Withdrawals next
        else if (tspltcntr(j) == '      WD') then
          do jj=1,nouts(j)
            qall = qall + qwdsav(jstsplt(j,jj))                                     ! sum up all the flows            !SR 06/29/2021
          end do
          do jj=1,nouts(j)                                                          ! set the flows and honor the maximum flow
            jwd = jstsplt(j,jj)
            qwd(jwd) = qwdfrac(jwd) * qall
            if (tsmaxflow(j,jj) > 0.0 .and. qwd(jwd) > tsmaxflow(j,jj)) qwd(jwd) = tsmaxflow(j,jj)
          end do
        end if
      end if
    end do
  end if

! Output some results.
  if (jday.ge.nxtstr) then
    nxtstr = nxtstr+tfrqtmp
    ifile=1949
    do jb=1,nbr
      if (nstr(jb) > 0) then
        ifile=ifile+1
        write (ifile,'(f10.4,",",<nstr(jb)>(f10.2,","),<nstr(jb)>(f10.2,","),<nstr(jb)>(f10.2,","))') jday,(tavg(i,jb),i=1,nstr(jb)),(qstr(i,jb),i=1,nstr(jb)),(estr(i,jb),i=1,nstr(jb))                   ! SW 8/28/2019
      end if
    end do
    if (nwd > 0) then
      ifile=ifile+1
      write (ifile,'(f10.4,<nwd>f10.2,<nwd>f10.2,<nwd>f10.2)') jday,(tavgw(i),i=1,nwd),(qwd(i),i=1,nwd),(ewd(i),i=1,nwd)
    end if

  ! computing reservoir volume and volume below 'tempcrit'
    volmc=0.0
    volm=0.0
    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        DO I=cus(jb),ds(jb)
          volm(jw) = volm(jw) +BH2(KT,I)*DLX(I)
          DO K=kt+1,kb(i)
            volm(jw) = volm(jw)+BH(K,I)*DLX(I)
          END DO
          do kk=1,tempn
            if(t2(kt,i).le.tempcrit(jw,kk))volmc(jw,kk) = volmc(jw,kk)+BH2(KT,I)*DLX(I)
            DO K=kt+1,kb(i)
              if(t2(k,i).le.tempcrit(jw,kk))volmc(jw,kk) = volmc(jw,kk)+BH(K,I)*DLX(I)
            END DO
          end do
        end do
      end do

      ifile=ifile+1
      write(ifile,5315)jday,volm(jw),(volmc(jw,kk), kk=1,tempn)
5315  format(f8.2,100(g12.4,g12.4))
    end do
  end if

! Check elevations and status of temperature control towers.
  if (tempc == '      ON' .and. jday .ge. nxttcd) then

  ! Update the temperature targets
    DO J=1,NUMTEMPC
      IF (DYNSEL(J) == '      ON') THEN
        SELD(J)=1009+J
        DO WHILE (JDAY >= NXSEL(J))
          tctemp(j) = TEMP2(J)
          IF(DYNSF(J))THEN
              READ (SELD(J),*) NXSEL(J),TEMP2(J)
          ELSE
              READ (SELD(J),'(1000F8.0)') NXSEL(J),TEMP2(J)
          ENDIF   
          !READ (SELD(J),'(1000F8.0)') NXSEL(J),TEMP2(J)
        END DO
      END IF
    END DO

    do j=1,numtempc

    ! Structures
      if (tcntr(j) == '      ST') then
        js = tcjs(j)                                 ! set structure index
        jb = tcjb(j)                                 ! set branch index
        DO JW=1,NWB                                  ! set waterbody index
          IF (JB >= BS(JW) .AND. JB <= BE(JW)) EXIT
        END DO

        if (tciseg(j) .eq. 0) then
          tcomp = tavg(js,jb)               !cb 9/8/06
        else if (tciseg(j) < 0) then
          tcomp = twdo(abs(tciseg(j)))      ! sw 11/26/10
        else

        ! Check to see if the monitoring segment tciseg is in the same branch and waterbody as the structure
          DO JJB=1,NBR
            IF (tciseg(j) >= US(JJB) .AND. tciseg(j) <= DS(JJB)) exit
          end do
          DO JJW=1,NWB
            IF (JjB >= BS(JJW) .AND. JJB <= BE(JJW)) EXIT
          END DO

          IF (tcklay(j) < 0) THEN
            K = INT(ABS(tcklay(j)))
          ELSE
            DO K=KTWB(JJW),KB(tciseg(j))
              IF (DEPTHB(K,tciseg(j)) > tcklay(j)) EXIT
            END DO
            K = MIN(K,KB(tciseg(j)))
          END IF
          tcomp = t2(k,tciseg(j))
        end if
        if (tcyearly(j) == '     OFF') then
          daytest = jday
        else
          daytest = real(jdayg) + jday - int(jday)
        end if
        if (daytest >= tctsrt(j) .and. daytest < tctend(j)) then
          estr(js,jb) = tcelev(j,ncountc(js,jb))                  ! initialize the structure elevation

          if (tcomp > tctemp(j) .and. tcnelev(j) > ncountc(js,jb)) then
          ! making sure that the next lower structure for a particular 'j' is found
            do nj=ncountc(js,jb)+1,tcnelev(j)
              if (tcelev(j,nj) < estr(js,jb)) then
                ncountc(js,jb) = nj
                estr(js,jb)    = tcelev(j,ncountc(js,jb))
                exit
              end if
            end do

          else if (tcomp < tctemp(j) .and. ncountc(js,jb) .gt. 1) then
          ! to prevent this happening at each time it checks it and hence oscillating back and forth - check the temp at the upper outlet also
            if (tciseg(j) > 0) then
              if (jb .eq. jjb) then
                wsel = elws(ds(jb)) - sina(jb) * dlx(ds(jb)) * 0.5                        ! compute water-surface elevation !SR 03/24/13
                do ks=ktwb(jw),kb(ds(jb))
                ! if (depthb(ks,tciseg(j)) > tcelev(j,ncountc(js,jb)-1)) exit             !??can't be right-- SR 03/24/13
                  if (wsel - depthb(ks,tciseg(j)) < tcelev(j,ncountc(js,jb)-1)) exit      !SR 03/24/13
                end do
                ks   = min(ks,kb(tciseg(j)))
                tmod = t2(ks,ds(jb))
              else
                tmod = t2(k,tciseg(j))
              end if
              if (tmod < tctemp(j) .and. tcelev(j,ncountc(js,jb)-1) < elws(ds(jb))) then
              ! making sure that the next upper structure for a particular 'j' is found
                do nj=ncountc(js,jb)-1,1,-1
                  if (tcelev(j,nj) > estr(js,jb)) then
                    ncountc(js,jb) = nj
                    estr(js,jb)    = tcelev(j,ncountc(js,jb))
                    exit
                  end if
                end do
              end if

            else if (tciseg(j) .eq. 0) then
            ! calculate the estimated outflow temperature at higher ports when tcomp < tctemp(j),
            ! and move up if higher port still meets criteria - this doesn't happen when tciseg < 0
              do nj=1,ncountc(js,jb)-1
                id = ds(jb)
                kt = ktwb(jw)
                call downstream_withdrawal_estimate(js,tempest,tcelev(j,nj))
                if (tempest < tctemp(j) .and. tcelev(j,nj) < elws(ds(jb))) then
                  ncountc(js,jb) = nj
                  estr(js,jb)    = tcelev(j,ncountc(js,jb))
                  exit
                end if
              end do
            end if
          end if
          if (tcelevcon(j) == '      ON' .and. tcnelev(j) > ncountc(js,jb) .and. estr(js,jb) > elws(ds(jb)-MINWL(J))) then
            ncountc(js,jb) = ncountc(js,jb)+1
            estr(js,jb)    = tcelev(j,ncountc(js,jb))
          end if
        end if

    ! Withdrawals
      else if (tcntr(j) == '      WD') then
        jwd = tcjs(j)
        if (tciseg(j) .eq. 0) then
        ! tcomp = tout(jb)
          tcomp = tavgw(tcjs(j))   !cb 9/8/06
        else if (tciseg(j) < 0) then
          tcomp = twdo(abs(tciseg(j)))
        else

        ! checking to see if the monitoring segment tciseg is in the same branch and water body as the withdrawal
          DO JJB=1,NBR
            IF (tciseg(j) >= US(JJB) .AND. tciseg(j) <= DS(JJB)) exit
          end do
          DO JJW=1,NWB
            IF (JjB >= BS(JJW) .AND. JJB <= BE(JJW)) EXIT
          END DO

          IF (tcklay(j) < 0) THEN
            K = INT(ABS(tcklay(j)))
          ELSE
            DO K=KTWB(JJW),KB(tciseg(j))
              IF (DEPTHB(K,tciseg(j)) > tcklay(j)) EXIT
            END DO
            K = MIN(K,KB(tciseg(j)))
          END IF
          tcomp = t2(k,tciseg(j))
        end if
        if (tcyearly(j) == '     OFF') then
          daytest = jday
        else
          daytest = real(jdayg) + jday - int(jday)
        end if
        if (daytest >= tctsrt(j) .and. daytest < tctend(j)) then
          ewd(jwd) = tcelev(j,ncountcw(jwd))                  ! initialize the withdrawal elevation

          if (tcomp > tctemp(j) .and. tcnelev(j) > ncountcw(jwd)) then
          ! making sure that the next lower structure for a particular 'j' is found
            do nj=ncountcw(jwd)+1,tcnelev(j)
              if (tcelev(j,nj) < ewd(jwd)) then
                ncountcw(jwd) = nj
                ewd(jwd)      = tcelev(j,ncountcw(jwd))
                exit
              end if
            end do

          else if (tcomp < tctemp(j) .and. ncountcw(jwd) .gt. 1) then
          ! to prevent this happening at each time it checks it and hence oscillating back and forth - check the temp at the upper outlet also
            if (tciseg(j) > 0) then
              tmod = t2(k,tciseg(j))
              if (tmod < tctemp(j) .and. tcelev(j,ncountcw(jwd)-1) < elws(iwd(jwd))) then
              ! making sure that the next upper structure for a particular 'j' is found
                do nj=ncountcw(jwd)-1,1,-1
                  if (tcelev(j,nj) > ewd(jwd)) then
                    ncountcw(jwd) = nj
                    ewd(jwd)      = tcelev(j,ncountcw(jwd))
                    exit
                  end if
                end do
              end if

            else if (tciseg(j) == 0) then
            ! calculate estimated outflow temperature at higher ports when tcomp < tctemp(j),
            ! and move up if higher port still meets criteria
              I = MAX(CUS(JBWD(JWD)),IWD(JWD))
              DO JJB=1,NBR
                IF (I >= US(JJB) .AND. I <= DS(JJB)) exit
              end do
              DO JJW=1,NWB
                IF (JJB >= BS(JJW) .AND. JJB <= BE(JJW)) EXIT
              END DO
              kt = ktwb(jjw)
              do nj=1,ncountcw(jwd)-1
                call lateral_withdrawal_estimate(jwd,tempest,tcelev(j,nj))
                if (tempest < tctemp(j) .and. tcelev(j,nj) < elws(iwd(jwd))) then
                  ncountcw(jwd) = nj
                  ewd(jwd)      = tcelev(j,ncountcw(jwd))
                  exit
                end if
              end do
            end if
          end if
          if (tcelevcon(j) == '      ON' .and. tcnelev(j) > ncountcw(jwd) .and. ewd(jwd) > elws(iwd(jwd))-MINWL(J)) then
            ncountcw(jwd) = ncountcw(jwd)+1
            ewd(jwd)      = tcelev(j,ncountcw(jwd))
          end if
        end if
      end if
    end do

    nxttcd = nxttcd + tcdfreq
  end if
return
		  
		  
ENTRY DEALLOCATE_SELECTIVEUSGS
  DEAllocate (tcnelev,tcjb,tcjs, tcelev,tctemp,tctend,tctsrt,ncountc,tciseg,tcklay,tcelevcon,elcontspl)
  DEAllocate (tspltjb,tspltt,nouts,jstsplt,kstrsplt,tcyearly, tcntr,tspltcntr)
  DEallocate (volm,ncountcw,qwdfrac,qstrfrac,MINWL)
  DEallocate (tempcrit,volmc,DYNSEL,SELD,NXSEL,TEMP2,TSYEARLY,TSTEND,TSTSRT)
  deallocate (tsdepth, tstype, tsminfrac, tsprior, tsminhead, tsmaxhead, tsmaxflow, no_flow)
  deallocate (tsdynsel, tsseld, nxtssel, tstemp2, ewdsav, estrsav, share_flow, wd_active, str_active)
  deallocate (nout0, nout1, nout2, minfrac1, minfrac2, maxfrac1, maxfrac2, splt2t, splt2e)
  DEALLOCATE (DYNSF, DYNTF)
RETURN

End Subroutine SelectiveUSGS


!***********************************************************************************************************************************
!**                                                S E T _ F L O W _ F R A C T I O N S                                            **
!***********************************************************************************************************************************
                                                                              ! Entire routine added/modified by S. Rounds, 06/26/13
Subroutine Set_Flow_Fractions
  Use Selective1USGS                                                          ! This routine sets the flow fractions, and then
  IMPLICIT NONE                                                               ! redistributes flow to other outlets in the same
  integer :: n, nj, nexcess, j, jb, jst, jwd, j2pref                          ! group when one or more outlets exceed their maximum
  real    :: excess_frac, addfrac                                             ! flow rates.  Excess flow that cannot be accommodated
Return                                                                        ! within each group will be discarded.


Entry Set_Flow_Fracs(j, jb, j2pref)
  excess_frac = 0.0                                                           ! Excess flow above maximum release rates
  nexcess     = 0                                                             ! Number of group 1 outlets exceeding maximum rates

! Set and rebalance release fractions for group 1 structures
  if (tspltcntr(j) == '      ST') then
    do n=1,ng1                                                                ! Find maxed-out outlets in group 1; set flows to max
      jst = jstsplt(j,nout1(n))
      qstrfrac(jst,jb) = minfrac1(n) + (qfrac1 - sum_minfrac1) / ng1
      if (qstrfrac(jst,jb) > maxfrac1(n)) then
        nexcess = nexcess + 1
        excess_frac = excess_frac + qstrfrac(jst,jb) - maxfrac1(n)
        qstrfrac(jst,jb) = maxfrac1(n)
      end if
    end do
    if (excess_frac > 0.0 .and. ng1 - nexcess > 0) then                       ! Redistribute excess flow to other outlets in group
      do nj=1,ng1                                                             ! Iterative process, in case others get maxed-out
        if (ng1 == nexcess .or. excess_frac <= 0.00001) exit
        addfrac = excess_frac / (ng1 - nexcess)
        do n=1,ng1
          jst = jstsplt(j,nout1(n))
          if (maxfrac1(n) - qstrfrac(jst,jb) > 0.00001) then
            if (qstrfrac(jst,jb) + addfrac > maxfrac1(n)) then
              nexcess = nexcess + 1
              excess_frac = excess_frac - (maxfrac1(n) - qstrfrac(jst,jb))
              qstrfrac(jst,jb) = maxfrac1(n)
            else
              excess_frac = excess_frac - addfrac
              qstrfrac(jst,jb) = qstrfrac(jst,jb) + addfrac
            end if
          end if
        end do
      end do
    end if

! Set and rebalance release fractions for group 1 withdrawals
  else
    do n=1,ng1                                                                ! Find maxed-out outlets in group 1; set flows to max
      jwd = jstsplt(j,nout1(n))
      qwdfrac(jwd) = minfrac1(n) + (qfrac1 - sum_minfrac1) / ng1
      if (qwdfrac(jwd) > maxfrac1(n)) then
        nexcess = nexcess + 1
        excess_frac = excess_frac + qwdfrac(jwd) - maxfrac1(n)
        qwdfrac(jwd) = maxfrac1(n)
      end if
    end do
    if (excess_frac > 0.0 .and. ng1 - nexcess > 0) then                       ! Redistribute excess flow to other outlets in group
      do nj=1,ng1                                                             ! Iterative process, in case others get maxed-out
        if (ng1 == nexcess .or. excess_frac <= 0.00001) exit
        addfrac = excess_frac / (ng1 - nexcess)
        do n=1,ng1
          jwd = jstsplt(j,nout1(n))
          if (maxfrac1(n) - qwdfrac(jwd) > 0.00001) then
            if (qwdfrac(jwd) + addfrac > maxfrac1(n)) then
              nexcess = nexcess + 1
              excess_frac = excess_frac - (maxfrac1(n) - qwdfrac(jwd))
              qwdfrac(jwd) = maxfrac1(n)
            else
              excess_frac = excess_frac - addfrac
              qwdfrac(jwd) = qwdfrac(jwd) + addfrac
            end if
          end if
        end do
      end do
    end if
  end if

Entry Set_Flow_Fracs2(j, jb, j2pref)                                          ! Separate entry just for group 2 outlets

  if (j2pref == 0) j2pref = 1                                                 ! Preferred outlet number, if not sharing
  excess_frac = 0.0                                                           ! Excess flow above maximum release rates
  nexcess     = 0                                                             ! Number of group 2 outlets exceeding maximum rates

! Set and rebalance release fractions for group 2 structures
  if (tspltcntr(j) == '      ST') then
    do n=1,ng2                                                                ! Find maxed-out outlets in group 2; set flows to max
      jst = jstsplt(j,nout2(n))
      if (.not. share_flow(j) .and. ng2 > 1) then                             ! Direct flow to preferred outlet if not shared
        if (n == j2pref) then
          qstrfrac(jst,jb) = qfrac2 - sum_minfrac2 + minfrac2(n)
        else
          qstrfrac(jst,jb) = minfrac2(n)
        end if
      else
        qstrfrac(jst,jb) = minfrac2(n) + (qfrac2 - sum_minfrac2) / ng2
      end if
      if (qstrfrac(jst,jb) > maxfrac2(n)) then
        nexcess = nexcess + 1
        excess_frac = excess_frac + qstrfrac(jst,jb) - maxfrac2(n)
        qstrfrac(jst,jb) = maxfrac2(n)
      end if
    end do
    if (excess_frac > 0.0 .and. ng2 - nexcess > 0) then                       ! Redistribute excess flow to other outlets in group
      do nj=1,ng2                                                             ! Iterative process, in case others get maxed-out
        if (ng2 == nexcess .or. excess_frac <= 0.00001) exit
        addfrac = excess_frac / (ng2 - nexcess)
        do n=1,ng2
          jst = jstsplt(j,nout2(n))
          if (maxfrac2(n) - qstrfrac(jst,jb) > 0.00001) then
            if (qstrfrac(jst,jb) + addfrac > maxfrac2(n)) then
              nexcess = nexcess + 1
              excess_frac = excess_frac - (maxfrac2(n) - qstrfrac(jst,jb))
              qstrfrac(jst,jb) = maxfrac2(n)
            else
              excess_frac = excess_frac - addfrac
              qstrfrac(jst,jb) = qstrfrac(jst,jb) + addfrac
            end if
          end if
        end do
      end do
    end if

! Set and rebalance release fractions for group 2 withdrawals
  else if (tspltcntr(j) == '      WD') then
    do n=1,ng2                                                                ! Find maxed-out outlets in group 2; set flows to max
      jwd = jstsplt(j,nout2(n))
      if (.not. share_flow(j) .and. ng2 > 1) then                             ! Direct flow to preferred outlet if not shared
        if (n == j2pref) then
          qwdfrac(jwd) = qfrac2 - sum_minfrac2 + minfrac2(n)
        else
          qwdfrac(jwd) = minfrac2(n)
        end if
      else
        qwdfrac(jwd) = minfrac2(n) + (qfrac2 - sum_minfrac2) / ng2
      end if
      if (qwdfrac(jwd) > maxfrac2(n)) then
        nexcess = nexcess + 1
        excess_frac = excess_frac + qwdfrac(jwd) - maxfrac2(n)
        qwdfrac(jwd) = maxfrac2(n)
      end if
    end do
    if (excess_frac > 0.0 .and. ng2 - nexcess > 0) then                       ! Redistribute excess flow to other outlets in group
      do nj=1,ng2                                                             ! Iterative process, in case others get maxed-out
        if (ng2 == nexcess .or. excess_frac <= 0.00001) exit
        addfrac = excess_frac / (ng2 - nexcess)
        do n=1,ng2
          jwd = jstsplt(j,nout2(n))
          if (maxfrac2(n) - qwdfrac(jwd) > 0.00001) then
            if (qwdfrac(jwd) + addfrac > maxfrac2(n)) then
              nexcess = nexcess + 1
              excess_frac = excess_frac - (maxfrac2(n) - qwdfrac(jwd))
              qwdfrac(jwd) = maxfrac2(n)
            else
              excess_frac = excess_frac - addfrac
              qwdfrac(jwd) = qwdfrac(jwd) + addfrac
            end if
          end if
        end do
      end do
    end if
  end if
Return
  End Subroutine Set_Flow_Fractions

!***********************************************************************************************************************************
! Selective Withdrawal and Forcast Computations for the Shasta Dam TCD                
! Developed by WaterCourse Engineering, 11/18/2019, 8/17/2021                       
! Incorporated into V5.0 and Updated by                                     
! Zhong Zhang, PSU                                                 
! 10/15/2024                                                        
!
! FORCAST input parameters:
! FCST	   ON/OFF switch
! LWRGATE	 Adjustment value to target temperature (oC)
! BEGIN_HI The begin julian day of high-gate open  restriction period
! END_HI	 The end julian day of high-gate open  restriction period
! N_TRIES	 The number of required attempts
! INGATE	 The low gate number to begin forecast blending
! MINSUB	 The level of minimum submersion above gate invert for opening gates (m)
   
!***********************************************************************************************************************************
!**                                               S E L E C T I V E  T C D  I N I T                                               **
!***********************************************************************************************************************************  
Module Selective1TCD
  INTEGER                                      :: numtempc, numtsplt, tempn, ng1, ng2
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: tcnelev, tcjb, tcjs, tciseg, kstrsplt, ncountcw, seld
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: tspltjb, tsseld, nouts, nout0, nout1, nout2
  INTEGER,       ALLOCATABLE, DIMENSION(:,:)   :: jstsplt, ncountc, tsprior
  REAL                                         :: nxtstr, nxttcd, nxtsplit, tcdfreq, tfrqtmp, tspltfreq
  REAL                                         :: sum_minfrac1, sum_minfrac2, qfrac1, qfrac2, tsconv
  REAL,          ALLOCATABLE, DIMENSION(:)     :: tctemp, tctend, tctsrt, tcklay, tspltt, volm, qwdfrac, tstend, tstsrt, nxsel,temp2
  REAL,          ALLOCATABLE, DIMENSION(:)     :: nxtssel, tstemp2, ewdsav, minfrac1, minfrac2, maxfrac1, maxfrac2, splt2t, splt2e
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: tcelev, tempcrit, qstrfrac, volmc
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: tsdepth, tsminfrac, tsminhead, tsmaxhead, tsmaxflow, estrsav
  CHARACTER(8)                                 :: tempc, tspltc
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: tcelevcon, tcyearly, tcntr, dynsel
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: tspltcntr, tsyearly, elcontspl, tsdynsel
  CHARACTER(5),  ALLOCATABLE, DIMENSION(:,:)   :: tstype
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: wd_active, share_flow
  LOGICAL,       ALLOCATABLE, DIMENSION(:,:)   :: no_flow, str_active
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: DYNSF, DYNTF
  
  !aeb    New variables for TCD operations
  INTEGER                                      :: maxouts, maxGates                     !aeb  "maxouts" dimensions outlet arrays (formerly limited to 10 outlets): maxGates dimensions maxGateNo
  INTEGER                                      :: BEGIN_HI, END_HI                      !aeb  Add start and end of gate restriction period
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: maxGateNo, minGateNo, nGates
  INTEGER,       ALLOCATABLE, DIMENSION(:,:)   :: tsminjday, tsmaxjday                  !aeb  Add min-max jday to restrict gate opening
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: tsminfrac2, jst_Volume
  LOGICAL                                      :: FORECASTING

  !aeb    New variables for checking level (gate) availability
  INTEGER                                      :: current_low_Gate, current_hi_Gate     !aeb Keep track of prior attempt to open low gate for blending (to record consecutive attempts)
  INTEGER                                      :: N_ATTEMPTS, iAttempt                  !aeb Keep track of prior attempt to open low gate for blending (to record consecutive attempts)
  INTEGER                                      :: TCD_DOWN_GATE                        
  INTEGER,                    DIMENSION(1:3670):: lowAttempt, hiAttempt                 !aeb This allocates blending frequency (= # attempts) of up to 10/day for a year
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: jst_TCDdown
  REAL                                         :: release_Volume, averaging_Interval    ! blended_Volume, not Blended_Volume
  REAL                                         :: LWR_GATE, MIN_SUBMERGENCE
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: minWSE, invertElev
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: offLimits_low, offLimits_hi
  LOGICAL                                      :: HI_GATE_RESTRICTIONS, attempt_offLimits
End Module Selective1TCD

Subroutine SelectiveInitTCD
  Use Selective1TCD; USE MAIN
  USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  IMPLICIT NONE
  
  integer      :: ifile, nj, N, JJ, jst       
  REAL         :: daytest
  CHARACTER(8) :: AID1, tsshare, FCST
  CHARACTER(1) :: CHAR1
  maxouts = 30          !aeb. Max number of structres per blending group
  maxGates = 10         !aeb. Max number of gates per blending group
  
  allocate (offLimits_low(0:maxGates), offLimits_hi(0:maxGates))
 
  offLimits_hi(:) = .FALSE.     !aeb  All levels (gates) available for hi blending at first
  offLimits_low(:) = .FALSE.    !aeb  All levels (gates) available for low blending at first.
  lowAttempt(:) = 0             !aeb  No days have been checked
  hiAttempt(:) = 0              !aeb  No days have been checked
   
  attempt_offLimits = .FALSE.   !aeb Gate combinations are OK at first
  HI_GATE_RESTRICTIONS = .FALSE.!aeb Start with no restrictions on setting hi_Gate
  iAttempt = 0                  !aeb  Initialize number of attempts to set gates
  LWR_GATE = 0.10
  TCD_DOWN_GATE = 4             !Define TCDdown as the 4th outlet of gate #4, the side gate
  current_low_Gate = 1          !aeb  Initialize low gate
  current_hi_Gate = 1           !aeb  Initialize hi gate
  
  !Debug gates attempts
  open(7777,file='Gate_attempts.csv',status='unknown')
  write(7777,*)'  JDay,',' Current,','  ,',' Attempt,','  ,'
  write(7777,*)'      ,','   hi,','  low,','   hi,','  low,'

  !Collect flows (check average of previous day's flows when forecasting)
  open(7778,file='Daily_flows.csv',status='unknown')
  write(7778,*)'    JDay,','  AvgGateQ,',' AvgNotBlendQ,','Interval(sec),'  
  
  ifile=1949
  tavg=0.0
  tavgw=0.0
  do jb=1,nbr
    if(nstr(jb) > 0)then
      ifile=ifile+1
      write (segnum,'(i0)') jb
      segnum = adjustl(segnum)
      l      = len_trim(segnum)

      IF(RESTART_IN)THEN
        OPEN  (ifile,file='str_br'//segnum(1:l)//'.csv',POSITION='APPEND')
        JDAY1=0.0
        REWIND (IFILE)
        READ   (IFILE,'(//)',END=13)
        DO WHILE (JDAY1 < JDAY)
          READ (IFILE,*,END=13) JDAY1              
        END DO
        BACKSPACE (IFILE)
13      JDAY1=0.0
      ELSE
        open  (ifile,file='str_br'//segnum(1:l)//'.csv',status='unknown')
        write (ifile,*)'Branch:,',jb,', # of structures:,',nstr(jb),' outlet temperatures'
        write (ifile,'("      JDAY,",<nstr(jb)>(6x,"T(C),"),<nstr(jb)>(3x,"Q(m3/s),"),<nstr(jb)>(4x,"ELEVCL,"))')
      ENDIF
    endif
  end do

  if(nwd > 0)then
    ifile=ifile+1
    IF(RESTART_IN)THEN
      OPEN  (ifile,file='wd_out.opt',POSITION='APPEND')
      JDAY1=0.0
      REWIND (IFILE)
      READ   (IFILE,'(/)',END=14)
      DO WHILE (JDAY1 < JDAY)
        READ (IFILE,'(F10.0)',END=14) JDAY1
      END DO
      BACKSPACE (IFILE)
14    JDAY1=0.0
    ELSE
      open  (ifile,file='wd_out.opt',status='unknown')
      write (ifile,*)'Withdrawals: # of withdrawals:',nwd,' outlet temperatures'
      write (ifile,'("      JDAY",<nwd>(6x,"T(C)"),<nwd>(3x,"Q(m3/s)"),<nwd>(4x,"ELEVCL"))')
    ENDIF
  end if

  !ZZ 8/2024, w2_selective.npt in a CSV format.
  open (NUNIT,file='w2_selective.npt',status='old')   
  READ (NUNIT,'(A)') AID1 
  IF (INDEX(AID1, "$") == 0) THEN
    write (w2err, '(A)') 'w2_selective.npt for the TCD selective withdrawal must be in a CSV format!'
  END IF
  DO J=1,2
    READ (NUNIT,*)
  END DO
  READ (NUNIT,*) AID1, tfrqtmp 
  DO J=1,2
    READ (NUNIT,*)
  END DO  
  READ (NUNIT,*) AID1, tempc, numtempc, tcdfreq
  tempc=ADJUSTR(tempc)

  IF(JDAY>TMSTRT)THEN   ! SW 8/10/2023  During restart this allows continuing output
      NXTSTR=JDAY          
      NXTTCD=JDAY   
      NXTSPLIT=JDAY
  ELSE
      NXTSTR=TMSTRT
      NXTTCD=TMSTRT   
      NXTSPLIT=TMSTRT
  ENDIF  
  
  allocate (tcnelev(numtempc),tcjb(numtempc),tcjs(numtempc),tcelev(numtempc,11),tctemp(numtempc),tctend(numtempc),tctsrt(numtempc))
  allocate (ncountc(nst,nbr),tciseg(numtempc),tcklay(numtempc),tcelevcon(numtempc))
  Allocate (tcyearly(numtempc), tcntr(numtempc))
  allocate (volm(nwb),ncountcw(nwd),qwdfrac(nwd),qstrfrac(nst,nbr),DYNSEL(numtempc),SELD(numtempc),NXSEL(numtempc),TEMP2(numtempc))
  ALLOCATE (DYNSF(numtempc))    
  DYNSF=.FALSE.
  
  ncountc=0
  DO J=1,2
    READ (NUNIT,*)
  END DO  
  do j=1,numtempc
    READ (NUNIT,*) AID1, tcntr(j),tcjb(j),tcjs(j),tcyearly(j),tctsrt(j),tctend(j),tctemp(j),tcnelev(j),(tcelev(j,n),n=1,tcnelev(j))
    tcntr(j)=ADJUSTR(tcntr(j))
    if(tcntr(j)=='      ST')then
      tcelev(j,tcnelev(j)+1)=ESTR(tcjs(j),tcjb(j))   ! always put the original elevation as the last elevation
    else
      tcelev(j,tcnelev(j)+1)=EWD(tcjs(j))            ! always put the original elevation as the last elevation
    endif
  end do 
  DO J=1,2
    READ (NUNIT,*)
  END DO  
  do j=1,numtempc
    READ (NUNIT,*) AID1, tciseg(j), tcklay(j), DYNSEL(J)
    DYNSEL(J)=ADJUSTR(DYNSEL(J))  
  end do
  DO J=1,2
    READ (NUNIT,*)
  END DO  
  do j=1,numtempc
    READ (NUNIT,*) AID1, tcelevcon(j)
    tcelevcon(j)=ADJUSTR(tcelevcon(j))   
  end do
  DO J=1,2
    READ (NUNIT,*)
  END DO  
  READ (NUNIT,*) AID1, tspltc, numtsplt, tspltfreq, tsconv, FCST, LWR_GATE
  tspltc=ADJUSTR(tspltc); FCST=ADJUSTR(FCST)
  FORECASTING = FCST == '      ON'
  if (LWR_GATE < 0.0 .or. LWR_GATE > 2.0) LWR_GATE = 0.10
    
  allocate (tsyearly(numtsplt), tstsrt(numtsplt), tstend(numtsplt), tspltjb(numtsplt), tspltt(numtsplt), nouts(numtsplt))
  allocate (jstsplt(numtsplt,maxouts), kstrsplt(maxouts), tspltcntr(numtsplt), elcontspl(numtsplt))                        
  allocate (tsdepth(numtsplt,maxouts), tstype(numtsplt,maxouts), tsminfrac(numtsplt,maxouts), tsprior(numtsplt,maxouts))   
  allocate (tsminhead(numtsplt,maxouts), tsmaxhead(numtsplt,maxouts), tsmaxflow(numtsplt,maxouts), no_flow(numtsplt,maxouts), share_flow(numtsplt))  
  allocate (tsdynsel(numtsplt), tsseld(numtsplt), nxtssel(numtsplt), tstemp2(numtsplt))
  allocate (nout0(maxouts), nout1(maxouts), nout2(maxouts), minfrac1(maxouts), maxfrac1(maxouts), minfrac2(maxouts), maxfrac2(maxouts), splt2t(maxouts), splt2e(maxouts))   
  allocate (ewdsav(nwd), wd_active(nwd), estrsav(nst,nbr), str_active(nst,nbr))
  allocate (tsminjday(numtsplt, maxouts), tsmaxjday(numtsplt, maxouts), minGateNo(numtsplt), maxGateNo(numtsplt), nGates(numtsplt), tsminfrac2(numtsplt,maxouts))     !aeb  Add min-max jday for gates to be open
  allocate (jst_Volume(numtsplt,maxouts), minWSE(numtsplt,4), invertElev(numtsplt,4), jst_TCDdown(numtsplt))   
  ALLOCATE (DYNTF(numtsplt))    
  DYNTF=.FALSE.
                       
  !Initialize variables
  jst_Volume = 0.0
  averaging_Interval = 0.0
  
  DO J=1,2
    READ (NUNIT,*)
  END DO  
  do j=1,numtsplt
    READ (NUNIT,*) AID1, tspltcntr(j), tspltjb(j), tsyearly(j), tstsrt(j), tstend(j), tspltt(j), tsdynsel(j), elcontspl(j), nouts(j), tsshare
    tspltcntr(j)=ADJUSTR(tspltcntr(j)); tsyearly(j)=ADJUSTR(tsyearly(j)); tsdynsel(j)=ADJUSTR(tsdynsel(j)); elcontspl(j)=ADJUSTR(elcontspl(j)); tsshare=ADJUSTR(tsshare)
    if (tspltc == '      ON') then
      if (nouts(j) < 2) then
        write (w2err, '(A,I0)') 'ERROR-- Less than two outlets specified for blending group ',j
        ERROR_OPEN = .TRUE.               ! will trigger the program to end when this subroutine is completed
        return
      else if (nouts(j) > maxouts) then   
        write (w2err, '(A,I0,A,I0)') 'ERROR-- More than ",maxouts," outlets specified for blending group ',j
        ERROR_OPEN = .TRUE.
        return
      end if
    end if
    share_flow(j) = tsshare == '      ON'
  end do
  DO J=1,2
    READ (NUNIT,*)
  END DO  
  do j=1,numtsplt
    READ (NUNIT,*) AID1, (jstsplt(j,n),n=1,nouts(j))
  end do
  DO J=1,2
    READ (NUNIT,*)
  END DO  
  do j=1,numtsplt
    READ (NUNIT,*) AID1, (tsdepth(j,n),n=1,nouts(j))
  end do
  DO J=1,2
    READ (NUNIT,*)
  END DO  
  do j=1,numtsplt
    READ (NUNIT,*) AID1, (tsminfrac(j,n),n=1,nouts(j))
  end do
  DO J=1,2
    READ (NUNIT,*)
  END DO  
  do j=1,numtsplt
    READ (NUNIT,*) AID1, (tsprior(j,n),n=1,nouts(j))
  end do
  DO J=1,2
    READ (NUNIT,*)
  END DO  
  do j=1,numtsplt
    READ (NUNIT,*) AID1, (tsminhead(j,n),n=1,nouts(j))
  end do
  DO J=1,2
    READ (NUNIT,*)
  END DO  
  do j=1,numtsplt
    READ (NUNIT,*) AID1, (tsmaxhead(j,n),n=1,nouts(j))
  end do
  DO J=1,2
    READ (NUNIT,*)
  END DO  
  do j=1,numtsplt
    READ (NUNIT,*) AID1, (tsmaxflow(j,n),n=1,nouts(j))
  end do   
  DO J=1,2
    READ (NUNIT,*)
  END DO  
  do j=1,numtsplt
    READ (NUNIT,*) AID1, (tsminfrac2(j,n),n=1,nouts(j))
  end do  
  DO J=1,2
    READ (NUNIT,*)
  END DO  
  do j=1,numtsplt       !aeb  Get restrictions on when gates can be open, min-max jday   
    READ (NUNIT,*) AID1, (tsminjday(j,n),n=1,nouts(j))
  end do    
  DO J=1,2
    READ (NUNIT,*)
  END DO  
  do j=1,numtsplt         
    READ (NUNIT,*) AID1, (tsmaxjday(j,n),n=1,nouts(j))
  end do 
  DO J=1,2
    READ (NUNIT,*)
  END DO
  
  ! FORCAST input parameters:
  READ (NUNIT,*) AID1, BEGIN_HI, END_HI, N_ATTEMPTS, current_low_Gate, MIN_SUBMERGENCE      
  if(N_ATTEMPTS < 0.0 .or. N_ATTEMPTS > 365) N_ATTEMPTS = 3                                 ! Default hr to set gates every day
  If(current_low_Gate < 1 .or. current_low_Gate > 4) current_low_Gate = 1                   ! If initial gate level is invalid, set to top gate
  current_hi_Gate = current_low_Gate

  ! Set gate (level) information
  maxGateNo= 0
  minGateNo= 9999
  invertElev= 9999.
  minWSE= 9999.
  jst_TCDdown= 0
  
  estrsav = estr    ! Save the original structure elevations
  ewdsav  = ewd     ! Save the original withdrawal elevations

  !***********************************************************************************************************************************
  do j=1,numtsplt
    do n=1,nouts(j)
      tstype(j,n) = "FIXED"
      if (tsdepth(j,n)   > 0.0) tstype(j,n) = "FLOAT"
      if (tsminfrac(j,n) > 1.0) tsminfrac(j,n) = 1.0    ! remove unrealistic input value
      if (tsminhead(j,n) < 0.0) tsminhead(j,n) = 0.0    
      if (tsmaxhead(j,n) < 0.0) tsmaxhead(j,n) = 0.0    
      if (tsmaxflow(j,n) < 0.0) tsmaxflow(j,n) = 0.0    
    end do
  end do  
  
  do j=1,numtsplt
    jb = tspltjb(j)
    do n=1,nouts(j)
      if (tsprior(j,n) == 5) then   ! Add TCDdown to side gate and record structure number
          tsprior(j,n) = TCD_DOWN_GATE
          jst_TCDdown(j) = jstsplt(j,n)
      end if
      if (tsprior(j,n) > 0) then    ! Get lowest and highest gate number for each blending group, and invert elev for each gate (level)                     
          jst = jstsplt(j,n)
          if (tsprior(j,n) > maxGateNo(j)) maxGateNo(j) = tsprior(j,n)            ! Get highest gate number for each blending group
          if (tsprior(j,n) < minGateNo(j)) minGateNo(j) = tsprior(j,n)            ! Get lowest gate number for each blending group
          if (invertElev(j, tsprior(j,n)) > estr(jst,jb)) invertElev(j, tsprior(j,n)) = estr(jst,jb)    !Get elevation of lowest outlet at gate level. 
      end if
      if (tsminfrac2(j,n) > 1.0) tsminfrac2(j,n) = 1.0  ! remove unrealistic input value
    end do
    nGates(j) = maxGateNo(j) - minGateNo(j) + 1         ! Calculate # of gates; gates should be contiguous. NOT CURRENTLY USED
   
    do n= 1,4
        minWSE(j,n)= invertElev(j,n) + MIN_SUBMERGENCE  ! Set minimum WSE required to allow each level to open alone
    end do
  end do  
  !***********************************************************************************************************************************
  
  if (tsconv <= 0.0) tsconv = 0.005   ! constrain the convergence criterion to be > 0.0 and <= 0.1
  if (tsconv >  0.1) tsconv = 0.1

  DO J=1,2
    READ (NUNIT,*)
  END DO
  READ (NUNIT,*) AID1, tempn
  allocate (tempcrit(nwb,tempn),volmc(nwb,tempn)) 
  DO J=1,2
    READ (NUNIT,*)
  END DO
  do j=1,tempn
    READ (NUNIT,*) AID1, (tempcrit(jw,j), jw=1,nwb) 
  end do
  close (NUNIT)

  do jw=1,nwb
    ifile=ifile+1
    write (segnum,'(i0)') jw
    segnum = adjustl(segnum)
    l      = len_trim(segnum)
    IF(RESTART_IN)THEN
      OPEN  (ifile,file='Volume_wb'//segnum(1:l)//'.opt',POSITION='APPEND')
      JDAY1=0.0
      REWIND (IFILE)
      READ   (IFILE,'(/)',END=15)
      DO WHILE (JDAY1 < JDAY)
        READ (IFILE,'(F10.0)',END=15) JDAY1
      END DO
      BACKSPACE (IFILE)
15    JDAY1=0.0
    ELSE
      open (ifile,file='Volume_wb'//segnum(1:l)//'.opt',status='unknown')
      write(ifile,'("    jday      Volume     ",<tempn>("Volcrit     "))')
    END IF
  end do

  if (tempc == '      ON') then
    do j=1,numtempc
      if (tcyearly(j) == '     OFF') then
        daytest=jday
      else
        daytest=real(jdayg)+jday-int(jday)
      end if
      if (daytest >= tctsrt(j) .and. daytest < tctend(j)) then

      ! initializing structure elevation
        if (tcntr(j) == '      ST') then
          jb = tcjb(j)                                ! set branch index
          js = tcjs(j)                                ! set structure index
          do nj=1,tcnelev(j)                          ! making sure that structure is below water surface
            if (tcelev(j,nj) < elws(ds(jb))) then
              ncountc(js,jb)=nj
              exit
            end if
          end do

      ! initializing withdrawal elevation
        else if (tcntr(j) == '      WD') then
          jwd = tcjs(j)                               ! set withdrawal index
          do nj=1,tcnelev(j)                          ! making sure that structure is below water surface
            if (tcelev(j,nj) < elws(iwd(jwd))) then
              ncountcw(jwd)=nj
              exit
            end if
          end do
        end if
      end if

    ! Open dynamic selective withdrawal files
      IF (DYNSEL(J) == '      ON') THEN
        WRITE (SEGNUM,'(I0)') J
        SEGNUM  = ADJUSTL(SEGNUM)
        L       = LEN_TRIM(SEGNUM)
        SELD(J) = 1009+J
        OPEN (SELD(J),FILE='dynselective'//SEGNUM(1:L)//'.npt',STATUS='OLD')      
        READ(SELD(J),'(A1)') CHAR1    
        IF(CHAR1=='$') DYNSF(J)=.TRUE.
        IF(DYNSF(J))THEN
          READ (SELD(J),'(/)')
          READ (SELD(J),*) NXSEL(J),TEMP2(J)
          tctemp(J)=TEMP2(J)
          READ (SELD(J),*) NXSEL(J),TEMP2(J)       
        ELSE
          READ (SELD(J),'(//1000F8.0)') NXSEL(J),TEMP2(J)
          tctemp(J)=TEMP2(J)
          READ (SELD(J),'(1000F8.0)') NXSEL(J),TEMP2(J)
        ENDIF       
      end if
    end do
  end if

! Open dynamic temperature target files for blending outlets
  if (tspltc == '      ON') then
    do j=1,numtsplt
      if (tsdynsel(j) == '      ON') then
        write (segnum,'(i0)') j
        segnum    = adjustl(segnum)
        L         = len_trim(segnum)
        tsseld(j) = 1009+numtempc+j
        open (tsseld(j),file='dynsplit_selective'//segnum(1:L)//'.npt',status='old')     
        READ(tsseld(j),'(A1)') CHAR1    
        IF(CHAR1=='$') DYNTF(J)=.TRUE.
        IF(DYNTF(J))THEN
          READ (tsseld(j),'(/)')
          READ (tsseld(j),*) nxtssel(j), tstemp2(j)
          tctemp(J)=TEMP2(J)
          READ (tsseld(j),*) nxtssel(j), tstemp2(j)       
        ELSE
          READ (tsseld(j),'(//2F8.0)') nxtssel(j), tstemp2(j)
          tspltt(j) = tstemp2(j)
          READ (tsseld(j),'(2F8.0)') nxtssel(j), tstemp2(j)
        ENDIF    
      end if
    end do
  end if

! Test to see if the user specified inconsistent inputs. If so, stop with an error message.
  if (tspltc == '      ON') then
    if (tspltfreq <= 0.0) then                                                                                       
      write (w2err, '(A)') 'w2_selective.npt TCD ERROR-- Update frequency for temperature blending must be greater than zero.'
      ERROR_OPEN = .TRUE.                                                                                             
    end if     
    do j=1,numtsplt
      do n=1,nouts(j)-1
        do nj=n+1,nouts(j)
          if (jstsplt(j,n) == jstsplt(j,nj)) then
            write (w2err, '(A,I0)') 'w2_selective.npt TCD ERROR-- Duplicate split outlet numbers in group ', j
            ERROR_OPEN = .TRUE.     ! will trigger the program to end when this subroutine is completed
          end if
        end do
      end do
    end do
    do j=1,numtsplt-1
      do jj=j+1,numtsplt
        if ((tstsrt(jj) >= tstsrt(j) .and. tstsrt(jj) <  tstend(j)) .or.  &
            (tstend(jj) >  tstsrt(j) .and. tstend(jj) <= tstend(j))) then
          if (tspltcntr(j) .eq. tspltcntr(jj) .and. (tspltcntr(j) .eq. '      WD' .or. tspltjb(j) == tspltjb(jj))) then
            do n=1,nouts(j)
              do nj=1,nouts(jj)
                if (jstsplt(j,n) == jstsplt(jj,nj)) then
                  write (w2err, '(A,I0,A)') 'w2_selective.npt TCD ERROR-- Split outlet number ', jstsplt(j,n), ' used in more than one group at a time.'
                  ERROR_OPEN = .TRUE.       ! will trigger the program to end when this subroutine is completed
                end if
              end do
            end do
          end if
        end if
      end do
    end do
    do j=1,numtsplt
      do n=1,nouts(j)
        if (tsprior(j,n) < -1) then
          tsprior(j,n) = -1           
          write (wrn, '(A,I0,A,I0,A)') 'w2_selective TCD: WARNING-- Priority input for outlet ', jstsplt(j,n), &
                                       ' in group ', j, ' reassigned to -1.'                                          
          WARNING_OPEN = .TRUE.                                                                                     
        end if               
        if (tsminhead(j,n) > 0.0 .and. tsmaxhead(j,n) > 0.0 .and. tsminhead(j,n) > tsmaxhead(j,n)) then
          write (wrn, '(A,I0,A,I0,A)') 'w2_selective TCD: WARNING-- Minimum and maximum head constraints for outlet ', jstsplt(j,n), ' in group ',  &
                                        j, ' are such that the outlet cannot ever be used.'
          WARNING_OPEN = .TRUE.
        end if
        if (tsdepth(j,n) > 0.0 .and. tsminhead(j,n) > 0.0 .and. tsdepth(j,n) < tsminhead(j,n)) then
          write (wrn, '(A,I0,A,I0,A)') 'w2_selective TCD: WARNING-- Depth of floating outlet ', jstsplt(j,n), ' in group ', j,  &
              ' is shallower than the minimum head constraint.  To honor the head constraint, no flow is possible for that outlet.'
          WARNING_OPEN = .TRUE.
        end if
      end do
    end do

    if (tempc == '      ON') then
      do j=1,numtsplt
        do jj=1,numtempc
          if ((tctsrt(jj) >= tstsrt(j) .and. tctsrt(jj) <  tstend(j)) .or.  &
              (tctend(jj) >  tstsrt(j) .and. tctend(jj) <= tstend(j))) then
            if (tspltcntr(j) .eq. tcntr(jj) .and. (tspltcntr(j) .eq. '      WD' .or. tspltjb(j) == tcjb(jj))) then
              do n=1,nouts(j)
                if (jstsplt(j,n) == tcjs(jj)) then
                  write (w2err, '(A,I0,A)') 'w2_selective TCD ERROR-- Outlet number ',tcjs(jj),' used in tower and blending group at same time.'
                  ERROR_OPEN = .TRUE.       ! will trigger the program to end when this subroutine is completed
                end if
              end do
            end if
          end if
        end do
      end do
    end if
  end if

  if (tempc == '      ON') then
    do j=1,numtempc-1
      do jj=j+1,numtempc
        if ((tctsrt(jj) >= tctsrt(j) .and. tctsrt(jj) <  tctend(j)) .or.  &
            (tctend(jj) >  tctsrt(j) .and. tctend(jj) <= tctend(j))) then
          if (tcntr(j) .eq. tcntr(jj) .and. (tcntr(j) .eq. '      WD' .or. tcjb(j) == tcjb(jj))) then
            if (tcjs(j) == tcjs(jj)) then
              write (w2err, '(A,I0,A)') 'w2_selective TCD: ERROR-- Tower outlet number ', tcjs(j), ' used more than once for overlapping dates.'
              ERROR_OPEN = .TRUE.       ! will trigger the program to end when this subroutine is completed
            end if
          end if
        end if
      end do
    end do
  end if

  return
End subroutine SelectiveInitTCD

!***********************************************************************************************************************************
!**                                                   S E L E C T I V E  T C D                                                    **
!***********************************************************************************************************************************
Subroutine SelectiveTCD
  Use Selective1TCD; USE MAIN
  USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  IMPLICIT NONE

  integer :: jj, jst, n, nj, num_noflow, ng0, prior1, prior2, ng1max, ng1min, num_left
  integer :: j2hi, j2lo, j2max, j2min, j2pref
  real    :: qall, elr, wsel, q_notblended, sum_minfrac0, sum_maxfrac1, sum_maxfrac2, sumfrac
  real    :: maxelev, minelev, blendfrac, excess_frac, addfrac, maxtemp, mintemp
  real    :: lastfrac, lastfrac2, ttarg, sumtemp, etemp, etemp1, etemp2, sumelev, elev1, elev2

  INTEGER :: JJW, KK, KS, IFILE, KSTR
  REAL    :: daytest, TCOMP, TEMPEST, TMOD

  ! Variables for bisection-type gate search
  REAL     ::	f_high, f_low, high_frac, tm1
  INTEGER  ::	max_iter, hi_Gate, low_Gate, low_Str, jStr_hi, jStr_low  !, blend_Gate

  ! New variables for blending
  REAL     :: orig_qstr(NST,NBR)  
  REAL     :: q_blended, q_Gates, gate_fraction, sum_wgt_Tw
  REAL     ::	q_hi, q_low, Tw_hi, Tw_low, q_temp
  REAL     :: theHour
  INTEGER  ::	iDay, min_BlendingGate, max_BlendingGate, lowBlendGate
                 
  ! Variables to put gate variables in arrays
  INTEGER  ::	iGate, iOut
  INTEGER  ::	nOutlets(maxGates), jstOut(maxGates,maxouts)
  REAL     ::	min_flow(maxouts), min_flow1(maxouts), min_flow2(maxouts), Tw_out(maxouts)
  REAL     ::	etemp_minQ(maxouts), etemp_minQ_1(maxouts), etemp_minQ_2(maxouts)
  
  ! Glossary
  !q_Gate is flow through gates (q_Gates = qall - q_notblended)
  !q_blended is actual blended flow (q_blended = q_Gates - gate_in_flows)
  !gate_fraction is sum total blended fraction from previous blending
  !sum_wgt_Tw is sum of flow-weighted Tw used in adjusting Tw target
  !nOutlets(iGate) holds number of structures associated with gate
  !jstsplit(j,jj) holds the place number of the structure in the list of all structures for this simulation (USGS variable)
  !jstOut(iGate,iOut) holds the jstsplit value for structure, iOut, in gate, iGate
  !min_flow(maxouts) holds the minimum flow for a structure, converted from input file or calculated from min fraction
  !Tw_out is the estimated release Tw assuming all blended flow through specified outlet
  !maxGateNo(j) is the highest gate number listed for blending group, j
  !nAttempts(maxGates) counts the number of consecutive attempts to use a gate; a condition for making gate available

  str_active = .FALSE.
  wd_active  = .FALSE.
  j2pref     = 1
  NWD = 0			    !	Ignore withdrawals
  orig_qstr = 0.0 ! Initialize array to temporarily store qstr values when forecast blending using interval averages

  ! Update some date variables and determine which outlets are being actively blended or adjusted.
  if (tspltc == '      ON') then
    do j=1,numtsplt
        if (tsyearly(j) == '     OFF') then
            daytest = jday
        else
            daytest = real(jdayg) + jday - int(jday)
        end if
        if (nxtsplit > tstsrt(j) .and. daytest <= tstsrt(j)) then
            nxtsplit = tstsrt(j)
        end if
        if (daytest >= tstsrt(j) .and. daytest < tstend(j)) then
            do jj=1,nouts(j)
                if (tspltcntr(j) == '      ST') then
                    str_active(jstsplt(j,jj),tspltjb(j)) = .TRUE.
                else if (tspltcntr(j) == '      WD') then
                    !aeb	Ignore withdrawals
                end if
            end do
            
            !ZZ, 10/24
            if (FORECASTING) then !aeb. Accumulate release volumes and averaging interval to average flow for forecast blending calculations
                do jj=1,nouts(j)
                    jst = jstsplt(j,jj)
                    jst_Volume(j,jj) = jst_Volume(j,jj) + qstr(jst,tspltjb(j)) * DLT
                end do
                averaging_Interval = averaging_Interval + DLT
            end if                       
        end if
    end do
  end if

  ! Reset elevations of outlets back to original values outside of control periods
  do jst=1,nst
      do jb=1,nbr
          if (.not. str_active(jst,jb)) estr(jst,jb) = estrsav(jst,jb)
      end do
  end do

  ! Check to see if it's time to update temperature targets and flow fractions for blended groups.
  if (tspltc=='      ON' .and. jday .ge. nxtsplit) then     ! At end of this section: nxtsplit = nxtsplit + tspltfreq

    ! Update the temperature targets
    do j=1,numtsplt
        if (tsdynsel(j) == '      ON') then
            do while (jday >= nxtssel(j))
                tspltt(j) = tstemp2(j)
                IF(DYNTF(J))THEN
                  read (tsseld(j),*) nxtssel(j), tstemp2(j)
                ELSE
                  read (tsseld(j),'(2F8.0)') nxtssel(j), tstemp2(j)
                END IF
            end do
        end if
    end do

    do j=1,numtsplt                                                            ! Cycle through all blending groups.
        qall    = 0.0                                                          ! sum up all the flows
        sumfrac = 0.0                                                          ! sum of flow fraction multipliers
        do jj=1,nouts(j)
            if (tspltcntr(j) == '      ST') then
                qall    = qall    + qstr(jstsplt(j,jj),tspltjb(j))              !aeb   This is for test below.
                sumfrac = sumfrac + qstrfrac(jstsplt(j,jj),tspltjb(j))          !aeb   qstrfrac() is either uninitialized or carried over from previous time step
            end if
        end do
        if (tsyearly(j) == '     OFF') then
            daytest = jday
        else
            daytest = real(jdayg) + jday - int(jday)
        end if

        ! Do blending calculations if date is in correct window or the window was just entered.
        ! If flows are zero and flow fractions are already initialized, then leave everything alone.
        if (daytest >= tstsrt(j) .and. daytest < tstend(j) .and.  &
            (qall > 0.0 .or. sumfrac < 0.001 .or. daytest < tstsrt(j) + tspltfreq)) then
  
          if (tspltcntr(j) == '      ST') then                                 ! Do structures only
            jb = tspltjb(j)                                                    ! set branch index
            do jw=1,nwb                                                        ! set waterbody index
                if (jb >= bs(jw) .and. jb <= be(jw)) exit
            end do
            no_flow(j,:) = .FALSE.
            num_noflow   = 0
            do jj=1,nouts(j)                                            
                jst          = jstsplt(j,jj)
                elr          = sina(jb) * dlx(ds(jb)) * 0.5                      
                wsel         = elws(ds(jb)) - elr                                ! compute water-surface elevation
                estr(jst,jb) = estrsav(jst,jb)                                   ! reset outlet elevation to original
                if (tstype(j,jj) == "FLOAT") then
                    estr(jst,jb) = wsel - tsdepth(j,jj)
                else if (estr(jst,jb) > wsel) then
                    if (elcontspl(j) == '     OFF') then
                        no_flow(j,jj) = .TRUE.                                    ! no flow-- high and dry
                        num_noflow    = num_noflow + 1
                    else
                        estr(jst,jb) = wsel                                       ! poor man's floating outlet
                    end if
                end if
                if (.not. no_flow(j,jj) .and. tsminhead(j,jj) > 0.0 .and. wsel - estr(jst,jb) < tsminhead(j,jj)) then
                    no_flow(j,jj) = .TRUE.                                         ! minimum head criterion not met -- no flow
                    num_noflow    = num_noflow + 1
                end if
                if (.not. no_flow(j,jj) .and. tsmaxhead(j,jj) > 0.0 .and. wsel - estr(jst,jb) > tsmaxhead(j,jj)) then
                    no_flow(j,jj) = .TRUE.                                         ! maximum head criterion exceeded -- no flow
                    num_noflow    = num_noflow + 1
                end if
                
                !	Check to see if we're within gate restriction period
                if(.not. no_flow(j,jj) .and. tsminjday(j,jj) > 0 .and. jday < tsminjday(j,jj)) then
                    no_flow(j,jj) = .TRUE.                                          !aeb	Too soon to open gate -- no flow
                    num_noflow    = num_noflow + 1
                end if
                if(.not. no_flow(j,jj) .and. tsmaxjday(j,jj) > 0 .and. jday > tsmaxjday(j,jj)) then
                    no_flow(j,jj) = .TRUE.                                          !aeb	Too late to open gate -- no flow
                    num_noflow    = num_noflow + 1
                end if

                do k=ktwb(jw),kb(ds(jb))
                    if (el(k,ds(jb))-elr < estr(jst,jb)) exit
                end do
                kstrsplt(jj)     = min(k-1,kb(ds(jb)))
                qstrfrac(jst,jb) = 0.0                                              ! initialize flow fractions
            end do

            !ZZ, 10/24
            if (FORECASTING) then              
                ! Keep track of original structure flows to re-set after establishing blending fractions from interval averages
                orig_qstr = qstr

                ! Set flows to interval averages for blending calculations
                qall = 0.0
                do jj=1,nouts(j)
                    jst = jstsplt(j,jj)
                    qstr(jst,jb) = jst_Volume(j,jj)/averaging_Interval
                    qall = qall + qstr(jst,jb)
                end do
            end if            
             
            ! Sum and set flow fractions for nonblended flows.  Outlets with a priority of -1 get used, but are not blended
            ng0 = 0
            q_notblended = 0.0
            do jj=1,nouts(j)
                jst = jstsplt(j,jj)
                if (.not. no_flow(j,jj) .and. tsprior(j,jj) == -1) then
                    ng0 = ng0 + 1
                    nout0(ng0) = jj
                    if (qstr(jst,jb) > tsmaxflow(j,jj) .and. tsmaxflow(j,jj) > 0.0) then  !aeb  Shouldnt we just set qstr = tsmaxflow?
                        q_notblended = q_notblended + tsmaxflow(j,jj)                     !aeb  qall should be reduced also
                        qstrfrac(jst,jb) = tsmaxflow(j,jj) / qall
                    else if (qall > 0.0) then
                        q_notblended = q_notblended + qstr(jst,jb)
                        qstrfrac(jst,jb) = qstr(jst,jb) / qall
                    end if
                end if
            end do  
            
            q_Gates = qall - q_notblended		!aeb	Total flow through gates
            
            !***********************************************************************************************************************************
            write(7778,'(4(f10.3,","))') jday, q_Gates, q_notblended, averaging_Interval
            ! Re-initialize release volumes and interval time
            jst_Volume = 0.0
            averaging_Interval = 0.0            
                        
            ! Begin the blending decisions.            
            ! No usable outlets.  All flow fractions remain at zero.
            if (nouts(j) == num_noflow) then
                write (wrn,'(A,I0,A,F0.3)') 'Warning-- All outlets dry or unusable for group ', j, ' at day ', jday
                WARNING_OPEN = .TRUE.

            ! Only nonblended outlets.
            else if (nouts(j) == ng0) then
                write (wrn,'(A,I0,A,F0.3)') 'Warning-- Only nonblended outlets present in group ', j, ' at day ', jday
                WARNING_OPEN = .TRUE.
            else    !aeb	Blend gates 
              
                !***********************************************************************************************************************************
                !ZZ 10/24 
                ! Initialize my arrays 
                nOutlets=0
                jstOut=0
                min_flow=0.0
                min_flow1=0.0
                min_flow2=0.0

                ! Convert minfracs to min flows.  These are only applied to blended flows
                do jj=1,nouts(j)
                    jst = jstsplt(j,jj)
                    if (tsminfrac(j,jj) >= 0.0) then
                        min_flow1(jst) = tsminfrac(j,jj)*q_Gates
                    else
                        min_flow1(jst) = -tsminfrac(j,jj)           !tsminfrac is given as a negative flow
                    end if
                    if (tsminfrac2(j,jj) >= 0.0) then
                        min_flow2(jst) = tsminfrac2(j,jj)*q_Gates
                    else
                        min_flow2(jst) = -tsminfrac2(j,jj)           !tsminfrac is given as a negative flow
                    end if
                end do
            
                ! Use priorities to identify outlets associated with gates.  Do this here instead of at INIT because "no_flow" may have changed   
                max_BlendingGate = 0
                min_BlendingGate = maxGateNo(j)
                do jj=1,nouts(j)           
                    if (.not. no_flow(j,jj) .and. tsprior(j,jj) > 0) then
                        iGate = tsprior(j,jj)
                        nOutlets(iGate) = nOutlets(iGate) + 1
                        jstOut(iGate,nOutlets(iGate)) = jstsplt(j,jj)
                        if (tsprior(j,jj) > max_BlendingGate) max_BlendingGate = tsprior(j,jj)    !aeb Get highest gate number for this blending group
                        if (tsprior(j,jj) < min_BlendingGate) min_BlendingGate = tsprior(j,jj)    !aeb Get lowest gate number for this blending group
                    endif                
                end do 

                ! Check current gate(s) in case they are no longer available
                if(nOutlets(current_low_Gate) == 0) then
                    do kk = min_BlendingGate, max_BlendingGate
                        if(nOutlets(kk) > 0) then
                            current_low_Gate = kk                   !aeb  Set highest gate with available outlets as low gate
                            current_hi_Gate = current_low_Gate      !aeb  Set hi gate equal to low gate
                            exit
                        end if
                    end do
                end if
                if(nOutlets(current_hi_Gate) == 0) current_hi_Gate = current_low_Gate   !aeb  Set hi gate equal to low gate              
                !***********************************************************************************************************************************
                           
                ! Set parameters for downstream_withdrawal_estimate
                id = ds(jb)                   ! needed for downstream_withdrawal_estimate
                kt = ktwb(jw)                 ! needed for downstream_withdrawal_estimate

                ! Get weighted blend of all nonblended release temperatures and adjust Tw target
                ttarg = tspltt(j)
                if (q_notblended > 0.0) then
                    sum_wgt_Tw = 0.0
                    do n=1,ng0
                        jj = nout0(n)
                        jst = jstsplt(j,jj)
                        if (qstr(jst,jb) > 0.0) then
                            call downstream_withdrawal_estimate(jst,etemp,estr(jst,jb))
                        else
                            etemp = 0.0             ! uses temperature at outlet elevation if no flow (typ.): etemp = t2(kstrsplt(nout0(n)),ds(jb))
                        end if
                        sum_wgt_Tw = sum_wgt_Tw + qstr(jst,jb) * etemp
                    end do
                    ttarg = (ttarg * qall - sum_wgt_Tw) / q_Gates  ! New temperature target for blended releases
                end if
                
!*********************************** Find gates that bracket target **********************************
                min_flow = 0.0
                etemp_minQ = 0.0
                etemp_minQ_1 = 0.0
                etemp_minQ_2 = 0.0
                Tw_out=0.0
                
                low_Gate = 0
                low_Str = 0
                hi_Gate = 0
                
                ! Initialize all gate flows to zero. (This removes carry-over flows in outlets that now have no flow.)
                do jj=1,nouts(j)
                    if (tsprior(j,jj) /= -1) then
                        jst = jstsplt(j,jj)
                        qstr(jst,jb) = 0.0 
                    endif                    
                end do
                                
               !*** Get Tw_out associated with min flows at each structure **********************************
                do iGate = min_BlendingGate, max_BlendingGate               
                    ! Put min flow1 into every outlet and get Tw estimate for each.  
                    do iOut = 1, nOutlets(iGate)    
                        jst = jstOut(iGate, iOut)
                        qstr(jst,jb) = min_flow1(jst)
                        if (qstr(jst,jb) > 0.0) then
                            call downstream_withdrawal_estimate(jst,etemp,estr(jst,jb))       !aeb Uses qstr(flow) not qstrfrac
                        else
                            etemp = 0.0         
                        end if 
                        etemp_minQ_1(jst) = etemp
                    end do
                        
                    if(min_BlendingGate /= max_BlendingGate) then          ! Might blend 2 gates
                        ! Put min flow2 into every outlet and get Tw estimate for each. For use when iGate > 1 and iOut == 1
                        do iOut = 1, nOutlets(iGate)
                            jst = jstOut(iGate, iOut)
                            qstr(jst,jb) = min_flow2(jst)
                            if (qstr(jst,jb) > 0.0) then
                                call downstream_withdrawal_estimate(jst,etemp,estr(jst,jb))       !aeb Uses qstr(flow) not qstrfrac
                            else
                                etemp = 0.0         
                            end if 
                            etemp_minQ_2(jst) = etemp
                        end do
                    end if
                end do

!*********************************** If Forecasting, establish gates and outlets for blending (at tspltfreq frequency) *************
                if (FORECASTING) then 
                    iDay = INT(JDay)                                       
                        ! Find low gate and low opening to blend            
                        do iGate = min_BlendingGate, max_BlendingGate
                            do iOut = 1, nOutlets(iGate)    ! Get weighted Tw for each opening in this gate, assuming min flows at all other openings.
                                sum_wgt_Tw = 0.0            ! Reset for each opening
                                q_blended = q_Gates         ! Reset for each opening
                                
                                if (min_BlendingGate == max_BlendingGate) then  ! Using one gate, so use min flow1 values
                                    min_flow = min_flow1
                                    etemp_minQ = etemp_minQ_1                           
                                else if(iGate == 1 .or. iOut > 1) then  ! Using one gate, so use min flow1 values
                                    min_flow = min_flow1
                                    etemp_minQ = etemp_minQ_1                               
                                else if(nOutlets(iGate-1) == 0) then    ! No outlets available in gate above.  Use min flow1 values
                                    min_flow = min_flow1
                                    etemp_minQ = etemp_minQ_1                                
                                else                              ! Special case: (2 gates .and. iGate > 1 .and. iOut == 1 .and. nOutlets(iGate-1) > 0)
                                    min_flow = min_flow2          ! If top of gate and there is a gate above, account for gate above used in blending
                                    etemp_minQ = etemp_minQ_2     ! Using two gates, so use min flow2 values                                   
                                    do kk = 1, nOutlets(iGate-1)  ! Remove min flows of gate above from blended flow and sum Tw effect 
                                        jst = jstOut(iGate-1, kk)
                                        q_blended = q_blended - min_flow(jst)
                                        sum_wgt_Tw = sum_wgt_Tw + etemp_minQ(jst)*min_flow(jst)
                                    end do
                                end if
                            
                                ! Remove min flows at this gate from blended flow
                                do kk = 1, nOutlets(iGate)
                                    jst = jstOut(iGate, kk)
                                    q_blended = q_blended - min_flow(jst)
                                end do
                                ! Blended flow is now reduced by sum of min flows at all openings and Tw_out from second gate (if there is one) is accounted for
                            
                                ! Now, find low gate and structure.  Calculate weighted Tw for each opening of this gate and test against Tw target
                                do kk = 1, nOutlets(iGate)                                
                                    ! Set flow at trial opening, iOut, to (blended + min) and estimate Tw.  Use min flow etemps for all other openings
                                    jst = jstOut(iGate, kk)
                                    if(kk == iOut) then
                                        qstr(jst,jb) = q_blended + min_flow(jst) 
                                        if (qstr(jst,jb) > 0.0) then
                                            call downstream_withdrawal_estimate(jst,etemp,estr(jst,jb))  ! Uses qstr(flow) not qstrfrac
                                        else
                                            etemp = 0.0         
                                        end if   
                                        sum_wgt_Tw = sum_wgt_Tw + etemp * qstr(jst,jb)
                                    else
                                        qstr(jst,jb) = min_flow(jst)    ! Set flow to min for debugging clarity
                                        sum_wgt_Tw = sum_wgt_Tw + etemp_minQ(jst) * min_flow(jst)
                                    end if
                                end do
                            
                                jst = jstOut(iGate, iOut)
                                Tw_out(jst) = sum_wgt_Tw / q_Gates      ! Weighted temperature of outlet, iOut, in iGate assuming it gets all but min flows
                                                                
                                ! If Tw_target is greater than Tw_out then we can meet target.  LWR_GATE raises target to bias to higher gates.
                                if(Tw_out(jst) <= (ttarg + LWR_GATE)) then 
                                    low_Gate = iGate
                                    low_Str = iOut
                                    exit
                                end if
                            end do
                            if(low_Str > 0) exit                        ! Opening found; etemp_minQ and min_flow are set for each opening
                        end do                        

                        ! Set gates
                        if(low_Str == 0) then                           ! No opening found. Target Tw is less than any opening. Use lowest available gate.  
                            low_Gate = max_BlendingGate
                        else if(low_Str == 1 .and. low_Gate > 1) then   ! Top of gate. Set higher gate if available for blending  
                            if(nOutlets(low_Gate - 1) > 0) hi_Gate = low_Gate - 1                           
                        end if  
                        if (hi_Gate == 0)  hi_Gate = low_Gate           ! No hi gate.  Blending will be done within a gate
                                                   
                        ! Write current and attempted gate combinations
                        write(7777,'(f8.3,",",4(I5,","))') jday, current_hi_Gate, current_low_Gate, hi_Gate, low_Gate

                        ! Record attempt
                        iAttempt = iAttempt + 1
                        lowAttempt(iAttempt) = low_Gate
                        hiAttempt(iAttempt) = hi_Gate
                        attempt_offLimits = .FALSE.
                                                
                        ! Check for consecutive attempts
                        if(iAttempt >= N_ATTEMPTS) then                           
                            do i = 1, N_ATTEMPTS-1
                                if(low_Gate /= lowAttempt(iAttempt-i) .or. hi_Gate /= hiAttempt(iAttempt-i)) then
                                    attempt_offLimits = .TRUE.
                                    exit
                                end if
                            end do
                        else
                            attempt_offLimits = .TRUE.
                        end if
                        
                        ! Check for off-limits gates      
                        if(offLimits_low(low_Gate) .or. offLimits_hi(hi_Gate)) attempt_offLimits = .TRUE.
                                                
                        ! Set current gate combination
                        if(attempt_offLimits) then     
                            ! Don't change current gates but ensure that current low gate still has outlets available
                            if(nOutlets(current_low_Gate) == 0) then
                                current_low_Gate = min(current_low_Gate + 1, max_BlendingGate)     ! This ensures that following "do loop" works (e.g., if max_BlendingGate < current_low_Gate + 1)
                                do kk = current_low_Gate, max_BlendingGate  ! Check for outlets
                                    if(nOutlets(kk) > 0) then
                                        current_low_Gate = kk               ! Set highest gate with outlets as low gate
                                        exit
                                    end if
                                end do
                            end if
                            
                            ! Ensure that current hi gate has outlets available
                            if(nOutlets(current_hi_Gate) == 0) current_hi_Gate = current_low_Gate
                        else
                            current_low_Gate = low_Gate
                            current_hi_Gate = hi_Gate
                        end if
                       
                        ! If one level open, check for min submergence
                        if(current_low_Gate == current_hi_Gate) then
                            if(wsel < minWSE(j,current_hi_Gate)) current_low_Gate= min(current_hi_Gate + 1, max_BlendingGate)
                        endif
                        
                        !*** Set off-limits gates based on final gate combination     
                        if(jday > BEGIN_HI .and. jday < END_HI) then                        ! Check to see if we're within gate restriction period
                            HI_GATE_RESTRICTIONS = .TRUE.
                            if (current_low_Gate > 1) then                                  ! For all gates below top gate:
                                offLimits_low(current_low_Gate-1) = .TRUE.                  ! Never allow gate above to be used exclusively
                                if (current_hi_Gate == current_low_Gate) offLimits_hi(current_low_Gate-1) = .TRUE.  !If only low gate selected, set gate above off limits for future blending
                            end if
                        else if (HI_GATE_RESTRICTIONS) then                                 ! On END_HI, reset all gates for this blending period to OK
                            do iGate = minGateNo(j), maxGateNo(j)
                                offLimits_low(iGate) = .FALSE.
                                offLimits_hi(iGate) = .FALSE.
                            end do
                            HI_GATE_RESTRICTIONS = .FALSE.
                        end if                                                
                end if         ! End of find gates that bracket target 

!*********************************** Set gates for blending (every tspltfreq interval) **********************************  
                if (FORECASTING) then 
                    low_Gate = current_low_Gate
                    hi_Gate = current_hi_Gate
                else                            
                    if(min_BlendingGate == max_BlendingGate) then
                        low_Gate = max_BlendingGate
                        hi_Gate = low_Gate
                    else
                        low_Gate = max_BlendingGate      ! Should only be 2 gates but, if more, always select bottom 2
                        hi_Gate  = low_Gate - 1          ! If here, min_BlendingGate /= max_BlendingGate and therefore low_Gate > 1
                        if(nOutlets(hi_Gate) == 0) hi_Gate  = low_Gate
                    end if
                end if
                
!*********************************** Given gates, set min flows and get Tw associated with min flows at each opening **********************************                        
                !*** Get Tw_out associated with min flows          
                etemp_minQ = 0.0
                q_blended = q_Gates
            
                if(hi_Gate == low_Gate) then          ! Blend within 1 gate
                    min_flow = min_flow1
                    etemp_minQ = etemp_minQ_1         !ZZ 10/24
                else
                    min_flow = min_flow2
                    etemp_minQ = etemp_minQ_2
                end if
            
                do iGate = hi_Gate, low_Gate          ! Put min flow into every outlet, get Tw estimate for each, reduce blended flow by min flows.            
                    do iOut = 1, nOutlets(iGate)    
                        jst = jstOut(iGate, iOut)
                        qstr(jst,jb) = min_flow(jst)
                        q_blended = q_blended - qstr(jst,jb)                  
                    end do                        
                end do
                ! Blended flow is now reduced by sum of min flows at all openings

!*********************************** Find low gate and low opening to blend **********************************
                lowBlendGate = 0     !ZZ 10/24
                low_Str = 0
                                    
                do iGate = hi_Gate, low_Gate            ! Cyle through all openings to find low gate and structure.           
                    do iOut = 1, nOutlets(iGate)        ! Get weighted Tw for each opening in this gate, assuming min flows at all other openings and test against Tw target
                        sum_wgt_Tw = 0.0
                            
                        do i = hi_Gate, low_Gate        !Set flow at trial opening, iOut, to (blended + min) and estimate Tw.  Use min flow etemps for all other openings       
                            do kk = 1, nOutlets(i)                          
                                jst = jstOut(i, kk)  

                                if(i == iGate .and. kk == iOut) then
                                    qstr(jst,jb) = q_blended + min_flow(jst) 
                                    if (qstr(jst,jb) > 0.0) then
                                        call downstream_withdrawal_estimate(jst,etemp,estr(jst,jb))       !aeb Uses qstr(flow) not qstrfrac
                                    else
                                        etemp = 0.0         
                                    end if   
                                    sum_wgt_Tw = sum_wgt_Tw + etemp * qstr(jst,jb)
                                else
                                    qstr(jst,jb) = min_flow(jst)    !Set flow to min for debugging clarity
                                    sum_wgt_Tw = sum_wgt_Tw + etemp_minQ(jst) * min_flow(jst)                                 
                                end if
                            end do
                        end do
                            
                        jst = jstOut(iGate, iOut)
                        Tw_out(jst) = sum_wgt_Tw / q_Gates          ! Weighted temperature of outlet, iOut, in iGate assuming it gets all but min flows
                        if(Tw_out(jst) <= (ttarg)) then             ! Don't use LWR_GATE to test here.  We've already chosen gates.  
                            lowBlendGate = iGate                    ! lowBlendGate is the gate containing the low-blending structure. It may not be low_Gate.
                            low_Str = iOut
                            exit
                        end if
                    end do
                    if(low_Str > 0) exit                            ! Opening found; etemp_minQ and min_flow are set for each opening
                end do                        

!*********************************** Set structures (jstr) to blend **********************************                                     
                if(low_Str == 0) then                           !No opening found. Target Tw is less than any opening. Use lowest available opening and no blending.                         
                    jStr_low = jstOut(low_Gate, nOutlets(low_Gate))
                    jStr_hi = 0                                 !No blending
                else 
                    jStr_low = jstOut(lowBlendGate, low_Str)
                    if(low_Str == 1) then                       !Top of gate
                        if(lowBlendGate == low_Gate .and. hi_Gate /= low_Gate) then
                            jStr_hi = jstOut(hi_Gate, nOutlets(hi_Gate))
                        else
                            jStr_hi = 0                         !No blending
                        end if
                    else                                        !Within gate. There is an opening at this gate (level) for blending.
                        jStr_hi = jstOut(lowBlendGate, low_Str - 1)
                    end if
                end if                    
                
                !ZZ 10/24
                if(jStr_low == jst_TCDdown(j)) then       !TCDdown outlet (lowest outlet in side gates) selected for blended release.  Make adjustments.
                    jStr_low = jstOut(TCD_DOWN_GATE, nOutlets(TCD_DOWN_GATE)-1)     !Select outlet above
                    jStr_hi = 0                                                     !No blending.
                end if                
                
                ! Now we have structures to blend
                ! End find gates and structures that bracket target
                                                   
                !aeb	Re-initialize all gate flows to zero for debuging clarity.
                do iGate = minGateNo(j), maxGateNo(j)
                    do iOut = 1, nOutlets(iGate)
                        jst = jstOut(iGate, iOut)
                        qstr(jst,jb) = 0.0                   
                    end do
                end do

                !aeb Set min flows for all openings in blending gates and remove them from blended flow.  Adjust Tw target.
                q_blended = q_Gates
                sum_wgt_Tw = 0.0
                
                ! Remove min flows at all outlets from blended flow. Get sum_wgt_Tw for all non-blended outlets in blending gates
                do iGate = hi_Gate, low_Gate
                    do kk = 1, nOutlets(iGate)
                        jst = jstOut(iGate, kk)
                        if(jst /= jStr_low .and. jst /= jStr_hi) then
                            q_blended = q_blended - min_flow(jst)
                            qstr(jst,jb) = min_flow(jst)
                            sum_wgt_Tw = sum_wgt_Tw + etemp_minQ(jst)*min_flow(jst)  !Check.
                        end if
                    end do
                end do
                ! Blended flow is now reduced by sum of min flows at all non-blending openings

                ! Calculate new temperature target for blended releases, factoring in min flows in non-blending openings
                ttarg = (ttarg * q_Gates - sum_wgt_Tw) / q_blended          

!*********************************** Apportion blending flow using a bisection method                
                If (jStr_hi > 0) then                                               !There is blending
                    q_blended = q_blended - min_flow(jStr_hi) -min_flow(jStr_low)   ! Remove min flows
                    max_iter = 6        ! Max iterations for blending search; max_iter=5 narrows search to a 3% interval
                    f_high = 1			                                        
                    f_low = 0
                
                    Do i=1, max_iter                                                ! Search using a bisection-type method
	                    high_frac = (f_high+f_low)/2		
	                    q_hi = q_blended*high_frac + min_flow(jStr_hi)                ! Add min flows back to blended flow
	                    q_low = q_blended*(1-high_frac) + min_flow(jStr_low)
                            
	                    ! Get Tw high	
                        if(q_hi > 0.) then
                            jst = jStr_hi
                            qstr(jst,jb) = q_hi
                            call downstream_withdrawal_estimate(jst,etemp,estr(jst,jb))     !aeb Uses qstr (flow)
                            Tw_hi = etemp                          
                        else
                            Tw_hi = 0.0
                        end if
                            
	                    ! Get Tw low
                        if(q_low > 0.) then
                            jst = jStr_low
                            qstr(jst,jb) = q_low
                            call downstream_withdrawal_estimate(jst,etemp,estr(jst,jb))     !aeb Uses qstr (flow)
                            Tw_low = etemp								  
                        else
                            Tw_low = 0.0
                        end if

                        ! Check for convergence to Tw target
	                      tm1 = (q_hi*Tw_hi + q_low*Tw_low)/(q_hi + q_low)	
                            
                        if(abs(ttarg-tm1) <= tsconv) then
                            exit
	                    else if(tm1 > ttarg) then		
		                    f_high= high_frac	
	                    else	! (Tw_trial <= Tw_target) 	
		                    f_low= high_frac	
	                    endif		
                    end do	
                    
                    ! If forecasting and solution not converged, check to see if there is a solution that is hotter but closer to Tw target
                    ! Two conditions may be better:
                    !   1) Blending between top of low gate and bottom of hi gate.  Try hi gate only with all blended flow through lowest structure of hi gate.
                    !   2) Blending between 1st and 2nd outlet of low gate. Try top outlet of low gate and min flows in hi gate and other low-gate outlets.
                    ! Don't do this here. Add back in at a later date if it works with new gate logic (e.g., gates are pre-selected at 8 am).
                    ! I don't think this is useful with current logic because it might work for daily setting of gates, but the rest of the day will be off a little anyway.
                    ! Also, new gate logic specifies gates (based on consecutive attempts and off-limits), so gates can't be eliminated at this point.

                    jst = jStr_hi
                    qstr(jst,jb) = q_hi
                    
                    jst = jStr_low
                    qstr(jst,jb) = q_low
                else
                    ! No blending.
                    jst = jStr_low
                    qstr(jst,jb) = q_blended                   
                end if

                ! Set final flow fractions.  All flow now is either unblended (untouched) or blended between 2 outlets (with min flows in associated outlets)
                do jj=1,nouts(j)
                    jst = jstsplt(j,jj)
                    qstrfrac(jst,jb) = qstr(jst,jb)/qall
                end do                
            end if  
            
            if (FORECASTING) qstr = orig_qstr                   ! Re-set flows to non-averaged values
            
          end if                                                !aeb End structures and withdrawals
      end if                                                    !aeb end of this blending group
    end do                                                      !aeb end cycle through each blending group, numtsplt

    nxtsplit = nxtsplit + tspltfreq                     
    ! After blending, both blended and non-blended flows for structures are set.  Non-blended flows have not been changed.
  end if
  
  ! Use the flow fractions to set flows in blended groups.    
  if (tspltc=='      ON') then                                            !aeb  Search for current blending period
      do j=1,numtsplt                                                     !aeb  We have daytest already, don't we?
          if (tsyearly(j) == '     OFF') then
              daytest = jday
          else
              daytest = real(jdayg) + jday - int(jday)
          end if
            
          if (daytest >= tstsrt(j) .and. daytest < tstend(j)) then        !aeb This selects only current blending period
              ! Do structures first
              if (tspltcntr(j) == '      ST') then                 
                  q_Gates = 0.0
                  gate_fraction = 0.0
                  qall = 0.0

                  ! Sum current flows for blended outlets and prior blended flow fractions
                  do jj=1,nouts(j)
                      jst = jstsplt(j,jj)          !ZZ 10/24
                      if(tsprior(j,jj) /= -1) then                          
                          q_Gates = q_Gates + qstr(jst,tspltjb(j))                    
                          gate_fraction = gate_fraction + qstrfrac(jst,tspltjb(j))      
                      endif
                  end do
                    
                  ! Set flows for blended outlets based on prior flow fractions and sum all flows. Honor the maximum flow
                  do jj=1,nouts(j)            
                      jst = jstsplt(j,jj)
                      if(tsprior(j,jj) /= -1) then                ! Set new flow only for blended outlets
                          jst = jstsplt(j,jj)
                          qstr(jst,tspltjb(j)) = 0.0
                          if(gate_fraction > 0.0) qstr(jst,tspltjb(j)) = qstrfrac(jst,tspltjb(j))*q_Gates/gate_fraction   ! In case gate_fraction == 0?
                      end if
                      qall = qall + qstr(jstsplt(j,jj),tspltjb(j))                        
                  end do   
                    
                  ! Set new flow fractions
                  do jj=1,nouts(j)
                      jst = jstsplt(j,jj)
                      qstrfrac(jst,tspltjb(j)) = 0.0
                      if (qall > 0.0) qstrfrac(jst,tspltjb(j)) = qstr(jst,tspltjb(j))/qall   ! In case qall == 0?
                      if (tsmaxflow(j,jj) > 0.0 .and. qstr(jst,tspltjb(j)) > tsmaxflow(j,jj)) qstr(jst,tspltjb(j)) = tsmaxflow(j,jj)
                  end do
              end if
          end if
      end do
  end if      ! End setting flows in blended groups        

  ! Output some results.
  if (jday.ge.nxtstr) then
        nxtstr = nxtstr+tfrqtmp
        ifile=1949
        do jb=1,nbr
            if (nstr(jb) > 0) then
                ifile=ifile+1
                write (ifile,'(f10.4,",",<nstr(jb)>(f10.2,","),<nstr(jb)>(f10.2,","),<nstr(jb)>(f10.2,","))')  &
                              jday,(tavg(i,jb),i=1,nstr(jb)),(qstr(i,jb),i=1,nstr(jb)),(estr(i,jb),i=1,nstr(jb))    
            end if
        end do

        ! computing reservoir volume and volume below 'tempcrit'
        volmc=0.0
        volm=0.0
        DO JW=1,NWB
            KT = KTWB(JW)
            DO JB=BS(JW),BE(JW)
                DO I=cus(jb),ds(jb)
                    volm(jw) = volm(jw) +BH2(KT,I)*DLX(I)
                    DO K=kt+1,kb(i)
                        volm(jw) = volm(jw)+BH(K,I)*DLX(I)
                    END DO
                    do kk=1,tempn
                        if(t2(kt,i).le.tempcrit(jw,kk))volmc(jw,kk) = volmc(jw,kk)+BH2(KT,I)*DLX(I)
                        DO K=kt+1,kb(i)
                            if(t2(k,i).le.tempcrit(jw,kk))volmc(jw,kk) = volmc(jw,kk)+BH(K,I)*DLX(I)
                        END DO
                    end do
                end do
            end do
            ifile=ifile+1
            write(ifile,'(f8.2,100(g12.4,g12.4))') jday,volm(jw),(volmc(jw,kk), kk=1,tempn)
        end do
  end if      
  return

  ENTRY DEALLOCATE_SELECTIVETCD
    DEAllocate (tcnelev,tcjb,tcjs, tcelev,tctemp,tctend,tctsrt,ncountc,tciseg,tcklay,tcelevcon,elcontspl)
    DEAllocate (tspltjb,tspltt,nouts,jstsplt,kstrsplt,tcyearly, tcntr,tspltcntr)
    DEallocate (volm,ncountcw,qwdfrac,qstrfrac)
    DEallocate (tempcrit,volmc,DYNSEL,SELD,NXSEL,TEMP2,TSYEARLY,TSTEND,TSTSRT)
    deallocate (tsdepth, tstype, tsminfrac, tsprior, tsminhead, tsmaxhead, tsmaxflow, no_flow)
    deallocate (tsdynsel, tsseld, nxtssel, tstemp2, ewdsav, estrsav, share_flow, wd_active, str_active)
    deallocate (nout0, nout1, nout2, minfrac1, minfrac2, maxfrac1, maxfrac2, splt2t, splt2e)
    DEALLOCATE (DYNSF, DYNTF)    !ZZ 8/24
  RETURN
    
End Subroutine SelectiveTCD  