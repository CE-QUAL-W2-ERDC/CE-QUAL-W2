Subroutine Plunge_Point
    ! SW 5/2024 Estimate of plunge point in a reservoir
  USE MAIN
  USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  
  IMPLICIT NONE
  REAL :: FRDN
  
  IF(NIT==0)THEN
      OPEN(2020,FILE='Plunge_Point_Estimate.csv',STATUS='unknown')
      write(2020,'(A)')'JDAY,I,JB,KT,RHO(KT I),RHO(KT I+1),T2(KT I),T2(KT I+1),W(KT I),U(KT I),U(KT I+1), Froude#'
  ENDIF
  
  !IF(JDAY>96.0)THEN
  !    PAUSE
  !ENDIF
  
  
  DO JW=1,NWB
      KT=KTWB(JW)
  DO JB=BS(JW),BE(JW)
      IF(BR_INACTIVE(JB))CYCLE
      DO I=CUS(JB),DS(JB)-1
          IF(T2(KT,I) <= 10.0)CYCLE
          
          IF(I==CUS(JB))THEN
                IF(RHO(KT,I) > RHO(KT,I+1) .AND.  U(KT,I) < U(KT-1,I) .AND. U(KB(I),I) > 0.0 .AND. T2(KT,I) < T2(KT,I+1) .and. u(kt,i) > 0.0)THEN                      !W(KT,I) > 0.0 .AND.
                              FRDN=u(kt,i)/(sqrt(((rho(kt,i)-rho(kt,i+1))/rho(kt,i))*g*depthm(kt,i)))
                    WRITE(2020,'(F10.3,",",I4,",",I4,",",I4,",",8(F12.4,","))')JDAY,I,JB,KT,RHO(KT,I),RHO(KT,I+1),T2(KT,I),T2(KT,I+1),W(KT,I),U(KT,I),U(KT,I+1),frdn
                ENDIF
          CYCLE    
          ENDIF
          IF(RHO(KT,I) > RHO(KT,I+1) .AND.  U(KT,I) > 0.0 .AND. U(KT,I+1)< 0.0 .AND. T2(KT,I) < T2(KT,I+1) .and. u(kt,i) > 0.0)THEN                      !W(KT,I) > 0.0 .AND.
                        FRDN=u(kt,i)/(sqrt((rho(kt,i)-rho(kt,i+1))/rho(kt,i))*g*depthm(kt,i))
              WRITE(2020,'(F10.3,",",I4,",",I4,",",I4,",",8(F12.4,","))')JDAY,I,JB,KT,RHO(KT,I),RHO(KT,I+1),T2(KT,I),T2(KT,I+1),W(KT,I),U(KT,I),U(KT,I+1),frdn
          ENDIF
      ENDDO
  ENDDO
  ENDDO
  
  RETURN
  END SUBROUTINE Plunge_Point