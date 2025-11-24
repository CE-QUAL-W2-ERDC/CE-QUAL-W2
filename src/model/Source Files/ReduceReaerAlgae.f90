Subroutine ReduceReaeration
! Reduce Reaeraton as a result of algae accumulation on the surface
! S. Wells August 2024

  USE MAIN; USE SCREENC
  USE GLOBAL;  USE KINETIC; USE AlgaeReduceGasTransfer
 
  
  IMPLICIT NONE
real :: algsum

! COMPUTE TOTAL ALGAE BIOMASS IN SURFACE LAYER


algsum=0.0
DO JA=1,NAL
      IF(ALG_CALC(JA))THEN
      IF(I_ALG(JA)==1)THEN
	  algsum=algsum+alg(kt,i,ja)
      ENDIF
	  endif
reaer(i)=reaer(i)*(1.-algsum/(algsum+khs_alg))
enddo

IF((IMM/IOUTFREQ)*IOUTFREQ == IMM)THEN
WRITE(ALGRED,'(F12.3,",",I5,",",f10.4,",",F10.4)')JDAY,I,algsum,(1.-algsum/(algsum+khs_alg))
ENDIF


return


End Subroutine ReduceReaeration