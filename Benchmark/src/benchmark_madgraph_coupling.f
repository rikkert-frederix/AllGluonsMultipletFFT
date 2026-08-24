      SUBROUTINE BENCHMARK_SET_STRONG_COUPLING(GS)
      IMPLICIT NONE
      DOUBLE PRECISION GS
      INCLUDE 'coupl.inc'

      G = GS
      GC_10 = -G
      GC_12 = DCMPLX(0D0,1D0)*G**2
      RETURN
      END

C     The uncalled convenience wrappers in MadGraph's generated matrix.f
C     reference these model routines.  The lightweight pure-gluon standalone
C     sets its only two couplings above and therefore needs no parameter card.
      SUBROUTINE SETPARA(PATH)
      IMPLICIT NONE
      CHARACTER*(*) PATH
      RETURN
      END

      SUBROUTINE UPDATE_AS_PARAM()
      IMPLICIT NONE
      INCLUDE 'coupl.inc'
      GC_10 = -G
      GC_12 = DCMPLX(0D0,1D0)*G**2
      RETURN
      END
