      SUBROUTINE GERK(F, NEQN, Y, T, TOUT, RELERR, ABSERR, IFLAG * GERROR, WORK, IWORK)
C
C     FEHLBERG FOURTH(FIFTH) ORDER RUNGE-KUTTA METHOD WITH
C     GLOBAL ERROR ASSESSMENT
C
C     WRITTEN BY H.A.WATTS AND L.F.SHAMPINE
C                   SANDIA LABORATORIES
C
C    GERK IS DESIGNED TO SOLVE SYSTEMS OF DIFFERENTIAL EQUATIONS
C    WHEN IT IS IMPORTANT TO HAVE A READILY AVAILABLE GLOBAL ERROR
C    ESTIMATE.  PARALLEL INTEGRATION IS PERFORMED TO YIELD TWO
C    SOLUTIONS ON DIFFERENT MESH SPACINGS AND GLOBAL EXTRAPOLATION
C    IS APPLIED TO PROVIDE AN ESTIMATE OF THE GLOBAL ERROR IN THE
C    MORE ACCURATE SOLUTION.
C
C    FOR IBM SYSTEM 360 AND 370 AND OTHER MACHINES OF SIMILAR
C    ARITHMETIC CHARACTERISTICS, THIS CODE SHOULD BE CONVERTED TO
C    DOUBLE PRECISION.
C
C*******************************************************************
C ABSTRACT
C*******************************************************************
C
C    SUBROUTINE  GERK  INTEGRATES A SYSTEM OF NEQN FIRST ORDER
C    ORDINARY DIFFERENTIAL EQUATIONS OF THE FORM
C             DY(I)/DT = F(T,Y(1),Y(2),...,Y(NEQN))
C              WHERE THE Y(I) ARE GIVEN AT T .
C    TYPICALLY THE SUBROUTINE IS USED TO INTEGRATE FROM T TO TOUT
C    BUT IT CAN BE USED AS A ONE-STEP INTEGRATOR TO ADVANCE THE
C    SOLUTION A SINGLE STEP IN THE DIRECTION OF TOUT. ON RETURN, AN
C    ESTIMATE OF THE GLOBAL ERROR IN THE SOLUTION AT T IS PROVIDED
C    AND THE PARAMETERS IN THE CALL LIST ARE SET FOR CONTINUING THE
C    INTEGRATION. THE USER HAS ONLY TO CALL GERK AGAIN (AND PERHAPS
C    DEFINE A NEW VALUE FOR TOUT). ACTUALLY, GERK  IS MERELY AN
C    INTERFACING ROUTINE WHICH ALLOCATES VIRTUAL STORAGE IN THE
C    ARRAYS WORK, IWORK AND CALLS SUBROUTINE GERKS FOR THE SOLUTION.
C    GERKS IN TURN CALLS SUBROUTINE FEHL WHICH COMPUTES AN APPROX-
C    IMATE SOLUTION OVER ONE STEP.
C
C    GERK  USES THE RUNGE-KUTTA-FEHLBERG (4,5)  METHOD DESCRIBED
C    IN THE REFERENCE
C    E.FEHLBERG , LOW-ORDER CLASSICAL RUNGE-KUTTA FORMULAS WITH
C                 STEPSIZE CONTROL , NASA TR R-315
C
C
C    THE PARAMETERS REPRESENT-
C      F -- SUBROUTINE F(T,Y,YP) TO EVALUATE DERIVATIVES
C           YP(I)=DY(I)/DT
C      NEQN -- NUMBER OF EQUATIONS TO BE INTEGRATED
C      Y(*) -- SOLUTION VECTOR AT T
C      T -- INDEPENDENT VARIABLE
C      TOUT -- OUTPUT POINT AT WHICH SOLUTION IS DESIRED
C      RELERR,ABSERR -- RELATIVE AND ABSOLUTE ERROR TOLERANCES FOR
C            LOCAL ERROR TEST. AT EACH STEP THE CODE REQUIRES THAT
C                 ABS(LOCAL ERROR) .LE. RELERR*ABS(Y) + ABSERR
C            FOR EACH COMPONENT OF THE LOCAL ERROR AND SOLUTION
C            VECTORS.
C      IFLAG -- INDICATOR FOR STATUS OF INTEGRATION
C      GERROR(*) -- VECTOR WHICH ESTIMATES THE GLOBAL ERROR AT T.
C                   THAT IS, GERROR(I) APPROXIMATES  Y(I)-TRUE
C                   SOLUTION(I).
C      WORK(*) -- ARRAY TO HOLD INFORMATION INTERNAL TO GERK WHICH
C            IS NECESSARY FOR SUBSEQUENT CALLS. MUST BE DIMENSIONED
C            AT LEAST  3+8*NEQN
C      IWORK(*) -- INTEGER ARRAY USED TO HOLD INFORMATION INTERNAL
C            TO GERK WHICH IS NECESSARY FOR SUBSEQUENT CALLS. MUST
C            BE DIMENSIONED AT LEAST  5
C
C
C*******************************************************************
C  FIRST CALL TO GERK
C*******************************************************************
C
C    THE USER MUST PROVIDE STORAGE IN HIS CALLING PROGRAM FOR THE
C    ARRAYS IN THE CALL LIST  -   Y(NEQN), WORK(3+8*NEQN), IWORK(5),
C    DECLARE F IN AN EXTERNAL STATEMENT, SUPPLY SUBROUTINE F(T,Y,YP)
C    AND INITIALIZE THE FOLLOWING PARAMETERS-
C
C      NEQN -- NUMBER OF EQUATIONS TO BE INTEGRATED.  (NEQN .GE. 1)
C      Y(*) -- VECTOR OF INITIAL CONDITIONS
C      T -- STARTING POINT OF INTEGRATION , MUST BE A VARIABLE
C      TOUT -- OUTPUT POINT AT WHICH SOLUTION IS DESIRED.
C            T=TOUT IS ALLOWED ON THE FIRST CALL ONLY,IN WHICH CASE
C            GERK RETURNS WITH IFLAG=2 IF CONTINUATION IS POSSIBLE.
C      RELERR,ABSERR -- RELATIVE AND ABSOLUTE LOCAL ERROR TOLERANCES
C            WHICH MUST BE NON-NEGATIVE BUT MAY BE CONSTANTS. WE CAN
C            USUALLY EXPECT THE GLOBAL ERRORS TO BE SOMEWHAT SMALLER
C            THAN THE REQUESTED LOCAL ERROR TOLERANCES. TO AVOID
C            LIMITING PRECISION DIFFICULTIES THE CODE ALWAYS USES
C            THE LARGER OF RELERR AND AN INTERNAL RELATIVE ERROR
C            PARAMETER WHICH IS MACHINE DEPENDENT.
C      IFLAG -- +1,-1  INDICATOR TO INITIALIZE THE CODE FOR EACH NEW
C            PROBLEM. NORMAL INPUT IS +1. THE USER SHOULD SET IFLAG=
C            -1 ONLY WHEN ONE-STEP INTEGRATOR CONTROL IS ESSENTIAL.
C            IN THIS CASE, GERK ATTEMPTS TO ADVANCE THE SOLUTION A
C            SINGLE STEP IN THE DIRECTION OF TOUT EACH TIME IT IS
C            CALLED. SINCE THIS MODE OF OPERATION RESULTS IN EXTRA
C            COMPUTING OVERHEAD, IT SHOULD BE AVOIDED UNLESS NEEDED.
C
C
C*******************************************************************
C  OUTPUT FROM GERK
C*******************************************************************
C
C      Y(*) -- SOLUTION AT T
C      T -- LAST POINT REACHED IN INTEGRATION.
C      IFLAG = 2 -- INTEGRATION REACHED TOUT.  INDICATES SUCCESSFUL
C                   RETURN AND IS THE NORMAL MODE FOR CONTINUING
C                   INTEGRATION.
C            =-2 -- A SINGLE SUCCESSFUL STEP IN THE DIRECTION OF
C                   TOUT HAS BEEN TAKEN. NORMAL MODE FOR CONTINUING
C                   INTEGRATION ONE STEP AT A TIME.
C            = 3 -- INTEGRATION WAS NOT COMPLETED BECAUSE MORE THAN
C                   9000 DERIVATIVE EVALUATIONS WERE NEEDED. THIS
C                   IS APPROXIMATELY 500 STEPS.
C            = 4 -- INTEGRATION WAS NOT COMPLETED BECAUSE SOLUTION
C                   VANISHED MAKING A PURE RELATIVE ERROR TEST
C                   IMPOSSIBLE. MUST USE NON-ZERO ABSERR TO CONTINUE.
C                   USING THE ONE-STEP INTEGRATION MODE FOR ONE STEP
C                   IS A GOOD WAY TO PROCEED.
C            = 5 -- INTEGRATION WAS NOT COMPLETED BECAUSE REQUESTED
C                   ACCURACY COULD NOT BE ACHIEVED USING SMALLEST
C                   ALLOWABLE STEPSIZE. USER MUST INCREASE THE ERROR
C                   TOLERANCE BEFORE CONTINUED INTEGRATION CAN BE
C                   ATTEMPTED.
C            = 6 -- GERK IS BEING USED INEFFICIENTLY IN SOLVING
C                   THIS PROBLEM. TOO MUCH OUTPUT IS RESTRICTING THE
C                   NATURAL STEPSIZE CHOICE. USE THE ONE-STEP
C                   INTEGRATOR MODE.
C            = 7 -- INVALID INPUT PARAMETERS
C                   THIS INDICATOR OCCURS IF ANY OF THE FOLLOWING IS
C                   SATISFIED - NEQN .LE. 0
C                               T=TOUT  AND  IFLAG .NE. +1 OR -1
C                               RELERR OR ABSERR .LT. 0.
C                               IFLAG .EQ. 0 OR .LT. -2 OR .GT. 7
C      GERROR(*) -- ESTIMATE OF THE GLOBAL ERROR IN THE SOLUTION AT T
C      WORK(*),IWORK(*) -- INFORMATION WHICH IS USUALLY OF NO
C                   INTEREST TO THE USER BUT NECESSARY FOR SUBSEQUENT
C                   CALLS.  WORK(1),...,WORK(NEQN) CONTAIN THE FIRST
C                   DERIVATIVES OF THE SOLUTION VECTOR Y AT T.
C                   WORK(NEQN+1) CONTAINS THE STEPSIZE H TO BE
C                   ATTEMPTED ON THE NEXT STEP.  IWORK(1) CONTAINS
C                   THE DERIVATIVE EVALUATION COUNTER.
C
C
C*******************************************************************
C  SUBSEQUENT CALLS TO GERK
C*******************************************************************
C
C    SUBROUTINE GERK RETURNS WITH ALL INFORMATION NEEDED TO CONTINUE
C    THE INTEGRATION. IF THE INTEGRATION REACHED TOUT, THE USER NEED
C    ONLY DEFINE A NEW TOUT AND CALL GERK AGAIN. IN THE ONE-STEP
C    INTEGRATOR MODE (IFLAG=-2) THE USER MUST KEEP IN MIND THAT EACH
C    STEP TAKEN IS IN THE DIRECTION OF THE CURRENT TOUT.  UPON
C    REACHING TOUT (INDICATED BY CHANGING IFLAG TO 2), THE USER MUST
C    THEN DEFINE A NEW TOUT AND RESET IFLAG TO -2 TO CONTINUE IN THE
C    ONE-STEP INTEGRATOR MODE.
C
C    IF THE INTEGRATION WAS NOT COMPLETED BUT THE USER STILL WANTS
C    TO CONTINUE (IFLAG=3 CASE), HE JUST CALLS GERK AGAIN.  THE
C    FUNCTION COUNTER IS THEN RESET TO 0 AND ANOTHER 9000 FUNCTION
C    EVALUATIONS ARE ALLOWED.
C
C    HOWEVER, IN THE CASE IFLAG=4, THE USER MUST FIRST ALTER THE
C    ERROR CRITERION TO USE A POSITIVE VALUE OF ABSERR BEFORE
C    INTEGRATION CAN PROCEED. IF HE DOES NOT,EXECUTION IS TERMINATED.
C
C    ALSO, IN THE CASE IFLAG=5, IT IS NECESSARY FOR THE USER TO
C    RESET IFLAG TO 2 (OR -2 WHEN THE ONE-STEP INTEGRATION MODE IS
C    BEING USED) AS WELL AS INCREASING EITHER ABSERR,RELERR OR BOTH
C    BEFORE THE INTEGRATION CAN BE CONTINUED. IF THIS IS NOT DONE,
C    EXECUTION WILL BE TERMINATED.  THE OCCURRENCE OF IFLAG=5
C    INDICATES A TROUBLE SPOT (SOLUTION IS CHANGING RAPIDLY,
C    SINGULARITY MAY BE PRESENT) AND IT OFTEN IS INADVISABLE TO
C    CONTINUE.
C
C    IF IFLAG=6 IS ENCOUNTERED, THE USER SHOULD USE THE ONE-STEP
C    INTEGRATION MODE WITH THE STEPSIZE DETERMINED BY THE CODE. IF
C    THE USER INSISTS UPON CONTINUING THE INTEGRATION WITH GERK IN
C    THE INTERVAL MODE, HE MUST RESET IFLAG TO 2 BEFORE CALLING GERK
C    AGAIN.  OTHERWISE,EXECUTION WILL BE TERMINATED.
C
C    IF IFLAG=7 IS OBTAINED, INTEGRATION CAN NOT BE CONTINUED UNLESS
C    THE INVALID INPUT PARAMETERS ARE CORRECTED.
C
C    IT SHOULD BE NOTED THAT THE ARRAYS WORK,IWORK CONTAIN
C    INFORMATION REQUIRED FOR SUBSEQUENT INTEGRATION. ACCORDINGLY,
C    WORK AND IWORK SHOULD NOT BE ALTERED.
C
C*******************************************************************
C
      DIMENSION Y(NEQN), GERROR(NEQN), WORK(1), IWORK(5)
      EXTERNAL F
C COMPUTE INDICES FOR THE SPLITTING OF THE WORK ARRAY
      K1M = NEQN + 1
      K1 = K1M + 1
      K2 = K1 + NEQN
      K3 = K2 + NEQN
      K4 = K3 + NEQN
      K5 = K4 + NEQN
      K6 = K5 + NEQN
      K7 = K6 + NEQN
      K8 = K7 + NEQN
C *******************************************************************
C      THIS INTERFACING ROUTINE MERELY RELIEVES THE USER OF A LONG
C      CALLING LIST VIA THE SPLITTING APART OF TWO WORKING STORAGE
C      ARRAYS. IF THIS IS NOT COMPATIBLE WITH THE USERS COMPILER,
C      HE MUST USE GERKS DIRECTLY.
C *******************************************************************
      CALL GERKS(F, NEQN, Y, T, TOUT, RELERR, ABSERR, IFLAG,
     * GERROR, WORK(1), WORK(K1M), WORK(K1), WORK(K2), WORK(K3),
     * WORK(K4), WORK(K5), WORK(K6), WORK(K7), WORK(K8),
     * WORK(K8+1), IWORK(1), IWORK(2), IWORK(3), IWORK(4), IWORK(5))
      RETURN
      END
      SUBROUTINE GERKS(F, NEQN, Y, T, TOUT, RELERR, ABSERR, IFLAG,      GER   10
     * GERROR, YP, H, F1, F2, F3, F4, F5, YG, YGP, SAVRE, SAVAE,
     * NFE, KOP, INIT, JFLAG, KFLAG)
C      FEHLBERG FOURTH(FIFTH) ORDER RUNGE-KUTTA METHOD WITH
C      GLOBAL ERROR ASSESSMENT
C *******************************************************************
C      GERKS INTEGRATES A SYSTEM OF FIRST ORDER ORDINARY DIFFERENTIAL
C      EQUATIONS AS DESCRIBED IN THE COMMENTS FOR GERK.  THE ARRAYS
C      YP,F1,F2,F3,F4,F5,YG AND YGP (OF DIMENSION AT LEAST NEQN) AND
C      THE VARIABLES H,SAVRE,SAVAE,NFE,KOP,INIT,JFLAG,AND KFLAG ARE
C      USED INTERNALLY BY THE CODE AND APPEAR IN THE CALL LIST TO
C      ELIMINATE LOCAL RETENTION OF VARIABLES BETWEEN CALLS.
C      ACCORDINGLY, THEY SHOULD NOT BE ALTERED. ITEMS OF POSSIBLE
C      INTEREST ARE
C          YP - DERIVATIVE OF SOLUTION VECTOR AT T
C          H  - AN APPROPRIATE STEPSIZE TO BE USED FOR THE NEXT STEP
C          NFE- COUNTER ON THE NUMBER OF DERIVATIVE FUNCTION
C               EVALUATIONS.
C *******************************************************************
      LOGICAL HFAILD, OUTPUT
      DIMENSION Y(NEQN), YP(NEQN), F1(NEQN), F2(NEQN), F3(NEQN),
     * F4(NEQN), F5(NEQN), YG(NEQN), YGP(NEQN), GERROR(NEQN)
      EXTERNAL F
C *******************************************************************
C   THE COMPUTER UNIT ROUNDOFF ERROR U IS THE SMALLEST POSITIVE VALUE
C   REPRESENTABLE IN THE MACHINE SUCH THAT  1.+ U .GT. 1.
C                   VALUES TO BE USED ARE
C               U = 9.5E-7          FOR IBM 360/370
C               U = 1.5E-8          FOR UNIVAC 1108
C               U = 7.5E-9          FOR PDP-10
C               U = 7.1E-15         FOR CDC 6000  SERIES
C               U = 2.2E-16         FOR IBM 360/370  DOUBLE PRECISION
      DATA U /       /
C *******************************************************************
C   REMIN IS A TOLERANCE THRESHOLD WHICH IS ALSO DETERMINED BY THE
C   INTEGRATION METHOD. IN PARTICULAR, A FIFTH ORDER METHOD WILL
C   GENERALLY NOT BE CAPABLE OF DELIVERING ACCURACIES NEAR LIMITING
C   PRECISION ON COMPUTERS WITH LONG WORDLENGTHS.
      DATA REMIN /3.E-11/
C *******************************************************************
C      THE EXPENSE IS CONTROLLED BY RESTRICTING THE NUMBER
C      OF FUNCTION EVALUATIONS TO BE APPROXIMATELY MAXNFE.
C      AS SET,THIS CORRESPONDS TO ABOUT 500 STEPS.
      DATA MAXNFE /9000/
C *******************************************************************
C      CHECK INPUT PARAMETERS
      IF (NEQN.LT.1) GO TO 10
      IF ((RELERR.LT.0.) .OR. (ABSERR.LT.0.)) GO TO 10
      MFLAG = IABS(IFLAG)
      IF ((MFLAG.GE.1) .AND. (MFLAG.LE.7)) GO TO 20
C INVALID INPUT
   10 IFLAG = 7
      RETURN
C IS THIS THE FIRST CALL
   20 IF (MFLAG.EQ.1) GO TO 70
C CHECK CONTINUATION POSSIBILITIES
      IF (T.EQ.TOUT) GO TO 10
      IF (MFLAG.NE.2) GO TO 30
C IFLAG = +2 OR -2
      IF (INIT.EQ.0) GO TO 60
      IF (KFLAG.EQ.3) GO TO 50
      IF ((KFLAG.EQ.4) .AND. (ABSERR.EQ.0.)) GO TO 40
      IF ((KFLAG.EQ.5) .AND. (RELERR.LE.SAVRE) .AND.
     * (ABSERR.LE.SAVAE)) GO TO 40
      GO TO 70
C IFLAG = 3,4,5,6, OR 7
   30 IF (IFLAG.EQ.3) GO TO 50
      IF ((IFLAG.EQ.4) .AND. (ABSERR.GT.0.)) GO TO 60
C INTEGRATION CANNOT BE CONTINUED SINCE USER DID NOT RESPOND TO
C THE INSTRUCTIONS PERTAINING TO IFLAG=4,5,6 OR 7
   40 STOP
C *******************************************************************
C      RESET FUNCTION EVALUATION COUNTER
   50 NFE = 0
      IF (MFLAG.EQ.2) GO TO 70
C RESET FLAG VALUE FROM PREVIOUS CALL
   60 IFLAG = JFLAG
C SAVE INPUT IFLAG AND SET CONTINUATION FLAG VALUE FOR SUBSEQUENT
C INPUT CHECKING
   70 JFLAG = IFLAG
      KFLAG = 0
C SAVE RELERR AND ABSERR FOR CHECKING INPUT ON SUBSEQUENT CALLS
      SAVRE = RELERR
      SAVAE = ABSERR
C RESTRICT RELATIVE ERROR TOLERANCE TO BE AT LEAST AS LARGE AS
C 32U+REMIN TO AVOID LIMITING PRECISION DIFFICULTIES ARISING
C FROM IMPOSSIBLE ACCURACY REQUESTS
      RER = AMAX1(RELERR,32.*U+REMIN)
      U26 = 26.*U
      DT = TOUT - T
      IF (MFLAG.EQ.1) GO TO 80
      IF (INIT.EQ.0) GO TO 90
      GO TO 110
C *******************************************************************
C      INITIALIZATION --
C                        SET INITIALIZATION COMPLETION INDICATOR,INIT
C                        SET INDICATOR FOR TOO MANY OUTPUT POINTS,KOP
C                        EVALUATE INITIAL DERIVATIVES
C                        COPY INITIAL VALUES AND DERIVATIVES FOR THE
C                              PARALLEL SOLUTION
C                        SET COUNTER FOR FUNCTION EVALUATIONS,NFE
C                        ESTIMATE STARTING STEPSIZE
   80 INIT = 0
      KOP = 0
      A = T
      CALL F(A, Y, YP)
      NFE = 1
      IF (T.NE.TOUT) GO TO 90
      IFLAG = 2
      RETURN
   90 INIT = 1
      H = ABS(DT)
      TOLN = 0.
      DO 100 K=1,NEQN
        YG(K) = Y(K)
        YGP(K) = YP(K)
        TOL = RER*ABS(Y(K)) + ABSERR
        IF (TOL.LE.0.) GO TO 100
        TOLN = TOL
        YPK = ABS(YP(K))
        IF (YPK*H**5.GT.TOL) H = (TOL/YPK)**0.2
  100 CONTINUE
      IF (TOLN.LE.0.) H = 0.
      H = AMAX1(H,U26*AMAX1(ABS(T),ABS(DT)))
C *******************************************************************
C      SET STEPSIZE FOR INTEGRATION IN THE DIRECTION FROM T TO TOUT
  110 H = SIGN(H,DT)
C TEST TO SEE IF GERK IS BEING SEVERELY IMPACTED BY TOO MANY
C OUTPUT POINTS
      IF (ABS(H).GT.2.*ABS(DT)) KOP = KOP + 1
      IF (KOP.NE.100) GO TO 120
      KOP = 0
      IFLAG = 6
      RETURN
  120 IF (ABS(DT).GT.U26*ABS(T)) GO TO 140
C IF TOO CLOSE TO OUTPUT POINT,EXTRAPOLATE AND RETURN
      DO 130 K=1,NEQN
        YG(K) = YG(K) + DT*YGP(K)
        Y(K) = Y(K) + DT*YP(K)
  130 CONTINUE
      A = TOUT
      CALL F(A, YG, YGP)
      CALL F(A, Y, YP)
      NFE = NFE + 2
      GO TO 230
C INITIALIZE OUTPUT POINT INDICATOR
  140 OUTPUT = .FALSE.
C TO AVOID PREMATURE UNDERFLOW IN THE ERROR TOLERANCE FUNCTION,
C SCALE THE ERROR TOLERANCES
      SCALE = 2./RER
      AE = SCALE*ABSERR
C *******************************************************************
C *******************************************************************
C      STEP BY STEP INTEGRATION
  150 HFAILD = .FALSE.
C SET SMALLEST ALLOWABLE STEPSIZE
      HMIN = U26*ABS(T)
C ADJUST STEPSIZE IF NECESSARY TO HIT THE OUTPUT POINT.
C LOOK AHEAD TWO STEPS TO AVOID DRASTIC CHANGES IN THE STEPSIZE
C AND THUS LESSEN THE IMPACT OF OUTPUT POINTS ON THE CODE.
      DT = TOUT - T
      IF (ABS(DT).GE.2.*ABS(H)) GO TO 170
      IF (ABS(DT).GT.ABS(H)) GO TO 160
C THE NEXT SUCCESSFUL STEP WILL COMPLETE THE INTEGRATION TO THE
C OUTPUT POINT
      OUTPUT = .TRUE.
      H = DT
      GO TO 170
  160 H = 0.5*DT
C *******************************************************************
C      CORE INTEGRATOR FOR TAKING A SINGLE STEP
C *******************************************************************
C      THE TOLERANCES HAVE BEEN SCALED TO AVOID PREMATURE UNDERFLOW
C      IN COMPUTING THE ERROR TOLERANCE FUNCTION ET.
C      TO AVOID PROBLEMS WITH ZERO CROSSINGS, RELATIVE ERROR IS
C      MEASURED USING THE AVERAGE OF THE MAGNITUDES OF THE SOLUTION
C      AT THE BEGINNING AND END OF A STEP.
C      THE ERROR ESTIMATE FORMULA HAS BEEN GROUPED TO CONTROL LOSS OF
C      SIGNIFICANCE.
C      TO DISTINGUISH THE VARIOUS ARGUMENTS, H IS NOT PERMITTED
C      TO BECOME SMALLER THAN 26 UNITS OF ROUNDOFF IN T.
C      PRACTICAL LIMITS ON THE CHANGE IN THE STEPSIZE ARE ENFORCED TO
C      SMOOTH THE STEPSIZE SELECTION PROCESS AND TO AVOID EXCESSIVE
C      CHATTERING ON PROBLEMS HAVING DISCONTINUITIES.
C      TO PREVENT UNNECESSARY FAILURES, THE CODE USES 9/10 THE
C      STEPSIZE IT ESTIMATES WILL SUCCEED.
C      AFTER A STEP FAILURE, THE STEPSIZE IS NOT ALLOWED TO INCREASE
C      FOR THE NEXT ATTEMPTED STEP.  THIS MAKES THE CODE MORE
C      EFFICIENT ON PROBLEMS HAVING DISCONTINUITIES AND MORE
C      EFFECTIVE IN GENERAL SINCE LOCAL EXTRAPOLATION IS BEING USED
C      AND THE ERROR ESTIMATE MAY BE UNRELIABLE OR UNACCEPTABLE WHEN
C      A STEP FAILS.
C *******************************************************************
C      TEST NUMBER OF DERIVATIVE FUNCTION EVALUATIONS.
C      IF OKAY,TRY TO ADVANCE THE INTEGRATION FROM T TO T+H
  170 IF (NFE.LE.MAXNFE) GO TO 180
C TOO MUCH WORK
      IFLAG = 3
      KFLAG = 3
      RETURN
C ADVANCE AN APPROXIMATE SOLUTION OVER ONE STEP OF LENGTH H
  180 CALL FEHL(F, NEQN, YG, T, H, YGP, F1, F2, F3, F4, F5, F1)
      NFE = NFE + 5
C COMPUTE AND TEST ALLOWABLE TOLERANCES VERSUS LOCAL ERROR
C ESTIMATES AND REMOVE SCALING OF TOLERANCES. NOTE THAT RELATIVE
C ERROR IS MEASURED WITH RESPECT TO THE AVERAGE MAGNITUDES OF THE
C OF THE SOLUTION AT THE BEGINNING AND END OF THE STEP.
      EEOET = 0.
      DO 200 K=1,NEQN
        ET = ABS(YG(K)) + ABS(F1(K)) + AE
        IF (ET.GT.0.) GO TO 190
C INAPPROPRIATE ERROR TOLERANCE
        IFLAG = 4
        KFLAG = 4
        RETURN
  190   EE = ABS((-2090.*YGP(K)+(21970.*F3(K)-15048.*F4(K)))
     *   +(22528.*F2(K)-27360.*F5(K)))
        EEOET = AMAX1(EEOET,EE/ET)
  200 CONTINUE
      ESTTOL = ABS(H)*EEOET*SCALE/752400.
      IF (ESTTOL.LE.1.) GO TO 210
C UNSUCCESSFUL STEP
C                   REDUCE THE STEPSIZE , TRY AGAIN
C                   THE DECREASE IS LIMITED TO A FACTOR OF 1/10
      HFAILD = .TRUE.
      OUTPUT = .FALSE.
      S = 0.1
      IF (ESTTOL.LT.59049.) S = 0.9/ESTTOL**0.2
      H = S*H
      IF (ABS(H).GT.HMIN) GO TO 170
C REQUESTED ERROR UNATTAINABLE AT SMALLEST ALLOWABLE STEPSIZE
      IFLAG = 5
      KFLAG = 5
      RETURN
C SUCCESSFUL STEP
C                    STORE ONE-STEP SOLUTION YG AT T+H
C                    AND EVALUATE DERIVATIVES THERE
  210 TS = T
      T = T + H
      DO 220 K=1,NEQN
        YG(K) = F1(K)
  220 CONTINUE
      A = T
      CALL F(A, YG, YGP)
      NFE = NFE + 1
C NOW ADVANCE THE Y SOLUTION OVER TWO STEPS OF
C LENGTH H/2 AND EVALUATE DERIVATIVES THERE
      HH = 0.5*H
      CALL FEHL(F, NEQN, Y, TS, HH, YP, F1, F2, F3, F4, F5, Y)
      TS = TS + HH
      A = TS
      CALL F(A, Y, YP)
      CALL FEHL(F, NEQN, Y, TS, HH, YP, F1, F2, F3, F4, F5, Y)
      A = T
      CALL F(A, Y, YP)
      NFE = NFE + 12
C CHOOSE NEXT STEPSIZE
C THE INCREASE IS LIMITED TO A FACT-6R OF 5
C IF STEP FAILURE HAS JUST OCCURRED, NEXT
C    STEPSIZE IS NOT ALLOWED TO INCREASE
      S = 5.
      IF (ESTTOL.GT.1.889568E-4) S = 0.9/ESTTOL**0.2
      IF (HFAILD) S = AMIN1(S,1.)
      H = SIGN(AMAX1(S*ABS(H),HMIN),H)
C *******************************************************************
C      END OF CORE INTEGRATOR
C *******************************************************************
C      SHOULD WE TAKE ANOTHER STEP
      IF (OUTPUT) GO TO 230
      IF (IFLAG.GT.0) GO TO 150
C *******************************************************************
C *******************************************************************
C      INTEGRATION SUCCESSFULLY COMPLETED
C      ONE-STEP MODE
      IFLAG = -2
      GO TO 240
C INTERVAL MODE
  230 T = TOUT
      IFLAG = 2
  240 DO 250 K=1,NEQN
        GERROR(K) = (YG(K)-Y(K))/31.
  250 CONTINUE
      RETURN
      END
      SUBROUTINE FEHL(F, NEQN, Y, T, H, YP, F1, F2, F3, F4, F5, S)      FEH   10
C      FEHLBERG FOURTH-FIFTH ORDER RUNGE-KUTTA METHOD
C *******************************************************************
C     FEHL INTEGRATES A SYSTEM OF NEQN FIRST ORDER
C     ORDINARY DIFFERENTIAL EQUATIONS OF THE FORM
C              DY(I)/DT=F(T,Y(1),---,Y(NEQN))
C     WHERE THE INITIAL VALUES Y(I) AND THE INITIAL DERIVATIVES
C     YP(I) ARE SPECIFIED AT THE STARTING POINT T. FEHL ADVANCES
C     THE SOLUTION OVER THE FIXED STEP H AND RETURNS
C     THE FIFTH ORDER (SIXTH ORDER ACCURATE LOCALLY) SOLUTION
C     APPROXIMATION AT T+H IN ARRAY S(I).
C     F1,---,F5 ARE ARRAYS OF DIMENSION NEQN WHICH ARE NEEDED
C     FOR INTERNAL STORAGE.
C     THE FORMULAS HAVE BEEN GROUPED TO CONTROL LOSS OF SIGNIFICANCE.
C     FEHL SHOULD BE CALLED WITH AN H NOT SMALLER THAN 13 UNITS OF
C     ROUNDOFF IN T SO THAT THE VARIOUS INDEPENDENT ARGUMENTS CAN BE
C     DISTINGUISHED.
C *******************************************************************
      DIMENSION Y(NEQN), YP(NEQN), F1(NEQN), F2(NEQN), F3(NEQN),
     * F4(NEQN), F5(NEQN), S(NEQN)
      CH = 0.25*H
      DO 10 K=1,NEQN
        F5(K) = Y(K) + CH*YP(K)
   10 CONTINUE
      CALL F(T+0.25*H, F5, F1)
      CH = 0.09375*H
      DO 20 K=1,NEQN
        F5(K) = Y(K) + CH*(YP(K)+3.*F1(K))
   20 CONTINUE
      CALL F(T+0.375*H, F5, F2)
      CH = H/2197.
      DO 30 K=1,NEQN
        F5(K) = Y(K) + CH*(1932.*YP(K)+(7296.*F2(K)-7200.*F1(K)))
   30 CONTINUE
      CALL F(T+12./13.*H, F5, F3)
      CH = H/4104.
      DO 40 K=1,NEQN
        F5(K) = Y(K) + CH*((8341.*YP(K)-845.*F3(K))+(29440.*F2(K)
     *   -32832.*F1(K)))
   40 CONTINUE
      CALL F(T+H, F5, F4)
      CH = H/20520.
      DO 50 K=1,NEQN
        F1(K) = Y(K) + CH*((-6080.*YP(K)+(9295.*F3(K)-5643.*F4(K)))
     *   +(41040.*F1(K)-28352.*F2(K)))
   50 CONTINUE
      CALL F(T+0.5*H, F1, F5)
C COMPUTE APPROXIMATE SOLUTION AT T+H
      CH = H/7618050.
      DO 60 K=1,NEQN
        S(K) = Y(K) + CH*((902880.*YP(K)+(3855735.*F3(K)-1371249.*
     *   F4(K)))+(3953664.*F2(K)+277020.*F5(K)))
   60 CONTINUE
      RETURN
      END
