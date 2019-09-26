	subroutine NELMIN ( fn, n, start, xmin, ynewlo, reqmin,
     .	step, konvge, kcount, icount, numres, ifault )
!! NELMIN minimizes a function using the Nelder-Mead algorithm.
!
!  Discussion:
!
!    This routine seeks the minimum value of a user-specified function.
!
!    Simplex function minimisation procedure due to Nelder and Mead (1965),
!    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
!    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
!    25, 97) and Hill(1978, 27, 380-2)
!
!    The function to be minimized must be defined by a function of
!    the form
!
!      function fn ( x, f )
!      real ( kind = 8 ) fn
!      real ( kind = 8 ) x(*)
!
!    and the name of this subroutine must be declared EXTERNAL in the
!    calling routine and passed as the argument FN.
!
!    This routine does not include a termination test using the
!    fitting of a quadratic surface.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by R ONeill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Nelder, Roger Mead,
!    A simplex method for function minimization,
!    Computer Journal,
!    Volume 7, 1965, pages 308-313.
!
!    R ONeill,
!    Algorithm AS 47:
!    Function Minimization Using a Simplex Procedure,
!    Applied Statistics,
!    Volume 20, Number 3, 1971, pages 338-345.
!
!  Parameters:
!
!    Input, external FN, the name of the function which evaluates
!    the function to be minimized.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!    0 < N is required.
!
!    Input/output, real ( kind = 8 ) START(N).  On input, a starting point
!    for the iteration.  On output, this data may have been overwritten.
!
!    Output, real ( kind = 8 ) XMIN(N), the coordinates of the point which
!    is estimated to minimize the function.
!
!    Output, real ( kind = 8 ) YNEWLO, the minimum value of the function.
!
!    Input, real ( kind = 8 ) REQMIN, the terminating limit for the variance
!    of the function values.  0 < REQMIN is required.
!
!    Input, real ( kind = 8 ) STEP(N), determines the size and shape of the
!    initial simplex.  The relative magnitudes of its elements should reflect
!    the units of the variables.
!
!    Input, integer ( kind = 4 ) KONVGE, the convergence check is carried out
!    every KONVGE iterations. 0 < KONVGE is required.
!
!    Input, integer ( kind = 4 ) KCOUNT, the maximum number of function
!    evaluations.
!
!    Output, integer ( kind = 4 ) ICOUNT, the number of function evaluations
!    used.
!
!    Output, integer ( kind = 4 ) NUMRES, the number of restarts.
!
!    Output, integer ( kind = 4 ) IFAULT, error indicator.
!    0, no errors detected.
!    1, REQMIN, N, or KONVGE has an illegal value.
!    2, iteration terminated because KCOUNT was exceeded without convergence.
!
	implicit none

	integer ( kind = 4 ) n

	real ( kind = 8 ), parameter :: ccoeff = 0.5D+00
	real ( kind = 8 ) del
	real ( kind = 8 ), parameter :: ecoeff = 2.0D+00
	real ( kind = 8 ), parameter :: eps = 0.001D+00
	real ( kind = 8 ), external :: fn
	integer ( kind = 4 ) i
	integer ( kind = 4 ) icount
	integer ( kind = 4 ) ifault
	integer ( kind = 4 ) ihi
	integer ( kind = 4 ) ilo
	integer ( kind = 4 ) j
	integer ( kind = 4 ) jcount
  	integer ( kind = 4 ) kcount
  	integer ( kind = 4 ) konvge
  	integer ( kind = 4 ) l
  	integer ( kind = 4 ) numres
  	real ( kind = 8 ) p(n,n+1)
  	real ( kind = 8 ) p2star(n)
  	real ( kind = 8 ) pbar(n)
  	real ( kind = 8 ) pstar(n)
  	real ( kind = 8 ), parameter :: rcoeff = 1.0D+00
  	real ( kind = 8 ) reqmin
  	real ( kind = 8 ) rq
  	real ( kind = 8 ) start(n)
  	real ( kind = 8 ) step(n)
  	real ( kind = 8 ) x
  	real ( kind = 8 ) xmin(n)
  	real ( kind = 8 ) y(n+1)
  	real ( kind = 8 ) y2star
  	real ( kind = 8 ) ylo
  	real ( kind = 8 ) ynewlo
  	real ( kind = 8 ) ystar
  	real ( kind = 8 ) z
  	
!
!  Check the input parameters.
!
  	if ( reqmin <= 0.0D+00 ) then
    		ifault = 1
    		return
  	end if

  	if ( n < 1 ) then
    		ifault = 1
    		return
  	end if

  	if ( konvge < 1 ) then
    		ifault = 1
    		return
  	end if
!
!  Initialization.
!
  	icount = 0
  	numres = 0
  	jcount = konvge
  	del = 1.0D+00
  	rq = reqmin * real ( n, kind = 8 )
!
!  Initial or restarted loop.
!
  	do

    		p(1:n,n+1) = start(1:n)
    		y(n+1) = fn ( start )
    		icount = icount + 1
!
!  Define the initial simplex.
!
    		do j = 1, n
      		x = start(j)
      		start(j) = start(j) + step(j) * del
      		p(1:n,j) = start(1:n)
      		y(j) = fn ( start )
      		icount = icount + 1
      		start(j) = x
    		end do
!
!  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
!  the vertex of the simplex to be replaced.
!
    		ilo = minloc ( y(1:n+1), 1 )
    		ylo = y(ilo)
!
!  Inner loop.
!
    		do while ( icount < kcount )
!
!  YNEWLO is, of course, the HIGHEST value???
!
      		ihi = maxloc ( y(1:n+1), 1 )
      		ynewlo = y(ihi)
!
!  Calculate PBAR, the centroid of the simplex vertices
!  excepting the vertex with Y value YNEWLO.
!
      		do i = 1, n
        			pbar(i) = ( sum ( p(i,1:n+1) ) - p(i,ihi) ) / 
     .        			real ( n, kind = 8 )
      		end do
!
!  Reflection through the centroid.
!
      		pstar(1:n) = pbar(1:n) + rcoeff * ( pbar(1:n) - p(1:n,ihi) )
      		ystar = fn ( pstar )
      		icount = icount + 1
!
!  Successful reflection, so extension.
!
      		if ( ystar < ylo ) then

        			p2star(1:n) = pbar(1:n) + 
     .        				ecoeff * ( pstar(1:n) - pbar(1:n) )
        			y2star = fn ( p2star )
        			icount = icount + 1
!
!  Retain extension or contraction.
!
        			if ( ystar < y2star ) then
          				p(1:n,ihi) = pstar(1:n)
          				y(ihi) = ystar
        			else
          				p(1:n,ihi) = p2star(1:n)
          				y(ihi) = y2star
        			end if
!
!  No extension.
!
      		else

        			l = 0
        			do i = 1, n + 1
          				if ( ystar < y(i) ) then
            				l = l + 1
          				end if
        			end do

        			if ( 1 < l ) then

          				p(1:n,ihi) = pstar(1:n)
          				y(ihi) = ystar
!
!  Contraction on the Y(IHI) side of the centroid.
!
        			else if ( l == 0 ) then

          				p2star(1:n) = pbar(1:n) + 
     .          				ccoeff * ( p(1:n,ihi) - pbar(1:n) )
          				y2star = fn ( p2star )
          				icount = icount + 1
!
!  Contract the whole simplex.
!
          				if ( y(ihi) < y2star ) then

            				do j = 1, n + 1
              					p(1:n,j) = ( p(1:n,j) +
     .         						p(1:n,ilo) ) * 0.5D+00
              					xmin(1:n) = p(1:n,j)
              					y(j) = fn ( xmin )
              					icount = icount + 1
            				end do

            				ilo = minloc ( y(1:n+1), 1 )
            				ylo = y(ilo)

            				cycle
!
!  Retain contraction.
!
          				else
            				p(1:n,ihi) = p2star(1:n)
            				y(ihi) = y2star
          				end if
!
!  Contraction on the reflection side of the centroid.
!
        			else if ( l == 1 ) then

          				p2star(1:n) = pbar(1:n) + 
     .					ccoeff * ( pstar(1:n) - pbar(1:n) )
          				y2star = fn ( p2star )
          				icount = icount + 1
!
!  Retain reflection?
!
          				if ( y2star <= ystar ) then
            				p(1:n,ihi) = p2star(1:n)
            				y(ihi) = y2star
          				else
            				p(1:n,ihi) = pstar(1:n)
            				y(ihi) = ystar
          				end if

        			end if

      		end if
!
!  Check if YLO improved.
!
      		if ( y(ihi) < ylo ) then
        			ylo = y(ihi)
        			ilo = ihi
      		end if

      		jcount = jcount - 1

      		if ( 0 < jcount ) then
        			cycle
      		end if
!
!  Check to see if minimum reached.
!
      		if ( icount <= kcount ) then

        			jcount = konvge

        			x = sum ( y(1:n+1) ) / real ( n + 1, kind = 8 )
        			z = sum ( ( y(1:n+1) - x )**2 )

        			if ( z <= rq ) then
          				exit
        			end if

      		end if

    		end do
!
!  Factorial tests to check that YNEWLO is a local minimum.
!
    		xmin(1:n) = p(1:n,ilo)
    		ynewlo = y(ilo)

    		if ( kcount < icount ) then
      		ifault = 2
      		exit
    		end if

    		ifault = 0

    		do i = 1, n
      		del = step(i) * eps
      		xmin(i) = xmin(i) + del
      		z = fn ( xmin )
      		icount = icount + 1
      		if ( z < ynewlo ) then
        			ifault = 2
        			exit
      		end if
      		xmin(i) = xmin(i) - del - del
      		z = fn ( xmin )
      		icount = icount + 1
      		if ( z < ynewlo ) then
        			ifault = 2
        			exit
      		end if
      		xmin(i) = xmin(i) + del
    		end do

    		if ( ifault == 0 ) then
      		exit
    		end if
!
!  Restart the procedure.
!
    		start(1:n) = xmin(1:n)
    		del = eps
    		numres = numres + 1

  	end do

  	return
	end
	
C	subroutine timestamp_nelmin ( )
C!*****************************************************************************80
C!
C!! TIMESTAMP prints the current YMDHMS date as a time stamp.
C!
C!  Example:
C!
C!    31 May 2001   9:45:54.872 AM
C!
C!  Licensing:
C!
C!    This code is distributed under the GNU LGPL license.
C!
C!  Modified:
C!
C!    18 May 2013
C!
C!  Author:
C!
C!    John Burkardt
C!
C!  Parameters:
C!
C!    None
C!
C  	implicit none
C
C  	character ( len = 8 ) ampm
C  	integer ( kind = 4 ) d
C  	integer ( kind = 4 ) h
C  	integer ( kind = 4 ) m
C  	integer ( kind = 4 ) mm
C  	character ( len = 9 ), parameter, dimension(12) :: month = (/ 
C     .	'January  ', 'February ', 'March    ', 'April    ',
C     .	'May      ', 'June     ', 'July     ', 'August   ',
C     .	'September', 'October  ', 'November ', 'December ' /)
C  	integer ( kind = 4 ) n
C  	integer ( kind = 4 ) s
C  	integer ( kind = 4 ) values(8)
C  	integer ( kind = 4 ) y
C
C  	call date_and_time ( values = values )
C
C  	y = values(1)
C  	m = values(2)
C  	d = values(3)
C  	h = values(5)
C  	n = values(6)
C  	s = values(7)
C  	mm = values(8)
C
C  	if ( h < 12 ) then
C    		ampm = 'AM'
C  	else if ( h == 12 ) then
C    		if ( n == 0 .and. s == 0 ) then
C      		ampm = 'Noon'
C    		else
C      		ampm = 'PM'
C    		end if
C  	else
C    		h = h - 12
C    		if ( h < 12 ) then
C      		ampm = 'PM'
C    		else if ( h == 12 ) then
C      		if ( n == 0 .and. s == 0 ) then
C        			ampm = 'Midnight'
C      		else
C        			ampm = 'AM'
C      		end if
C    		end if
C  	end if
C
C  	write(*, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )
C     .d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )
C
C  	return
C	end	



C*======================================================OPTIMIZATION CODE
C      REAL*8 FUNCTION PRAXIS(T0,MACHEP,H0,N,PRIN,X,F,FMIN)
CC                             LAST MODIFIED 3/1/73
C      REAL*8 T0,MACHEP,H0,X(N),F,FMIN
C      EXTERNAL F
C      INTEGER PRIN
CC
CC     PRAXIS RETURNS THE MINIMUM OF THE FUNCTION F(X,N) OF N VARIABLES
CC     USING THE PRINCIPAL AXIS METHOD.  THE GRADIENT OF THE FUNCTION IS
CC     NOT REQUIRED.
CC
CC     FOR A DESCRIPTION OF THE ALGORITHM, SEE CHAPTER SEVEN OF
CC     "ALGORITHMS FOR FINDING ZEROS AND EXTREMA OF FUNCTIONS WITHOUT
CC     CALCULATING DERIVATIVES" BY RICHARD P BRENT.
CC
CC     THE PARAMETERS ARE:
CC     T0       IS A TOLERANCE.  PRAXIS ATTEMPTS TO RETURN PRAXIS=F(X)
CC              SUCH THAT IF X0 IS THE TRUE LOCAL MINIMUM NEAR X, THEN
CC              NORM(X-X0) < T0 + SQUAREROOT(MACHEP)*NORM(X).
CC     MACHEP   IS THE MACHINE PRECISION, THE SMALLEST NUMBER SUCH THAT
CC              1 + MACHEP > 1.  MACHEP SHOULD BE 16.**-13 (ABOUT
CC              2.22D-16) FOR REAL*8 ARITHMETIC ON THE IBM 360.
CC     H0       IS THE MAXIMUM STEP SIZE.  H0 SHOULD BE SET TO ABOUT THE
CC              MAXIMUM DISTANCE FROM THE INITIAL GUESS TO THE MINIMUM.
CC              (IF H0 IS SET TOO LARGE OR TOO SMALL, THE INITIAL RATE OF
CC              CONVERGENCE MAY BE SLOW.)
CC     N        (AT LEAST TWO) IS THE NUMBER OF VARIABLES UPON WHICH
CC              THE FUNCTION DEPENDS.
CC     PRIN     CONTROLS THE PRINTING OF INTERMEDIATE RESULTS.
CC              IF PRIN=0, NOTHING IS PRINTED.
CC              IF PRIN=1, F IS PRINTED AFTER EVERY N+1 OR N+2 LINEAR
CC              MINIMIZATIONS.  FINAL X IS PRINTED, BUT INTERMEDIATE X IS
CC              PRINTED ONLY IF N IS AT MOST 4.
CC              IF PRIN=2, THE SCALE FACTORS AND THE PRINCIPAL VALUES OF
CC              THE APPROXIMATING QUADRATIC FORM ARE ALSO PRINTED.
CC              IF PRIN=3, X IS ALSO PRINTED AFTER EVERY FEW LINEAR
CC              MINIMIZATIONS.
CC              IF PRIN=4, THE PRINCIPAL VECTORS OF THE APPROXIMATING
CC              QUADRATIC FORM ARE ALSO PRINTED.
CC     X        IS AN ARRAY CONTAINING ON ENTRY A GUESS OF THE POINT OF
CC              MINIMUM, ON RETURN THE ESTIMATED POINT OF MINIMUM.
CC     F(X,N)   IS THE FUNCTION TO BE MINIMIZED.  F SHOULD BE A REAL*8
CC              FUNCTION DECLARED EXTERNAL IN THE CALLING PROGRAM.
CC     FMIN     IS AN ESTIMATE OF THE MINIMUM, USED ONLY IN PRINTING
CC              INTERMEDIATE RESULTS.
CC     THE APPROXIMATING QUADRATIC FORM IS
CC              Q(X') = F(X,N) + (1/2) * (X'-X)-TRANSPOSE * A * (X'-X)
CC     WHERE X IS THE BEST ESTIMATE OF THE MINIMUM AND A IS
CC              INVERSE(V-TRANSPOSE) * D * INVERSE(V)
CC     (V(*,*) IS THE MATRIX OF SEARCH DIRECTIONS; D(*) IS THE ARRAY
CC     OF SECOND DIFFERENCES).  IF F HAS CONTINUOUS SECOND DERIVATIVES
CC     NEAR X0, A WILL TEND TO THE HESSIAN OF F AT X0 AS X APPROACHES X0.
CC
CC     IT IS ASSUMED THAT ON FLOATING-POINT UNDERFLOW THE RESULT IS SET
CC     TO ZERO.
CC     THE USER SHOULD OBSERVE THE COMMENT ON HEURISTIC NUMBERS AFTER
CC     THE INITIALIZATION OF MACHINE DEPENDENT NUMBERS.
CC
C      LOGICAL ILLC
C      INTEGER NL,NF,KL,KT,KTM
C      REAL*8 S,SL,DN,DMIN,FX,F1,LDS,LDT,T,H,SF,DF,QF1,QD0,QD1,QA,QB,QC
C      REAL*8 M2,M4,SMALL,VSMALL,LARGE,VLARGE,SCBD,LDFAC,T2,DNI,VALUE
C      REAL*8 RANDOM,DSQRT,DABS
CC
CC.....IF N>20 OR IF N<20 AND YOU NEED MORE SPACE, CHANGE '20' TO THE
CC     LARGEST VALUE OF N IN THE NEXT CARD, IN THE CARD 'IDIM=20', AND
CC     IN THE DIMENSION STATEMENTS IN SUBROUTINES MINFIT,MIN,FLIN,QUAD.
CC
C      REAL*8 D(20),Y(20),Z(20),Q0(20),Q1(20),V(20,20)
C      COMMON /GLOBAL/ FX,LDT,DMIN,NF,NL
C     .       /Q/ V,Q0,Q1,QA,QB,QC,QD0,QD1,QF1
CC
CC.....INITIALIZATION.....
CC     MACHINE DEPENDENT NUMBERS:
CC
C      SMALL=MACHEP*MACHEP
C      VSMALL=SMALL*SMALL
C      LARGE=1.D0/SMALL
C      VLARGE=1.D0/VSMALL
C      M2=DSQRT(MACHEP)
C      M4=DSQRT(M2)
CC
CC     HEURISTIC NUMBERS:
CC     IF THE AXES MAY BE BADLY SCALED (WHICH IS TO BE AVOIDED IF
CC     POSSIBLE), THEN SET SCBD=10.  OTHERWISE SET SCBD=1.
CC     IF THE PROBLEM IS KNOWN TO BE ILL-CONDITIONED, SET ILLC=TRUE.
CC     OTHERWISE SET ILLC=FALSE.
CC     KTM IS THE NUMBER OF ITERATIONS WITHOUT IMPROVEMENT BEFORE THE
CC     ALGORITHM TERMINATES.  KTM=4 IS VERY CAUTIOUS; USUALLY KTM=1
CC     IS SATISFACTORY.
CC
C      SCBD=1.D0
C      ILLC=.FALSE.
C      KTM=1
CC
C      LDFAC=0.01D0
C      IF (ILLC) LDFAC=0.1D0
C      KT=0
C      NL=0
C      NF=1
C      FX=F(X,N)
C      QF1=FX
C      T=SMALL+DABS(T0)
C      T2=T
C      DMIN=SMALL
C      H=H0
C      IF (H.LT.100*T) H=100*T
C      LDT=H
CC.....THE FIRST SET OF SEARCH DIRECTIONS V IS THE IDENTITY MATRIX.....
C      DO 20 I=1,N
C           DO 10 J=1,N
C10              V(I,J)=0.0D0
C20         V(I,I)=1.D0
C      D(1)=0.D0
C      QD0=0.D0
C      DO 30 I=1,N
C           Q0(I)=X(I)
C30         Q1(I)=X(I)
C      IF (PRIN.GT.0) CALL PRINT(N,X,PRIN,FMIN)
CC
CC.....THE MAIN LOOP STARTS HERE.....
C40    SF=D(1)
C      D(1)=0.D0
C      S=0.D0
CC
CC.....MINIMIZE ALONG THE FIRST DIRECTION V(*,1).
CC     FX MUST BE PASSED TO MIN BY VALUE.
C      VALUE=FX
C      CALL MIN(N,1,2,D(1),S,VALUE,.FALSE.,F,X,T,MACHEP,H)
C      IF (S.GT.0.D0) GO TO 50
C           DO 45 I=1,N
C45              V(I,1)=-V(I,1)
C50    IF (SF.GT.0.9D0*D(1).AND.0.9D0*SF.LT.D(1)) GO TO 70
C           DO 60 I=2,N
C60              D(I)=0.D0
CC
CC.....THE INNER LOOP STARTS HERE.....
C70    DO 170 K=2,N
C           DO 75 I=1,N
C75              Y(I)=X(I)
C           SF=FX
C           IF (KT.GT.0) ILLC=.TRUE.
C80         KL=K
C           DF=0.D0
CC
CC.....A RANDOM STEP FOLLOWS (TO AVOID RESOLUTION VALLEYS).
CC     PRAXIS ASSUMES THAT RANDOM RETURNS A RANDOM NUMBER UNIFORMLY
CC     DISTRIBUTED IN (0,1).
CC
C           IF(.NOT.ILLC) GO TO 95
C                DO 90 I=1,N
C                     S=(0.1D0*LDT+T2*(10**KT))*(RANDOM(N)-0.5D0)
C                     Z(I)=S
C                     DO 85 J=1,N
C85                        X(J)=X(J)+S*V(J,I)
C90              CONTINUE
C                FX=F(X,N)
C                NF=NF+1
CC
CC.....MINIMIZE ALONG THE "NON-CONJUGATE" DIRECTIONS V(*,K),...,V(*,N)
CC
C95         DO 105 K2=K,N
C                SL=FX
C                S=0.D0
C                VALUE=FX
C                CALL MIN(N,K2,2,D(K2),S,VALUE,.FALSE.,F,X,T,MACHEP,H)
C                IF (ILLC) GO TO 97
C                     S=SL-FX
C                     GO TO 99
C97              S=D(K2)*((S+Z(K2))**2)
C99              IF (DF.GT.S) GO TO 105
C                     DF=S
C                     KL=K2
C105        CONTINUE
C           IF (ILLC.OR.(DF.GE.DABS((100*MACHEP)*FX))) GO TO 110
CC
CC.....IF THERE WAS NOT MUCH IMPROVEMENT ON THE FIRST TRY, SET
CC     ILLC=TRUE AND START THE INNER LOOP AGAIN.....
CC
C           ILLC=.TRUE.
C           GO TO 80
C110        IF (K.EQ.2.AND.PRIN.GT.1) CALL VCPRNT(1,D,N)
CC
CC.....MINIMIZE ALONG THE "CONJUGATE" DIRECTIONS V(*,1),...,V(*,K-1)
CC
C           KM1=K-1
C           DO 120 K2=1,KM1
C           S=0
C           VALUE=FX
C           CALL MIN(N,K2,2,D(K2),S,VALUE,.FALSE.,F,X,T,MACHEP,H)
C120        CONTINUE
C           F1=FX
C           FX=SF
C           LDS=0
C           DO 130 I=1,N
C                SL=X(I)
C                X(I)=Y(I)
C                SL=SL-Y(I)
C                Y(I)=SL
C130             LDS=LDS+SL*SL
C           LDS=DSQRT(LDS)
C           IF (LDS.LE.SMALL) GO TO 160
CC
CC.....DISCARD DIRECTION V(*,KL).
CC     IF NO RANDOM STEP WAS TAKEN, V(*,KL) IS THE "NON-CONJUGATE"
CC     DIRECTION ALONG WHICH THE GREATEST IMPROVEMENT WAS MADE.....
CC
C           KLMK=KL-K
C           IF (KLMK.LT.1) GO TO 141
C           DO 140 II=1,KLMK
C                I=KL-II
C                DO 135 J=1,N
C135                  V(J,I+1)=V(J,I)
C140             D(I+1)=D(I)
C141        D(K)=0
C           DO 145 I=1,N
C145             V(I,K)=Y(I)/LDS
CC
CC.....MINIMIZE ALONG THE NEW "CONJUGATE" DIRECTION V(*,K), WHICH IS
CC     THE NORMALIZED VECTOR:  (NEW X) - (0LD X).....
CC
C           VALUE=F1
C           CALL MIN(N,K,4,D(K),LDS,VALUE,.TRUE.,F,X,T,MACHEP,H)
C           IF (LDS.GT.0.D0) GO TO 160
C                LDS=-LDS
C                DO 150 I=1,N
C150                  V(I,K)=-V(I,K)
C160        LDT=LDFAC*LDT
C           IF (LDT.LT.LDS) LDT=LDS
C           IF (PRIN.GT.0) CALL PRINT(N,X,PRIN,FMIN)
C           T2=0.D0
C           DO 165 I=1,N
C165             T2=T2+X(I)**2
C           T2=M2*DSQRT(T2)+T
CC
CC.....SEE WHETHER THE LENGTH OF THE STEP TAKEN SINCE STARTING THE
CC     INNER LOOP EXCEEDS HALF THE TOLERANCE.....
CC
C           IF (LDT.GT.(0.5*T2)) KT=-1
C           KT=KT+1
C           IF (KT.GT.KTM) GO TO 400
C170   CONTINUE
CC.....THE INNER LOOP ENDS HERE.
CC
CC     TRY QUADRATIC EXTRAPOLATION IN CASE WE ARE IN A CURVED VALLEY.
CC
C171   CALL QUAD(N,F,X,T,MACHEP,H)
C      DN=0.D0
C      DO 175 I=1,N
C           D(I)=1.D0/DSQRT(D(I))
C           IF (DN.LT.D(I)) DN=D(I)
C175   CONTINUE
C      IF (PRIN.GT.3) CALL MAPRNT(1,V,IDIM,N)
C      DO 180 J=1,N
C           S=D(J)/DN
C           DO 180 I=1,N
C180             V(I,J)=S*V(I,J)
CC
CC.....SCALE THE AXES TO TRY TO REDUCE THE CONDITION NUMBER.....
CC
C      IF (SCBD.LE.1.D0) GO TO 200
C           S=VLARGE
C           DO 185 I=1,N
C                SL=0.D0
C                DO 182 J=1,N
C182                  SL=SL+V(I,J)*V(I,J)
C                Z(I)=DSQRT(SL)
C                IF (Z(I).LT.M4) Z(I)=M4
C                IF (S.GT.Z(I)) S=Z(I)
C185        CONTINUE
C           DO 195 I=1,N
C                SL=S/Z(I)
C                Z(I)=1.D0/SL
C                IF (Z(I).LE.SCBD) GO TO 189
C                     SL=1.D0/SCBD
C                     Z(I)=SCBD
C189             DO 190 J=1,N
C190                  V(I,J)=SL*V(I,J)
C195        CONTINUE
CC
CC.....CALCULATE A NEW SET OF ORTHOGONAL DIRECTIONS BEFORE REPEATING
CC     THE MAIN LOOP.
CC     FIRST TRANSPOSE V FOR MINFIT:
CC
C200   DO 220 I=2,N
C           IM1=I-1
C           DO 210 J=1,IM1
C                S=V(I,J)
C                V(I,J)=V(J,I)
C210             V(J,I)=S
C220   CONTINUE
CC
CC.....CALL MINFIT TO FIND THE SINGULAR VALUE DECOMPOSITION OF V.
CC     THIS GIVES THE PRINCIPAL VALUES AND PRINCIPAL DIRECTIONS OF THE
CC     APPROXIMATING QUADRATIC FORM WITHOUT SQUARING THE CONDITION
CC     NUMBER.....
CC
C      CALL MINFIT(IDIM,N,MACHEP,VSMALL,V,D)
CC
CC.....UNSCALE THE AXES.....
CC
C      IF (SCBD.LE.1.D0) GO TO 250
C           DO 230 I=1,N
C                S=Z(I)
C                DO 225 J=1,N
C225                  V(I,J)=S*V(I,J)
C230        CONTINUE
C           DO 245 I=1,N
C                S=0.D0
C                DO 235 J=1,N
C235                  S=S+V(J,I)**2
C                S=DSQRT(S)
C                D(I)=S*D(I)
C                S=1/S
C                DO 240 J=1,N
C240                  V(J,I)=S*V(J,I)
C245        CONTINUE
CC
C250   DO 270 I=1,N
C           DNI=DN*D(I)
C           IF (DNI.GT.LARGE) GO TO 265
C                IF (DNI.LT.SMALL) GO TO 260
C                     D(I)=1/(DNI*DNI)
C                     GO TO 270
C260             D(I)=VLARGE
C                GO TO 270
C265        D(I)=VSMALL
C270   CONTINUE
CC
CC.....SORT THE EIGENVALUES AND EIGENVECTORS.....
CC
C      CALL SORT(IDIM,N,D,V)
C      DMIN=D(N)
C      IF (DMIN.LT.SMALL) DMIN=SMALL
C      ILLC=.FALSE.
C      IF (M2*D(1).GT.DMIN) ILLC=.TRUE.
C      IF (PRIN.GT.1.AND.SCBD.GT.1.D0) CALL VCPRNT(2,Z,N)
C      IF (PRIN.GT.1) CALL VCPRNT(3,D,N)
C      IF (PRIN.GT.3) CALL MAPRNT(2,V,IDIM,N)
CC.....THE MAIN LOOP ENDS HERE.....
CC
C      GO TO 40
CC
CC.....RETURN.....
CC
C400   IF (PRIN.GT.0) CALL VCPRNT(4,X,N)
C      PRAXIS=FX
C      RETURN
C      END
C      SUBROUTINE MINFIT(M,N,MACHEP,TOL,AB,Q)
C      IMPLICIT REAL*8 (A-H,O-Z)
C      REAL*8 MACHEP
C      DIMENSION AB(M,N),Q(N)
CC...AN IMPROVED VERSION OF MINFIT (SEE GOLUB AND REINSCH, 1969)
CC   RESTRICTED TO M=N,P=0.
CC   THE SINGULAR VALUES OF THE ARRAY AB ARE RETURNED IN Q AND AB IS
CC   OVERWRITTEN WITH THE ORTHOGONAL MATRIX V SUCH THAT U.DIAG(Q) = AB.V,
CC   WHERE U IS ANOTHER ORTHOGONAL MATRIX.
C      DIMENSION E(20)
CC...HOUSEHOLDER'S REDUCTION TO BIDIAGONAL FORM...
C      IF (N.EQ.1) GO TO 200
C      EPS = MACHEP
C      G = 0.D0
C      X = 0.D0
C      DO 11 I=1,N
C         E(I) = G
C         S = 0.D0
C         L = I + 1
C         DO 1 J=I,N
C1           S = S + AB(J,I)**2
C         G = 0.D0
C         IF (S.LT.TOL) GO TO 4
C            F = AB(I,I)
C            G = DSQRT(S)
C            IF (F.GE.0.D0) G = -G
C            H = F*G - S
C            AB(I,I)=F-G
C            IF (L.GT.N) GO TO 4
C            DO 3 J=L,N
C               F = 0.D0
C               DO 2 K=I,N
C2                 F = F + AB(K,I)*AB(K,J)
C               F = F/H
C               DO 3 K=I,N
C3                 AB(K,J) = AB(K,J) + F*AB(K,I)
C4        Q(I) = G
C         S = 0.D0
C         IF (I.EQ.N) GO TO 6
C         DO 5 J=L,N
C5           S = S + AB(I,J)*AB(I,J)
C6        G = 0.D0
C         IF (S.LT.TOL) GO TO 10
C            IF (I.EQ.N) GO TO 16
C            F = AB(I,I+1)
C16          G = DSQRT(S)
C            IF (F.GE.0.D0) G = -G
C            H = F*G - S
C            IF (I.EQ.N) GO TO 10
C            AB(I,I+1) = F - G
C            DO 7 J=L,N
C7              E(J) = AB(I,J)/H
C            DO 9 J=L,N
C               S = 0.D0
C               DO 8 K=L,N
C8                 S = S + AB(J,K)*AB(I,K)
C               DO 9 K=L,N
C9                 AB(J,K) = AB(J,K) + S*E(K)
C10       Y = DABS(Q(I)) + DABS(E(I))
C11       IF (Y.GT.X) X = Y
CC...ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS...
C      AB(N,N) = 1.D0
C      G = E(N)
C      L = N
C      DO 25 II=2,N
C         I = N - II + 1
C         IF (G.EQ.0.D0) GO TO 23
C         H = AB(I,I+1)*G
C         DO 20 J=L,N
C20          AB(J,I) = AB(I,J)/H
C         DO 22 J=L,N
C            S = 0.D0
C            DO 21 K=L,N
C21             S = S + AB(I,K)*AB(K,J)
C            DO 22 K=L,N
C22             AB(K,J) = AB(K,J) + S*AB(K,I)
C23       DO 24 J=L,N
C            AB(I,J) = 0.D0
C24          AB(J,I) = 0.D0
C         AB(I,I) = 1.D0
C         G = E(I)
C25       L = I
CC...DIAGONALIZATION OF THE BIDIAGONAL FORM...
C100   EPS = EPS*X
C      DO 150 KK=1,N
C         K = N - KK + 1
C         KT = 0
C101      KT = KT + 1
C         IF (KT.LE.30) GO TO 102
C            E(K) = 0.D0
C            WRITE (6,1000)
C1000        FORMAT (' QR FAILED')
C102      DO 103 LL2=1,K
C            L2 = K - LL2 + 1
C            L = L2
C            IF (DABS(E(L)).LE.EPS) GO TO 120
C            IF (L.EQ.1) GO TO 103
C            IF (DABS(Q(L-1)).LE.EPS) GO TO 110
C103         CONTINUE
CC...CANCELLATION OF E(L) IF L>1...
C110      C = 0.D0
C         S = 1.D0
C         DO 116 I=L,K
C            F = S*E(I)
C            E(I) = C*E(I)
C            IF (DABS(F).LE.EPS) GO TO 120
C            G = Q(I)
CC...Q(I) = H = DSQRT(G*G + F*F)...
C            IF (DABS(F).LT.DABS(G)) GO TO 113
C            IF (F) 112,111,112
C111         H = 0.D0
C            GO TO 114
C112         H = DABS(F)*DSQRT(1 + (G/F)**2)
C            GO TO 114
C113         H = DABS(G)*DSQRT(1 + (F/G)**2)
C114         Q(I) = H
C            IF (H.NE.0.D0) GO TO 115
C               G = 1.D0
C               H = 1.D0
C115         C = G/H
C116         S = -F/H
CC...TEST FOR CONVERGENCE...
C120      Z = Q(K)
C         IF (L.EQ.K) GO TO 140
CC...SHIFT FROM BOTTOM 2*2 MINOR...
C         X = Q(L)
C         Y = Q(K-1)
C         G = E(K-1)
C         H = E(K)
C         F = ((Y - Z)*(Y + Z) + (G - H)*(G + H))/(2*H*Y)
C         G = DSQRT(F*F + 1.0D0)
C         TEMP = F - G
C         IF (F.GE.0.D0) TEMP = F + G
C         F = ((X - Z)*(X + Z) + H*(Y/TEMP - H))/X
CC...NEXT QR TRANSFORMATION...
C         C = 1.D0
C         S = 1.D0
C         LP1 = L + 1
C         IF (LP1.GT.K) GO TO 133
C         DO 132 I=LP1,K
C            G = E(I)
C            Y = Q(I)
C            H = S*G
C            G = G*C
C            IF (DABS(F).LT.DABS(H)) GO TO 123
C            IF (F) 122,121,122
C121         Z = 0.D0
C            GO TO 124
C122         Z = DABS(F)*DSQRT(1 + (H/F)**2)
C            GO TO 124
C123         Z = DABS(H)*DSQRT(1 + (F/H)**2)
C124         E(I-1) = Z
C            IF (Z.NE.0.D0) GO TO 125
C               F = 1.D0
C               Z = 1.D0
C125         C = F/Z
C            S = H/Z
C            F = X*C + G*S
C            G = -X*S + G*C
C            H = Y*S
C            Y = Y*C
C            DO 126 J=1,N
C               X = AB(J,I-1)
C               Z = AB(J,I)
C               AB(J,I-1) = X*C + Z*S
C126            AB(J,I) = -X*S + Z*C
C            IF (DABS(F).LT.DABS(H)) GO TO 129
C            IF (F) 128,127,128
C127         Z = 0.D0
C            GO TO 130
C128         Z = DABS(F)*DSQRT(1 + (H/F)**2)
C            GO TO 130
C129         Z = DABS(H)*DSQRT(1 + (F/H)**2)
C130         Q(I-1) = Z
C            IF (Z.NE.0.D0) GO TO 131
C               F = 1.D0
C               Z = 1.D0
C131         C = F/Z
C            S = H/Z
C            F = C*G + S*Y
C132         X = -S*G + C*Y
C133      E(L) = 0.D0
C         E(K) = F
C         Q(K) = X
C         GO TO 101
CC...CONVERGENCE:  Q(K) IS MADE NON-NEGATIVE...
C140      IF (Z.GE.0.D0) GO TO 150
C         Q(K) = -Z
C         DO 141 J=1,N
C141         AB(J,K) = -AB(J,K)
C150      CONTINUE
C      RETURN
C200   Q(1) = AB(1,1)
C      AB(1,1) = 1.D0
C      RETURN
C      END
C      SUBROUTINE MIN(N,J,NITS,D2,X1,F1,FK,F,X,T,MACHEP,H)
C      IMPLICIT REAL*8 (A-H,O-Z)
C      EXTERNAL F
C      LOGICAL FK
C      REAL*8 MACHEP,X(N),LDT
C      DIMENSION V(20,20),Q0(20),Q1(20)
C      COMMON /GLOBAL/ FX,LDT,DMIN,NF,NL
C     .       /Q/ V,Q0,Q1,QA,QB,QC,QD0,QD1,QF1
CC...THE SUBROUTINE MIN MINIMIZES F FROM X IN THE DIRECTION V(*,J) UNLESS
CC   J IS LESS THAN 1, WHEN A QUADRATIC SEARCH IS MADE IN THE PLANE
CC   DEFINED BY Q0,Q1,X.
CC   D2 IS EITHER ZERO OR AN APPROXIMATION TO HALF F".
CC   ON ENTRY, X1 IS AN ESTIMATE OF THE DISTANCE FROM X TO THE MINIMUM
CC   ALONG V(*,J) (OR, IF J=0, A CURVE).  ON RETURN, X1 IS THE DISTANCE
CC   FOUND.
CC   IF FK=.TRUE., THEN F1 IS FLIN(X1).  OTHERWISE X1 AND F1 ARE IGNORED
CC   ON ENTRY UNLESS FINAL FX IS GREATER THAN F1.
CC   NITS CONTROLS THE NUMBER OF TIMES AN ATTEMPT WILL BE MADE TO HALVE
CC   THE INTERVAL.
C      LOGICAL DZ
C      REAL*8 M2,M4
C      SMALL = MACHEP**2
C      M2 = DSQRT(MACHEP)
C      M4 = DSQRT(M2)
C      SF1 = F1
C      SX1 = X1
C      K = 0
C      XM = 0.D0
C      FM = FX
C      F0 = FX
C      DZ = D2.LT.MACHEP
CC...FIND THE STEP SIZE...
C      S = 0.D0
C      DO 1 I=1,N
C1        S = S + X(I)**2
C      S = DSQRT(S)
C      TEMP = D2
C      IF (DZ) TEMP = DMIN
C      T2 = M4*DSQRT(DABS(FX)/TEMP + S*LDT) + M2*LDT
C      S = M4*S + T
C      IF (DZ.AND.T2.GT.S) T2 = S
C      T2 = DMAX1(T2,SMALL)
C      T2 = DMIN1(T2,.01D0*H)
C      IF (.NOT.FK.OR.F1.GT.FM) GO TO 2
C      XM = X1
C      FM = F1
C2     IF (FK.AND.DABS(X1).GE.T2) GO TO 3
C      TEMP=1.D0
C      IF (X1.LT.0.D0) TEMP=-1.D0
C      X1=TEMP*T2
C      F1 = FLIN(N,J,X1,F,X,NF)
C3     IF (F1.GT.FM) GO TO 4
C      XM = X1
C      FM = F1
C4     IF (.NOT.DZ) GO TO 6
CC...EVALUATE FLIN AT ANOTHER POINT AND ESTIMATE THE SECOND DERIVATIVE...
C      X2 = -X1
C      IF (F0.GE.F1) X2 = 2.D0*X1
C      F2 = FLIN(N,J,X2,F,X,NF)
C      IF (F2.GT.FM) GO TO 5
C         XM = X2
C         FM = F2
C5     D2 = (X2*(F1 - F0)-X1*(F2 - F0))/((X1*X2)*(X1 - X2))
CC...ESTIMATE THE FIRST DERIVATIVE AT 0...
C6     D1 = (F1 - F0)/X1 - X1*D2
C      DZ = .TRUE.
CC...PREDICT THE MINIMUM...
C      IF (D2.GT.SMALL) GO TO 7
C         X2 = H
C         IF (D1.GE.0.D0) X2 = -X2
C         GO TO 8
C7        X2 = (-.5D0*D1)/D2
C8     IF (DABS(X2).LE.H) GO TO 11
C         IF (X2) 9,9,10
C9        X2 = -H
C         GO TO 11
C10       X2 = H
CC...EVALUATE F AT THE PREDICTED MINIMUM...
C11    F2 = FLIN(N,J,X2,F,X,NF)
C      IF (K.GE.NITS.OR.F2.LE.F0) GO TO 12
CC...NO SUCCESS, SO TRY AGAIN...
C         K = K + 1
C         IF (F0.LT.F1.AND.(X1*X2).GT.0.D0) GO TO 4
C         X2 = 0.5D0*X2
C         GO TO 11
CC...INCREMENT THE ONE-DIMENSIONAL SEARCH COUNTER...
C12    NL = NL + 1
C      IF (F2.LE.FM) GO TO 13
C      X2 = XM
C      GO TO 14
C13    FM = F2
CC...GET A NEW ESTIMATE OF THE SECOND DERIVATIVE...
C14    IF (DABS(X2*(X2 - X1)).LE.SMALL) GO TO 15
C         D2 = (X2*(F1-F0) - X1*(FM-F0))/((X1*X2)*(X1 - X2))
C         GO TO 16
C15       IF (K.GT.0) D2 = 0.D0
C16    IF (D2.LE.SMALL) D2 = SMALL
C      X1 = X2
C      FX = FM
C      IF (SF1.GE.FX) GO TO 17
C         FX = SF1
C         X1 = SX1
CC...UPDATE X FOR LINEAR BUT NOT PARABOLIC SEARCH...
C17    IF (J.EQ.0) RETURN
C      DO 18 I=1,N
C18       X(I) = X(I) + X1*V(I,J)
C      RETURN
C      END
C      REAL*8 FUNCTION FLIN (N,J,L,F,X,NF)
C      IMPLICIT REAL*8 (A-H,O-Z)
C      REAL*8 L,X(N)
C      DIMENSION V(20,20),Q0(20),Q1(20)
CC...FLIN IS THE FUNCTION OF ONE REAL VARIABLE L THAT IS MINIMIZED
CC   BY THE SUBROUTINE MIN...
C      COMMON /Q/ V,Q0,Q1,QA,QB,QC,QD0,QD1,QF1
C      DIMENSION T(20)
C      IF (J .EQ. 0) GO TO 2
CC...THE SEARCH IS LINEAR...
C      DO 1 I=1,N
C1        T(I) = X(I) + L*V(I,J)
C      GO TO 4
CC...THE SEARCH IS ALONG A PARABOLIC SPACE CURVE...
C2     QA = (L*(L - QD1))/(QD0*(QD0 + QD1))
C      QB = ((L + QD0)*(QD1 - L))/(QD0*QD1)
C      QC = (L*(L + QD0))/(QD1*(QD0 + QD1))
C      DO 3 I=1,N
C3        T(I) = (QA*Q0(I) + QB*X(I)) + QC*Q1(I)
CC...THE FUNCTION EVALUATION COUNTER NF IS INCREMENTED...
C4     NF = NF + 1
C      FLIN = F(T,N)
C      RETURN
C      END
C      SUBROUTINE SORT(M,N,D,V)
C      REAL*8 D(N),V(M,N)
CC...SORTS THE ELEMENTS OF D(N) INTO DESCENDING ORDER AND MOVES THE
CC   CORRESPONDING COLUMNS OF V(N,N).
CC   M IS THE ROW DIMENSION OF V AS DECLARED IN THE CALLING PROGRAM.
C      REAL*8 S
C      IF (N.EQ.1) RETURN
C      NM1 = N - 1
C      DO 3 I = 1,NM1
C         K=I
C         S = D(I)
C         IP1 = I + 1
C         DO 1 J = IP1,N
C            IF (D(J) .LE. S) GO TO 1
C            K = J
C            S = D(J)
C1           CONTINUE
C         IF (K .LE. I) GO TO 3
C         D(K) = D(I)
C         D(I) = S
C         DO 2 J = 1,N
C            S = V(J,I)
C            V(J,I) = V(J,K)
C2           V(J,K) = S
C3        CONTINUE
C      RETURN
C      END
C      SUBROUTINE QUAD(N,F,X,T,MACHEP,H)
C      IMPLICIT REAL*8 (A-H,O-Z)
C      EXTERNAL F
CC...QUAD LOOKS FOR THE MINIMUM OF F ALONG A CURVE DEFINED BY Q0,Q1,X...
C      REAL*8 X(N),MACHEP,LDT,L
C      DIMENSION V(20,20),Q0(20),Q1(20)
C      COMMON /GLOBAL/ FX,LDT,DMIN,NF,NL
C     .       /Q/ V,Q0,Q1,QA,QB,QC,QD0,QD1,QF1
C      S = FX
C      FX = QF1
C      QF1 = S
C      QD1 = 0.D0
C      DO 1 I=1,N
C         S = X(I)
C         L = Q1(I)
C         X(I) = L
C         Q1(I) = S
C1        QD1 = QD1 + (S-L)**2
C      QD1 = DSQRT(QD1)
C      L = QD1
C      S = 0.D0
C      IF (QD0. LE. 0.D0 .OR. QD1 .LE. 0.D0 .OR. NL .LT. 3*N*N) GO TO 2
C      VALUE=QF1
C      CALL MIN(N,0,2,S,L,VALUE,.TRUE.,F,X,T,MACHEP,H)
C      QA = (L*(L-QD1))/(QD0*(QD0+QD1))
C      QB = ((L+QD0)*(QD1-L))/(QD0*QD1)
C      QC = (L*(L+QD0))/(QD1*(QD0+QD1))
C      GO TO 3
C2     FX = QF1
C      QA = 0.D0
C      QB = QA
C      QC = 1.D0
C3     QD0 = QD1
C      DO 4 I=1,N
C         S = Q0(I)
C         Q0(I) = X(I)
C4        X(I) = (QA*S + QB*X(I)) + QC*Q1(I)
C      RETURN
C      END
C      SUBROUTINE VCPRNT(OPTION,V,N)
C      REAL*8 V(N)
C      INTEGER OPTION
C      GO TO (1,2,3,4),OPTION
C1     WRITE (6,101) (V(I),I=1,N)
C      RETURN
C2     WRITE (6,102) (V(I),I=1,N)
C      RETURN
C3     WRITE (6,103) (V(I),I=1,N)
C      RETURN
C4     WRITE (6,104) (V(I),I=1,N)
C      RETURN
C101   FORMAT (/' THE SECOND DIFFERENCE ARRAY D(*) IS:'/
C     .        (E32.14,4E25.14))
C102   FORMAT (/' THE SCALE FACTORS ARE:'/(E32.14,4E25.14))
C103   FORMAT (/' THE APPROXIMATING QUADRATIC FORM HAS THE PRINCIPAL VALU
C     .ES:'/(E32.14,4E25.14))
C104   FORMAT (/' X IS:',E26.14/(E32.14))
C      END
C      SUBROUTINE PRINT(N,X,PRIN,FMIN)
C      IMPLICIT REAL*8 (A-H,O-Z)
C      INTEGER PRIN
C      REAL*8 X(N),LN,LDT
C      COMMON /GLOBAL/ FX,LDT,DMIN,NF,NL
C      WRITE (6,101) NL,NF,FX
C      IF (FX.LE.FMIN) GO TO 1
C      LN = DLOG10(FX-FMIN)
C      WRITE (6,102) FMIN,LN
C      GO TO 2
C1     WRITE (6,103) FMIN
C2     IF (N.GT.4.AND.PRIN.LE.2) RETURN
C      WRITE (6,104) (X(I),I=1,N)
C      RETURN
C101   FORMAT (/' AFTER',I6,
C     . ' LINEAR SEARCHES, THE FUNCTION HAS BEEN EVALUATED',I6,
C     . ' TIMES.  THE SMALLEST VALUE FOUND IS F(X) = ',E21.14)
C102   FORMAT (' LOG (F(X)-',E21.14,') = ',E21.14)
C103   FORMAT (' LOG (F(X)-',E21.14,') IS UNDEFINED.')
C104   FORMAT (' X IS:',E26.14/(E32.14))
C      END
C      SUBROUTINE MAPRNT(OPTION,V,M,N)
C      REAL*8 V(M,N)
C      INTEGER OPTION,LOW,UPP
CC...THE SUBROUTINE MAPRNT PRINTS THE COLUMNS OF THE NXN MATRIX V
CC   WITH A HEADING AS SPECIFIED BY OPTION.
CC   M IS THE ROW DIMENSION OF V AS DECLARED IN THE CALLING PROGRAM...
C      LOW = 1
C      UPP = 5
C      GO TO (1,2),OPTION
C1     WRITE (6,101)
C101   FORMAT (/' THE NEW DIRECTIONS ARE:')
C      GO TO 3
C2     WRITE (6,102)
C102   FORMAT (' AND THE PRINCIPAL AXES:')
C3     IF (N.LT.UPP) UPP = N
C      DO 4 I=1,N
C4        WRITE (6,104) (V(I,J),J=LOW,UPP)
C      LOW = LOW + 5
C      IF (N.LT.LOW) RETURN
C      UPP = UPP + 5
C      WRITE (6,103)
C      GO TO 3
C103   FORMAT (' ')
C104   FORMAT (E32.14,4E25.14)
C      END
C      REAL*8 FUNCTION RANDOM(NAUGHT)
C      REAL*8 RAN1,RAN3(127),HALF
C      INTEGER RAN2,Q,R
C      LOGICAL INIT
C      DATA INIT/.FALSE./
C      IF (INIT) GO TO 3
C      R = MOD(NAUGHT,8190) + 1
C      RAN2 = 128
C      DO 2 I=1,127
C         RAN2 = RAN2 - 1
C         RAN1 = -2.D0**55
C         DO 1 J=1,7
C            R = MOD(1756*R,8191)
C            Q = R/32
C1           RAN1 = (RAN1 + Q)*(1.0D0/256)
C2        RAN3(RAN2) = RAN1
C      INIT = .TRUE.
C3     IF (RAN2.EQ.1) RAN2 = 128
C      RAN2 = RAN2 - 1
C      RAN1 = RAN1 + RAN3(RAN2)
C      HALF = .5D0
C      IF (RAN1.GE.0.D0) HALF = -HALF
C      RAN1 = RAN1 + HALF
C      RAN3(RAN2) = RAN1
C      RANDOM = RAN1 + .5D0
C      RETURN
C      END
