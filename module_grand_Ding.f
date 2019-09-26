*     Routines needed for leukemia modeling:
*     [PARAMETER_PRAXIS]                  Parameter fitting using PRAXIS
*     [PARAMETER_NELMIN]                  Parameter fitting using NELMIN
*     [PARAMETER_NELMIN_SENSITIVE]        Parameter fitting using NELMIN
*                                         when the result is very
*                                         sensitive to the parameters,
*                                         therefore an extra step is
*                                         taken to fine-tune the initial
*                                         guess; avoid using this if
*                                         possible, since it can cause
*                                         bias
*     [STOCHASTIC]                        Producing stochastic
*                                         trajectories for leukemia with
*                                         chemotherapy, using a hybrid
*                                         SSA - tau leaping algorithm
      MODULE module_grand
      use module_patient
      implicit none
*------------------Global parameters for parameter fitting by any method
*	Bounds on parameters to be fitted
      double precision,parameter	:: min_a_l	= 0.7D0
      double precision,parameter	:: max_a_l	= 1.0D0
      double precision,parameter	:: min_p_l	= 0.3D0
      double precision,parameter	:: max_p_l	= 0.7D0
*	Parameter set being considered in derivative
      double precision,protected	:: PARAMETERS(1:3,1:N_clones)
*	Cancer clonal death rate, assumed to be the same for all clones
      double precision,parameter	:: CANCER_DEATH_RATE	= 0.2D0
*----------------------Global parameters for parameter fitting by PRAXIS
*	Tolerance for PRAXIS
      double precision,parameter	:: T0		= 1.0D-4
*	Machine precision
      double precision,parameter	:: MACHEP	= 2.22D-16
*	Control of the printing of intermediate results. 0 prints nothing,
*	1 prints F after ever N+1 or N+2 linear minimizations and the
*	final result. 2 also prints the scale factors and the principal
*	values of the approximating quadratic form. 3 also prints result
*	after every few linear minimizations. 4 also prints the principal
*	vectors of the approximating quadratic form
      integer,parameter			:: PRIN	= 1
*----------------------Global parameters for parameter fitting by NELMIN
*	Number of parameter sets needed
      integer,parameter			:: N_TOTAL_PARAMETERS	= 100
*	Bounds on parameters to be fitted
      double precision,parameter	:: min_a_l_nelmin	= 0.7D0
      double precision,parameter	:: max_a_l_nelmin	= 1.0D0
      double precision,parameter	:: min_p_l_nelmin	= 0.3D0
      double precision,parameter	:: max_p_l_nelmin	= 0.7D0
*	The terminating limit for the variance of the function values
      double precision,parameter	:: REQMIN	= 1000000.0D0
*	The convergence check which is carried out every KONVGE iterations
*	(has to be >0)
      integer,parameter			:: KONVGE	= 10
*	Maximum number of function evaluations
      integer,parameter			:: KCOUNT	= 1000
*----------------------Global parameters for parameter fitting by NELMIN
*--------------------------------------------------in the sensitive case
*----------------------(additional to the parameters in the normal case)
*	Bounds on number of times correcting the initial guess
      integer,parameter			:: BOUND_CORRECTION	= 500
*----------------------------Global parameters for NEW parameter fitting
*----------------------(additional to the parameters in the normal case)
*     The "center" of the parameter space;
*     the goal is to achieve a parameter set as close to this center as
*     possible
      double precision,parameter    :: center_a_l     = 0.9D0
      double precision,parameter    :: center_p_l     = 0.5D0
*------------------------------------------------Other global parameters
*	A big integer that FORTRAN can handle (up to 1D18)
      integer(kind=8),parameter	:: BIG 	= 1D6
*	Number zero
      double precision,parameter	:: ZERO	= 0.0D0

      double precision              :: MRD,BM_cellularity

      save

      CONTAINS
*=========================COMPUTE DERIVATIVES IN THE DETERMINISTIC MODEL
      SUBROUTINE derivative(t,state,der)
      implicit none
*	Input: state of the system and current time
      double precision	:: state(1:N_total),t
*	Output: derivative according to the deterministic model
      double precision	:: der(1:N_total)
*	Parameters for the algorithm
      double precision	:: sum_total,sum_BM,sum_blood,S,D,var
      double precision	:: a_l,p_l,d_l,k_chemo_c,k_chemo_l
      integer		:: i,j
*--------------------------------------------------Compute the feedbacks
      sum_total	= sum(state)
      sum_blood	= state(2)
      do i=1,N_clones
      	sum_blood	= sum_blood+state(2*i+2)
      enddo
      sum_BM	= sum_total-sum_blood
      S		= 1.0D0/(1.0D0+k_feedback*sum_blood)
      D		= C1_feedback*max(ZERO,sum_BM-C2_feedback*C_1_BAR)
*-------------------------Compute for hematopoietic clone in normal time
      if (state(1).lt.1.0D0) then
      	der(1)	= 0.0D0
      else
      	der(1)	= (2*a_c*S-1.0D0)*p_c*state(1)-D*state(1)
      endif
      der(2)	= 2.0D0*(1-a_c*S)*p_c*state(1)-d_c*state(2)
*-------------------------------Compute for cancer clones in normal time
      do i=1,N_clones
      	a_l	= PARAMETERS(1,i)
      	p_l	= PARAMETERS(2,i)
      	d_l	= PARAMETERS(3,i)
      	j	= 2*i+1
      	if (state(j).lt.1.0D0) then
      		der(j)	= 0.0D0
      	else
      		der(j)	= (2*a_l*S-1.0D0)*p_l*state(j)-D*state(j)
      	endif
      	j	= 2*i+2
            der(j)      = 2.0D0*(1-a_l*S)*p_l*state(j-1)-d_l*state(j)
      enddo
*---------------------------------Add in the terms for chemotherapy time
      do j=1,CHEMO_TOTAL
      	if ((t.ge.CHEMO_TIME(j,1)).and.(t.le.CHEMO_TIME(j,2))) then
      		k_chemo_c	= CHEMO_STRE(j,1)
      		k_chemo_l	= CHEMO_STRE(j,2)
      		if (state(1).ge.1.0D0) then
      			der(1)= der(1)-k_chemo_c*p_c*state(1)
      		endif
      		do i=1,N_clones
      			p_l	= PARAMETERS(2,i)
      			if (state(2*i+1).ge.1.0D0) then
      				der(2*i+1)=der(2*i+1)-
     .                                   k_chemo_l*p_l*state(2*i+1)
      			endif
      		enddo
      		goto 203
      	endif
      enddo
203   END SUBROUTINE
*===========COMPUTE THE FINAL STATE ACCORDING TO THE DETERMINISTIC MODEL
      SUBROUTINE deterministic(final_state,para)
      implicit none
*	Input: parameters for the cancer clones
      double precision			:: para(1:3,1:N_clones)
*	Output: trajectory of the system state and time accordingly, and
*	the position where data ends
      double precision			:: final_state(1:N_total)
*	Parameters for the algorithm
      integer,parameter			:: neqn	= N_total
      double precision,parameter	:: tout	= T_RELAPSE
      double precision			:: RELERR,ABSERR
      integer				:: iflag,iwork(5)
      double precision			:: state(N_total),t
      double precision			:: work(100+21*neqn)
      integer				:: i
*--------------------------------Initialize the parameters for using ode
      RELERR	= 1.0D-3
      ABSERR	= 1.0D-3
      t		= 0.0D0
      iflag		= 1
      PARAMETERS	= para
      state(1)		= diag_BM_sum*diag_clonal(N_clones+1)/100.0D0
      state(2)		= diag_normal
      do i=1,N_clones
      	state(2*i+1)= diag_BM_sum*diag_clonal(i)/100.0D0
      	state(2*i+2)= diag_blast*diag_clonal(i)/100.0D0
      enddo
*---------------------------------
160   call ode(derivative,neqn,state,t,tout,
     .	   RELERR,ABSERR,iflag,work,iwork)
     	if (t.lt.tout) then
     		goto 160
     	endif
     	if (t.ne.tout) then
     		print*,'Something is wrong with SUBROUTINE deterministic'
     		return
     	endif
      final_state	= state
      END SUBROUTINE
*===================================COMPUTE THE MINIMAL RESIDUAL DISEASE
*===================================ACCORDING TO THE DETERMINISTIC MODEL
      SUBROUTINE deterministic_mrd(final_state,para)
      implicit none
*	Input: parameters for the cancer clones
      double precision			:: para(1:3,1:N_clones)
*	Output: trajectory of the system state and time accordingly, and
*	the position where data ends
      double precision			:: final_state(1:N_total)
*	Parameters for the algorithm
      integer,parameter			:: neqn	= N_total
      double precision,parameter	:: tout	= CHEMO_TIME(2,2)
      double precision			:: RELERR,ABSERR
      integer				:: iflag,iwork(5)
      double precision			:: state(N_total),t
      double precision			:: work(100+21*neqn)
      integer				:: i
*--------------------------------Initialize the parameters for using ode
      RELERR	= 1.0D-3
      ABSERR	= 1.0D-3
      t		= 0.0D0
      iflag		= 1
      PARAMETERS	= para
      state(1)		= diag_BM_sum*diag_clonal(N_clones+1)/100.0D0
      state(2)		= diag_normal
      do i=1,N_clones
      	state(2*i+1)= diag_BM_sum*diag_clonal(i)/100.0D0
      	state(2*i+2)= diag_blast*diag_clonal(i)/100.0D0
      enddo
*---------------------------------
160   call ode(derivative,neqn,state,t,tout,
     .	   RELERR,ABSERR,iflag,work,iwork)
     	if (t.lt.tout) then
     		goto 160
     	endif
     	if (t.ne.tout) then
     		print*,'Something is wrong with SUBROUTINE deterministic'
     		return
     	endif
      final_state	= state
      END SUBROUTINE
*========COMPUTE THE ERROR OF A PARAMETER SET, TO BE MINIMIZED BY PRAXIS
      FUNCTION error_praxis(para_set,N)
      implicit none
*	Input: parameters for the cancer clones, and number of variables
*	in the parameter set (always 2*N_clones, for a_l and p_l)
      double precision			:: para_set(1:N)
      integer				:: N
*	Output: infinity-norm error of the percentages at relapse compared
*	to the real percentages
      double precision			:: error_praxis
*	Parameters for the algorithm
      double precision			:: rel_clonal_com(1:N_clones+1)
      double precision			:: difference(1:N_clones+1)
      double precision			:: para(1:3,1:N_clones)
      double precision			:: final_state(N_total)
      double precision			:: sum_BM
      integer				:: i
*-----------------------------------------------------------------Set up
      if (N.ne.(2*N_clones)) then
      	print*,'Something fishy with the parameter set'
      	return
      endif
      if ((minval(para_set(1:N_clones)).lt.min_a_l).or.
     .    (maxval(para_set(1:N_clones)).gt.max_a_l).or.
     .    (minval(para_set(N_clones+1:2*N_clones)).lt.min_p_l).or.
     .    (maxval(para_set(N_clones+1:2*N_clones)).gt.max_p_l)) then
            error_praxis	= 100.0D0
     		goto 197
     	endif
      do i=1,N_clones
*		Renewal rates a_l
      	para(1,i)	= para_set(i)
*		Proliferation rates p_l
      	para(2,i)	= para_set(i+N_clones)
*		Death rates d_l
      	para(3,i)	= CANCER_DEATH_RATE
      enddo
*----------------------------------------------------------Check the MRD
      call deterministic_mrd(final_state,para)
      sum_BM		= final_state(1)
      do i=1,N_clones
            sum_BM	= sum_BM+final_state(2*i+1)
      enddo
      MRD   = (sum_BM-final_state(1))/sum_BM
      if (MRD>0.05) then
            error_praxis	= 100.0D0
            goto 197
      endif
*-----------------------------------------------Check the BM cellularity
C      BM_cellularity    = 100.0D0*sum_BM/full_BM_sum
C      if (BM_cellularity<15) then
C            error_praxis	= 100.0D0
C            goto 197
C      elseif (BM_cellularity>20) then
C            error_praxis	= 100.0D0
C            goto 197
C      endif
*-----------Compute the population at relapse according to the ODE model
      call deterministic(final_state,para)
*---------------------------------Compute the error of the parameter set
      sum_BM		= final_state(1)
      do i=1,N_clones
      	sum_BM	= sum_BM+final_state(2*i+1)
      enddo
      rel_clonal_com(N_clones+1)	= 100.0D0*final_state(1)/sum_BM
      do i=1,N_clones
      	rel_clonal_com(i)		= 100.0D0*final_state(2*i+1)/sum_BM
     	enddo
     	difference		= abs(rel_clonal_com-rel_clonal)
     	error_praxis	= maxval(difference)
197   END FUNCTION
*======================================PARAMETER FITTING BY USING PRAXIS
      SUBROUTINE PARAMETER_PRAXIS
      implicit none
*	Parameters for the algorithm
      integer,parameter		:: N = 2*N_clones
      double precision		:: fmin
      double precision		:: para_in_fitting(1:2*N_clones)
      double precision		:: para_initial(1:2*N_clones)
      double precision		:: para_final(1:2*N_clones)
      double precision		:: result,PRAXIS
      double precision		:: H0(1:2*N_clones)

      double precision		:: K

      para_initial(1)	= 0.8
      para_initial(2)	= 0.8
      para_initial(3)	= 0.5
      para_initial(4)	= 0.5

      para_in_fitting	= para_initial

      H0(1)				= 1.0D-3
      H0(2)				= 1.0D-3
      H0(3)				= 1.0D-3
      H0(4)				= 1.0D-3

      print*,'---------------------------------------------------------'
      print*,'Initial guess:'
      print*,para_initial
      result = PRAXIS(T0,MACHEP,H0,N,PRIN,
     .		    para_in_fitting,error_praxis,fmin)
      print*,'Final guess:'
      print*,para_in_fitting

      END SUBROUTINE
*========COMPUTE THE ERROR OF A PARAMETER SET, TO BE MINIMIZED BY NELMIN
      FUNCTION error_nelmin(para_set)
      implicit none
*	Input: parameters for the cancer clones, and number of variables
*	in the parameter set (always 2*N_clones, for a_l and p_l)
      double precision			:: para_set(1:2*N_clones)
*	Output: infinity-norm error of the percentages at relapse compared
*	to the real percentages
      double precision			:: error_nelmin
*	Parameters for the algorithm
      integer				:: N
*----------------------------------Use error_praxis to compute the error
      N     = 2*N_clones
      error_nelmin	= error_praxis(para_set,N)
      END FUNCTION
*======================================PARAMETER FITTING BY USING NELMIN
      SUBROUTINE PARAMETER_NELMIN
      implicit none
*	Output: No direct output, but all the parameter sets found here
*	are printed out in a *.txt file

*	Error indicator of NELMIN. 0 means no error, 1 means REQMIN, N, or
*	KONVGE has an illegal value, 2 means iteration terminated because
*	KCOUNT was exceeded without convergence.
      integer			:: ifault
*	Parameters for the algorithm
      integer,parameter		:: N = 2*N_clones
      double precision		:: para_initial(1:2*N_clones)
      double precision		:: para_final(1:2*N_clones)
      double precision		:: step(1:2*N_clones)
      double precision		:: ynewlo
      integer			:: icount,numres,i,index
      double precision		:: current_error
      double precision		:: ZBQLU01,ran
      integer			:: n_parameters_found
*-------------------Determines the size and shape of the initial simplex
      step(1)		= 0.0001
      step(2)		= 0.0001
      step(3)		= 0.0001
      step(4)		= 0.0001
*-------------------------------------------Use NELMIN to fit parameters
      n_parameters_found	= 0
      index				= 0
      open(UNIT=14,FILE='parameters_'//str_ID//'.txt',STATUS='replace')
      do while (n_parameters_found.lt.N_TOTAL_PARAMETERS)
*		Randomize the initial guess for the parameter set
      	call ZBQLINI(0)
      	do i=1,N_clones
      		ran	= ZBQLU01(1)
      			para_initial(i)	= min_a_l_nelmin+
     .				  	ran*(max_a_l_nelmin-min_a_l_nelmin)
      		ran	= ZBQLU01(2)
      		para_initial(i+N_clones)	= min_p_l_nelmin+
     .				  	ran*(max_p_l_nelmin-min_p_l_nelmin)
      	enddo
*		Repeat NELMIN again and again until error can't be lower
      	current_error	= error_nelmin(para_initial)
297         call NELMIN(error_nelmin,N,para_initial,para_final,ynewlo,
     .			REQMIN,step,KONVGE,KCOUNT,icount,numres,ifault)
     		if (ifault.ne.0) then
     			print*,'NELMIN returns ifault ',ifault
     			continue
     		endif
      	if (ynewlo.lt.(current_error-1.0D-1)) then
      		current_error	= ynewlo
      		goto 297
      	endif
      	print*,n_parameters_found,index,error_nelmin(para_final)
      	index	= index+1
*		In case another parameter set is found
      	if (error_nelmin(para_final).lt.1.0D0) then
      		n_parameters_found	= n_parameters_found+1
      		write(UNIT=14,FMT=*) para_final
      	endif
      enddo
      close(UNIT=14)
      END SUBROUTINE
*====================FIND THE PERCENTAGES AT RELAPSE FOR A PARAMETER SET
      SUBROUTINE percentage_relapse(para_set,perc)
      implicit none
*	Input: parameters for the cancer clones, and number of variables
*	in the parameter set (always 2*N_clones, for a_l and p_l)
      double precision			:: para_set(1:2*N_clones)
*	Output: percentages at relapse of the mitotic compartments of each
*	cancer clone and the normal clone (the last number)
      double precision			:: perc(1:N_clones+1)
*	Parameters for the algorithm
      double precision			:: final_state(N_total)
      double precision			:: para(1:3,1:N_clones)
      integer				:: i
      double precision			:: sum_BM
*---------------Compute the state at relapse for the given parameter set
      do i=1,N_clones
*		Renewal rates a_l
      	para(1,i)	= para_set(i)
*		Proliferation rates p_l
      	para(2,i)	= para_set(i+N_clones)
*		Death rates d_l
      	para(3,i)	= CANCER_DEATH_RATE
      enddo
      call deterministic(final_state,para)
*---------------Compute the percentages at relapse for the parameter set
      sum_BM		= final_state(1)
      do i=1,N_clones
      	sum_BM	= sum_BM+final_state(2*i+1)
      enddo
      perc(N_clones+1)	= 100.0D0*final_state(1)/sum_BM
      do i=1,N_clones
      	perc(i)		= 100.0D0*final_state(2*i+1)/sum_BM
     	enddo
      END SUBROUTINE
*======================================PARAMETER FITTING BY USING NELMIN
*=======================IN THE CASE WHERE THE RESULT IS OVERLY SENSITIVE
*======================================================TO THE PARAMETERS
      SUBROUTINE PARAMETER_NELMIN_SENSITIVE
      implicit none
*	Output: No direct output, but all the parameter sets found here
*	are printed out in a *.txt file

*	Error indicator of NELMIN. 0 means no error, 1 means REQMIN, N, or
*	KONVGE has an illegal value, 2 means iteration terminated because
*	KCOUNT was exceeded without convergence.
      integer			:: ifault
*	Parameters for the algorithm
      integer,parameter		:: N = 2*N_clones
      double precision		:: para_initial(1:2*N_clones)
      double precision		:: para_final(1:2*N_clones)
      double precision		:: step(1:2*N_clones)
      double precision		:: ynewlo
      integer			:: icount,numres,i,index
      double precision		:: current_error
      double precision		:: ZBQLU01,ran
      integer			:: n_parameters_found

      double precision		:: change_a_l,change_p_l
      double precision		:: rel_clonal_comp(1:N_clones+1)
      integer			:: clonal_marker(1:N_clones+1)
      integer			:: no_correction
*--------------------Determine the size and shape of the initial simplex
      step(1)		= 0.01
      step(2)		= 0.01
      step(3)		= 0.01
      step(4)		= 0.01
*---------------------Determine the step for educating the initial guess
      change_a_l		= 0.0005
      change_p_l		= 0.0005
*-------------------------------------------Use NELMIN to fit parameters
      n_parameters_found	= 0
      index				= 0
      no_correction		= 0
      open(UNIT=14,FILE='parameters_'//str_ID//'.txt',STATUS='replace')
      do while (n_parameters_found.lt.N_TOTAL_PARAMETERS)
*		Randomize the initial guess for the parameter set
      	call ZBQLINI(0)
      	do i=1,N_clones
      		ran	= ZBQLU01(1)
      			para_initial(i)	= min_a_l_nelmin+
     .				  	ran*(max_a_l_nelmin-min_a_l_nelmin)
      		ran	= ZBQLU01(2)
      		para_initial(i+N_clones)	= min_p_l_nelmin+
     .				  	ran*(max_p_l_nelmin-min_p_l_nelmin)
      	enddo
*		Repeat NELMIN again and again until error can't be lower
296         current_error	= error_nelmin(para_initial)
297         call NELMIN(error_nelmin,N,para_initial,para_final,ynewlo,
     .			REQMIN,step,KONVGE,KCOUNT,icount,numres,ifault)
     		if (ifault.ne.0) then
C     			print*,'NELMIN returns ifault ',ifault
     			continue
     		endif
      	if (ynewlo.lt.(current_error-1.0D-1)) then
      		current_error	= ynewlo
      		goto 297
      	endif
      	print*,n_parameters_found,index,
     .		 no_correction,error_nelmin(para_final)
      	index	= index+1
*		Check the result
      	if (error_nelmin(para_final).lt.1.0D0) then
*			In case another parameter set is found
      		n_parameters_found	= n_parameters_found+1
      		write(UNIT=14,FMT=*) para_final
      		no_correction		= 0
      	elseif ((no_correction.gt.BOUND_CORRECTION).and.
     .		  (current_error.ge.50.0D0)) then
      		no_correction	= 0
      	elseif ((no_correction.le.BOUND_CORRECTION).and.
     .		  (current_error.lt.50.0D0)) then
*			Make a more educated initial guess for NELMIN
*			first find the % at relapse for the current guess
      		call percentage_relapse(para_initial,rel_clonal_comp)
      		print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      		print*,rel_clonal_comp
*			Mark the clones whose percentages are smaller/larger
*			than the real percentages
      		do i=1,N_clones+1
      			if (abs(rel_clonal_comp(i)-rel_clonal(i))
     .			    .lt.0.5D0) then
      					clonal_marker(i)	= -1
      			elseif(rel_clonal_comp(i).lt.rel_clonal(i)) then
      				clonal_marker(i)	= 0
      			else
      				clonal_marker(i)	= 1
      			endif
      		enddo
*			Correct the initial guess
      		do i=1,N_clones
      			if (clonal_marker(i).eq.0) then
					para_initial(i+N_clones) =
     .						para_initial(i+N_clones)+
     .						change_p_l
     				elseif (clonal_marker(i).eq.1) then
					para_initial(i+N_clones) =
     .						para_initial(i+N_clones)-
     .						change_p_l
     				endif
      		enddo
      		no_correction	= no_correction+1
      		goto 296
      	endif
      enddo
      close(UNIT=14)
      END SUBROUTINE
*==================================================TAU-LEAPING TIME STEP
*===========================================FOR THE STOCHASTIC ALGORITHM
      SUBROUTINE tau_step(tau_time,vec_major,len_major,
     .                    para_p_l,para_d_l,D,Time,epsilon)
      implicit none
*     Input: vector of major species, parameter sets of leukemic clones,
*     current time, and error control for hybrid tau-leaping algorithm
      integer,dimension(2*N_clones+2)           :: vec_major
      integer                                   :: len_major
      double precision,dimension(N_clones)      :: para_p_l
      double precision,dimension(N_clones)      :: para_d_l
      double precision                          :: D
      double precision                          :: Time
      double precision                          :: epsilon
*     Output: tau-leaping time step
      double precision        :: tau_time
*     Parameters for the algorithm
      integer                 :: i,j,k,species
      double precision        :: rate1,rate2,k_nor,k_leu,p_l,d_l
*---------------------------------------------Find tau-leaping time step
      tau_time    = -1.0D0
      k_nor       = 1.0D0
      k_leu       = 1.0D0
      do i=1,CHEMO_TOTAL
            if((Time.ge.CHEMO_TIME(i,1)).and.
     .         (Time.le.CHEMO_TIME(i,2))) then
                  k_nor       = CHEMO_STRE(i,1)
                  k_leu       = CHEMO_STRE(i,2)
            endif
      enddo
      do i=1,len_major
            species     = vec_major(i)
*           Normal mitotic cells
            if(species.eq.1) then
                  rate1 = 1/D
                  rate2 = 1/(p_c*k_nor)
                  if(tau_time.le.0.0D0) then
                        tau_time    = min(rate1,rate2)
                  else
                        tau_time    = min(tau_time,rate1,rate2)
                  endif
*           Normal mature cells
            elseif(species.eq.2) then
                  rate1 = 1/d_c
                  if(tau_time.le.0.0D0) then
                        tau_time    = rate1
                  else
                        tau_time    = min(tau_time,rate1)
                  endif
            else
            do j=1,N_clones
                  p_l   = para_p_l(j)
                  d_l   = para_d_l(j)
*                 Leukemic mitotic cells
                  if(species.eq.2*j+1) then
                        rate1 = 1/D
                        rate2 = 1/(p_l*k_leu)
                        if(tau_time.le.0.0D0) then
                              tau_time    = min(rate1,rate2)
                        else
                              tau_time    = min(tau_time,rate1,rate2)
                        endif
*                 Leukemic mature cells in blood
                  elseif(species.eq.2*j+2) then
                        rate1 = 1/d_l
                        if(tau_time.le.0.0D0) then
                              tau_time    = rate1
                        else
                              tau_time    = min(tau_time,rate1)
                        endif
                  endif
            enddo
            endif
      enddo
      tau_time    = epsilon*tau_time
      END SUBROUTINE
*=================================STOICHIOMETRIC MATRIX AND INDEX MATRIX
*===========================================FOR THE STOCHASTIC ALGORITHM
      SUBROUTINE stoichiometric_matrix(Stoi,Index_mat)
      implicit none
*     Output: stoichiometric matrix, index matrix (1st row: reactions,
*     2nd row: species, 3rd row: clones)
      integer,dimension(2*(N_clones+1),5*(N_clones+1))      :: Stoi
      integer,dimension(3,5*(N_clones+1))                   :: Index_mat
*     Parameters for the algorithm
      integer     :: column,i,j
*-----------------------Build the stoichiometric matrix and index matrix
      do i=1,2*(N_clones+1)
            do j=1,5*(N_clones+1)
                  Stoi(i,j)         = 0
            enddo
      enddo
      do i=1,3
            do j=1,5*(N_clones+1)
                  Index_mat(i,j)    = 0
            enddo
      enddo
*     Asymmetric division of normal mitotic cells
      Stoi(1,1)               = 1
      Index_mat(1,1) = 1; Index_mat(2,1) = 1; Index_mat(3,1) = 0
      Stoi(2,2)               = 1
      Index_mat(1,2) = 2; Index_mat(2,2) = 1; Index_mat(3,2) = 0
      Stoi(1,3)               = -1
      Stoi(2,3)               = 2
      Index_mat(1,3) = 3; Index_mat(2,3) = 1; Index_mat(3,3) = 0
*     Death of normal mitotic cells and mature cells
      Stoi(1,4)               = -1
      Index_mat(1,4) = 4; Index_mat(2,4) = 1; Index_mat(3,4) = 0
      Stoi(2,5)               = -1
      Index_mat(1,5) = 5; Index_mat(2,5) = 2; Index_mat(3,5) = 0
*     Of leukemic cells...
      column                  = 6
      do i=1,N_clones
*           Asymmetric division of leukemic mitotic cells
            Stoi(2*i+1,column)      = 1
            Index_mat(1,column)     = column
            Index_mat(2,column)     = 2*i+1
            Index_mat(3,column)     = i
            column                  = column+1
            Stoi(2*i+2,column)      = 1
            Index_mat(1,column)     = column
            Index_mat(2,column)     = 2*i+1
            Index_mat(3,column)     = i
            column                  = column+1
            Stoi(2*i+1,column)      = -1
            Stoi(2*i+2,column)      = 2
            Index_mat(1,column)     = column
            Index_mat(2,column)     = 2*i+1
            Index_mat(3,column)     = i
            column                  = column+1
*           Death of leukemic mitotic cells
            Stoi(2*i+1,column)      = -1
            Index_mat(1,column)     = column
            Index_mat(2,column)     = 2*i+1
            Index_mat(3,column)     = i
            column                  = column+1
*           Death of leukemic mature cells in BM and in blood
            Stoi(2*i+2,column)      = -1
            Index_mat(1,column)     = column
            Index_mat(2,column)     = 2*i+2
            Index_mat(3,column)     = i
            column                  = column+1
      enddo
      END SUBROUTINE
*======================================================PROPENSITY VECTOR
*===========================================FOR THE STOCHASTIC ALGORITHM
      SUBROUTINE propensity(P,Time,vec_propensity,S,D,
     .                      para_a_l,para_p_l,para_d_l)
      implicit none
*     Input: current population vector, time, feedback parameters, and
*     leukemic parameters
      integer*8,dimension(2*(N_clones+1))       :: P
      double precision                          :: Time
      double precision                          :: S,D
      double precision,dimension(N_clones)      :: para_a_l
      double precision,dimension(N_clones)      :: para_p_l
      double precision,dimension(N_clones)      :: para_d_l
*     Output: propensity vector
      double precision,dimension(5*(N_clones+1)):: vec_propensity
*     Parameters for the algorithm
      integer           :: pos,i,j
      double precision  :: a_l,p_l,d_l,k_nor,k_leu
*----------------------------Build the propensity vector for normal time
*     Propensities for reactions involving the normal clone
      vec_propensity(1) = P(1)*p_c*(a_c*S)**2
      vec_propensity(2) = P(1)*p_c*2.0D0*(a_c*S)*(1-a_c*S)
      vec_propensity(3) = P(1)*p_c*(1-a_c*S)**2
      vec_propensity(4) = P(1)*D
      vec_propensity(5) = P(2)*d_c
*     Propensities for reactions involving the leukemic clones
      pos         = 6
      do j=1,N_clones
            a_l   = para_a_l(j)
            p_l   = para_p_l(j)
            d_l   = para_d_l(j)
            vec_propensity(pos)   = P(2*j+1)*p_l*(a_l*S)**2
            pos   = pos+1
            vec_propensity(pos)   = P(2*j+1)*p_l*2.0D0*(a_l*S)*(1-a_l*S)
            pos   = pos+1
            vec_propensity(pos)   = P(2*j+1)*p_l*(1-a_l*S)**2
            pos   = pos+1
            vec_propensity(pos)   = P(2*j+1)*D
            pos   = pos+1
            vec_propensity(pos)   = P(2*j+2)*d_l
            pos   = pos+1
      enddo
*--------------------Change the propensity vector if during chemotherapy
      do i=1,CHEMO_TOTAL
            if((Time.ge.CHEMO_TIME(i,1)).and.
     .         (Time.le.CHEMO_TIME(i,2))) then
                  k_nor       = CHEMO_STRE(i,1)
                  k_leu       = CHEMO_STRE(i,2)
                  vec_propensity(4) = vec_propensity(4)+P(1)*k_nor*p_c
                  do j=1,N_clones
                        p_l   = para_p_l(j)
                        vec_propensity(5*j+4) =
     .                      vec_propensity(5*j+4) + P(2*j+1)*k_leu*p_l
                  enddo
            endif
      enddo
      END SUBROUTINE
*================================================STOCHASTIC TRAJECTORIES
*=======================USING A HYBRID SSA - TAU LEAPING - ODE ALGORITHM
      SUBROUTINE STOCHASTIC(N_Traj,minsavestep)
      implicit none
*     Input: number of trajectories per parameter set, and minimal time
*     step for storing
      integer                       :: N_Traj
      double precision              :: minsavestep
*     Parameters for the algorithm
*     *---threshold between major and minor species
      integer*8,parameter           :: N_C=100
*     *---error control for hybrid tau-leaping algorithm
      double precision,parameter    :: epsilon=0.02D0
*     *---population vector at diagnosis
      integer*8,dimension(2*(N_clones+1))       :: P0
*     *---population and propensity vectors and time during integration
      integer*8,dimension(2*(N_clones+1))       :: P,P_temp,P_old
      double precision,dimension(5*(N_clones+1)):: vec_propensity
      double precision                          :: Time,SavedTime
*     *---integer versions of the patient parameters
      integer*8                     :: C_1_BAR_int,diag_BM_sum_int
      integer*8                     :: diag_normal_int,diag_blast_int
*     *---patient parameter file name
      character(len=30)             :: str_file
*     *---number of parameter sets
      integer,parameter             :: NO_PARAMETERS=3
*     *---matrix of all parameter sets
      double precision,dimension(NO_PARAMETERS,2*N_clones)
     .                              :: ALL_PARAMETERS
*     *---stoichiometric matrix and index matrix (1st row: reactions,
*     2nd row: species, 3rd row: clones)
      integer,dimension(2*(N_clones+1),5*(N_clones+1)):: Stoi
      integer,dimension(3,5*(N_clones+1))             :: Index_mat
*     *---a specific parameter set
      double precision,dimension(N_clones)      :: para_a_l
      double precision,dimension(N_clones)      :: para_p_l
      double precision,dimension(N_clones)      :: para_d_l
*     *---others
      integer                       :: i,j,k,column,index_parameter
      integer*8                     :: Sum_blood,Sum_BM
      integer,dimension(2*(N_clones+1))         :: vec_minor,vec_major
      double precision,dimension(5*(N_clones+1)):: a_mj,a_mn
      integer                       :: len_minor,len_major,species
      integer                       :: reaction
      double precision              :: S,D,tau_time,e_time,Timestep
      double precision              :: sum_mn,sum_SSA,tic,toc,clock
      double precision              :: ZBQLEXP,ZBQLNOR,ZBQLU01,lambda,r1
      integer*8                     :: r_num,ZBQLPOI
      character(len=10)             :: str_1,str_2,str_3,str_4
      character(len=3)              :: str_temp_1,str_temp_2


      double precision              :: PB_blast,PB_blast_old
      double precision              :: T_RELAPSE_computed
*----------------------------Change all cell count parameters to integer
      C_1_BAR_int       = nint(C_1_BAR,8)
      diag_BM_sum_int   = nint(diag_BM_sum,8)
      diag_normal_int   = nint(diag_normal,8)
      diag_blast_int    = nint(diag_blast,8)
*--------------------------------Find the population vector at diagnosis
      P0(1)            = nint(diag_BM_sum*diag_clonal(N_clones+1)/100,8)
      P0(2)            = diag_normal_int
      do j=1,N_clones
            P0(2*j+1)  = nint(diag_BM_sum*diag_clonal(j)/100,8)
            P0(2*j+2)  = nint(diag_blast*diag_clonal(j)/100,8)
      enddo
*-----------------------------------------------Input all parameter sets
      do j=1,NO_PARAMETERS
            if (j==1) then
                  str_file      = 'parameters_'//str_ID//'_left.txt'
            elseif (j==2) then
                  str_file      = 'parameters_'//str_ID//'_centered.txt'
            else
                  str_file      = 'parameters_'//str_ID//'_right.txt'
            endif
            open(unit=2,file=str_file,status='old',action='read')
            read(2,*) ALL_PARAMETERS(j,1:2*N_clones)
            close(unit=2)
      enddo
*-----------------------Build the stoichiometric matrix and index matrix
      call stoichiometric_matrix(Stoi,Index_mat)
*-----------------Loop computing the trajectories for all parameter sets
      do index_parameter=1,NO_PARAMETERS
*           Set up the specific leukemic clones' parameters
            para_a_l    = ALL_PARAMETERS(index_parameter,1:N_clones)
            para_p_l    = ALL_PARAMETERS(index_parameter,
     .                                            N_clones+1:2*N_clones)
            do i=1,N_clones
                  para_d_l    = CANCER_DEATH_RATE
            enddo
*           Set up the random number generator
            call ZBQLINI(0)
*           Loop computing the trajectories for this parameter set
            do i=1,N_Traj
                  PB_blast_old      = -100.0D0
                  write(str_temp_1,'(I3)') index_parameter
                  write(str_temp_2,'(I3)') i
                  print*
                  print*,'Patient ',str_ID,
     .            ' --- parameter set ',str_temp_1,
     .            ' --- trajectory ',str_temp_2
*                 Set up the initial data
                  P           = P0
                  Time        = 0.0D0
*                 Open file for storing
                  if(index_parameter.lt.10) then
                        str_2 = '(I1)'
                  elseif(index_parameter.lt.100) then
                        str_2 = '(I2)'
                  elseif(index_parameter.lt.1000) then
                        str_2 = '(I3)'
                  elseif(index_parameter.lt.10000) then
                        str_2 = '(I4)'
                  else
                        str_2 = '(I5)'
                  endif
                  write(str_1,str_2) index_parameter
                  if(i.lt.10) then
                        str_4 = '(I1)'
                  elseif(i.lt.100) then
                        str_4 = '(I2)'
                  elseif(i.lt.1000) then
                        str_4 = '(I3)'
                  elseif(i.lt.10000) then
                        str_4 = '(I4)'
                  else
                        str_4 = '(I5)'
                  endif
                  write(str_3,str_4) i
                  str_file    = 'stochastic_'//str_ID//'_'//trim(str_1)
     .                          //'_'//trim(str_3)//'.txt'

                  open(unit=4,file=str_file,action='write')
*                 Store the initial point
                  write(unit=4,fmt=*) Time,P
                  SavedTime   = Time
*-----------------Main program for one SSA
                  tic   = clock()

                  T_RELAPSE_computed      = T_RELAPSE

100               do while(Time.lt.T_RELAPSE_computed)
*                       Set up parameters
                        Sum_blood    = P(2)
                        do j=1,N_clones
                              Sum_blood    = Sum_blood+P(2*j+2)
                        enddo
                        Sum_BM = sum(P)-Sum_blood
                        S      = 1.0D0/(1.0D0+k_feedback*Sum_blood)
                        D      = C1_feedback*max(0.0D0,
     .                                       Sum_BM-C2_feedback*C_1_BAR)
*                       Find indices of minor and major species
                        vec_minor(1)      = 0;
                        len_minor         = 0;
                        vec_major(1)      = 0;
                        len_major         = 0;
                        do j=1,2*N_clones+2
                              if(P(j).lt.N_C) then
                                    len_minor   = len_minor+1
                                    vec_minor(len_minor)    = j
                              else
                                    len_major   = len_major+1
                                    vec_major(len_major)    = j
                              endif
                        enddo
*                       Find the propensities of all reactions
                        call propensity(P,Time,vec_propensity,
     .                                  S,D,para_a_l,para_p_l,para_d_l)
*-----------------------Find the time step
*                       Find the propensities of major/minor reactions
                        do k=1,5*(N_clones+1)
                              a_mj(k)     = 0.0D0
                              a_mn(k)     = 0.0D0
                        enddo
                        do j=1,len_major
                              species     = vec_major(j)
                              do k=1,5*(N_clones+1)
                                    if(Index_mat(2,k).eq.species) then
                                          a_mj(k)    = vec_propensity(k)
                                    endif
                              enddo
                        enddo
                        do j=1,len_minor
                              species     = vec_minor(j)
                              do k=1,5*(N_clones+1)
                                    if(Index_mat(2,k).eq.species) then
                                          a_mn(k)    = vec_propensity(k)
                                    endif
                              enddo
                        enddo
*                       Time step for major reactions
                        call tau_step(tau_time,vec_major,len_major,
     .                                para_p_l,para_d_l,D,Time,epsilon)
*                       Time step for minor reactions
                        e_time      = ZBQLEXP(1/sum(a_mn))
*                       Real time step
                        Timestep    = min(tau_time,e_time)
                        if(Timestep.ge.minsavestep) then
                              Timestep    = minsavestep
                        endif
*                       Store data if the end point has been reached
                        if((Time+Timestep).ge.T_RELAPSE_computed) then
                              Time  = T_RELAPSE_computed
                              write(unit=4,fmt=*) Time,P
                        endif
*-----------------------Simulate major reactions
101                     continue
*                       Initialize the temporary population vector
                        do j=1,2*(N_clones+1)
                              P_temp(j)   = P(j)
                              P_old(j)    = P(j)
                        enddo
*                       For every reactions...
                        do j=1,5*(N_clones+1)
*                             use either Poisson distribution (if lambda
*                             is small) or normal distribution (if
*                             lambda is big) to find number of reactions
                              lambda      = a_mj(j)*Timestep
                              if(lambda.le.1.0D8) then
                                    r_num = ZBQLPOI(lambda)
                              else
                                    r_num = ZBQLNOR(lambda,sqrt(lambda))
                              endif
*                             update the temporary population vector
                              do k=1,2*(N_clones+1)
                                    P_temp(k) =P_temp(k)+r_num*Stoi(k,j)
                              enddo
                        enddo
*                       Redo tau-leaping if any species is negative
                        if(minval(P_temp).lt.0) then
                              Timestep    = 0.5D0*Timestep
                              goto 101
                        else
                              P     = P_temp
                        endif
*-----------------------Simulate minor reactions
                        if(e_time.le.Timestep) then
*                             Find the reaction that happened
                              sum_mn      = sum(a_mn)
                              r1          = ZBQLU01(1)
                              reaction    = 1
                              do j=1,5*(N_clones+1)
                                    sum_SSA     = sum(a_mn(1:j))
                                    if(r1.gt.sum_SSA/sum_mn) then
                                          reaction    = reaction+1
                                    endif
                              enddo
*                             Update the population vector
                              do k=1,2*(N_clones+1)
                                    P(k)  = P(k)+Stoi(k,reaction)
                              enddo
*                             Redo tau-leaping and SSA if any species is
*                             negative
                              if(minval(P).lt.0) then
                                    Timestep    = 0.5D0*Timestep
                                    do k=1,2*(N_clones+1)
                                          P(k)  = P_old(k)
                                    enddo
                                    goto 101
                              endif
                        endif
*                       Store the data
                        Time  = Time+Timestep
                        if(Time.gt.SavedTime+minsavestep) then
                              SavedTime   = Time
                              write(unit=4,fmt=*) Time,P
                        endif
                  enddo
*-----------------Compute the MRD at the end of simulation and continue
*-----------------if MRD<5%
*                 Compute the MRD
                  PB_blast    = 0.0D0
                  do j=1,N_clones
                        PB_blast    = PB_blast+1.0D0*P(2*j+2)
                  enddo
                  PB_blast    = 100.0D0*PB_blast/(PB_blast+1.0D0*P(2));
                  print*,'MRD at ',T_RELAPSE_computed,' is ',PB_blast
*                 If continued simulation doesn't change MRD, stop
                  if(abs(PB_blast_old-PB_blast)<1.0D-3) then
                        goto 102
                  endif
*                 Continue simulation if MRD<5%
                  if(PB_blast.lt.5.0D0) then
                        T_RELAPSE_computed = T_RELAPSE_computed+100.0D0
                        PB_blast_old       = PB_blast
                        goto 100
                  endif
*                 End simulation
102               close(unit=4)
                  toc   = clock()
                  print*,'     time cost = ',toc-tic
            enddo
      enddo
      END SUBROUTINE
*=======================================================================
      END MODULE
