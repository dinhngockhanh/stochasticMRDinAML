*	Important data from each patient
      MODULE module_patient
      implicit none
*	TCGA ID of the patient
      character(len=6),parameter	:: str_ID 	= '573988'
*	Number of cancer clones
      integer,parameter			:: N_clones	= 2
*	Number of species
      integer,parameter			:: N_total	= 2 + N_clones*2
*	Equilibrium hematopoietic mitotic cell numbers
      double precision,parameter	:: weight 	= 80.7D0
      double precision,parameter	:: C_1_BAR 	= 2.0D0*(10.0**9)*weight
*	Total number of cells in bone marrow at primary/diagnosis
      double precision,parameter	:: p_cellularity	= 0.9D0
      double precision,parameter	:: full_BM_sum  = 4.6D0*(10.0D0**11)
      double precision,parameter	:: diag_BM_sum    = p_cellularity*
     .								  full_BM_sum
*	Number of cells in blood at primary/diagnosis
      double precision,parameter	:: diag_normal = 1.22D0*(10.0D0**10)
      double precision,parameter	:: diag_blast  = 7.60D0*(10.0D0**9)
*	Clonal percentages in bone marrow at diagnosis: (last element is
*	for normal clone)
      double precision,parameter	:: diag_clonal(1:N_clones+1)	=
     .				(/72.0D0,3.0D0,25.0D0/)
*	Clonal percentages in bone marrow at relapse: (last element is
*	for normal clone)
      double precision,parameter	:: rel_clonal(1:N_clones+1)	=
     .				(/15.0D0,39.0D0,46.0D0/)
*	Time and strength of chemotherapy
      integer,parameter			:: CHEMO_TOTAL	= 2
      double precision,parameter	:: CHEMO_TIME(1:CHEMO_TOTAL,1:2) =
     .		         reshape((/0.0D0,37.0D0,
     .					 7.0D0,42.0D0/),shape(CHEMO_TIME))
      double precision,parameter	:: CHEMO_STRE(1:CHEMO_TOTAL,1:2) =
     .		       reshape((/1*1.0D0,1*1.0D0,
     .				     5*1.0D0,5*1.0D0/),shape(CHEMO_TIME))
C      double precision,parameter	:: CHEMO_STRE(1:CHEMO_TOTAL,1:2) =
C     .		          reshape((/3.9D0,6.5D0,
C     .					  3.9D0,6.5D0/),shape(CHEMO_TIME))
C     	double precision,parameter	:: CHEMO_STRE(1:CHEMO_TOTAL,1:2) =
C     .		          reshape((/3.5D0,6.5D0,
C     .					  3.5D0,6.5D0/),shape(CHEMO_TIME))
*	Time until relapse
      double precision,parameter	:: T_RELAPSE	= 365.0D0
*----------Other important parameters that are the same for all patients
*	Parameters of hematopoietic clone:
      double precision,parameter	:: a_c = 0.87D0
      double precision,parameter	:: p_c = 0.45D0
      double precision,parameter	:: d_c = 2.3D0
*	Constant for negative feedback of hematopoietic clone, as in
*	feedback = 1/(1+k*c2) where c2 is number of post-mitotic
*	hematopoietic cells number:
      double precision,parameter	:: k_feedback	= 1.0D-12
*	Constants for negative feedback of overcrowding bone marrow, as in
*	feedback = C1*max(0,x-C2*c_1_d_bar) where x is total number of
*	cells in bone marrow (including mitotic hematopoietic and mitotic
*	cancer cells):
      double precision,parameter	:: C2_feedback 	= 1.0D0
      double precision,parameter	:: C1_feedback	= 1.0D-12

      save

      END MODULE
