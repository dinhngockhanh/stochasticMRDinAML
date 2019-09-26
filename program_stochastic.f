      PROGRAM program_stochastic
      use module_grand
      implicit none
*-----------------------------------------------------Control parameters
*     Number of trajectories per parameter set to produce
      integer,parameter             :: N_Traj = 1000
*     Minimal time step for storing
      double precision,parameter    :: minsavestep = 1.0D0
*------------------------------------Compute the stochastic trajectories
      call STOCHASTIC(N_Traj,minsavestep)

      END
