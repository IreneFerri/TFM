!  
!   Functions Module for the Landau Wang algorithm on the alpha-3-states
!   model
!    
! -------------------------------------------------------------------------
!                              landau_wang.f90
! --------------------------------------------------------------------------
!
! Monte Carlo simulator of a system of spins that can take 3 values 
! corresponding to the vectors (-1,0),(0,alpha) and (1,0) 
! respectively, where alpha is a real parameter. Energy is defined as the
! negative sum of dot products over all pairs (fully connected graph).
! The program proposes sequencial changes and uses the histogram entropic
! Landau-Wang algorithm to visit all states and find the first order transition
! of the system for alpha above the tricritical point.
! After obtaining the states density array it uses it to calculate the
! canonical probability, the partition function and the average energy for 
! different temperatures and stores them in a file
! 
! -----------------------------------------------------------------------------
      module my_func
      implicit none
      contains

!  FIND INDEX  -----------------------------------------------------------
!  It seeks for the index corresponding to the tuple num1, num0, num_1
        function find_index(N, n0, n1)
          integer*16 :: find_index, N, n1, n0, counter, counter_2
!
          find_index = 0
          counter = N + 1
          counter_2 = n0
          do while (counter_2 > 0)
            find_index = find_index + counter ! Combinations of 'N-counter' agents in 2 states (2 degrees of freedom)
            counter = counter - 1
            counter_2 = counter_2 - 1
          enddo
          find_index = find_index + n1 + 1
          return
        end function find_index
! -----------------------------------------------------------
! --------------------------------------------------------------------------
!  SYSTEM EVOLVES ----------------------------------------------------------
!  It proposes a single flip in one randomly chosen agent and 
!  returns the new energy, and number of agents in each state associated 
!  to this change
        subroutine system_evolves(N, alpha2, old_1, old0, old1, & ! Input
                                  DE, new_1, new0, new1)   ! Output
          integer*16 :: N
          integer*16 :: agent, b
          real*16 :: DE, C1, C2, alpha2
          integer*16 :: old1, old_1, old0, new1, new_1, new0
          real :: r1279
!
          C1 = -4.0/2.0
          C2 = -1.0/2.0
          agent = mod(int(N*r1279()), N) + 1  ! Pick a random spin
          b = mod(int(2*r1279()), 2)   ! Change proposal
          if (agent <= old_1) then 
            if (b == 0) then  ! ----------------------------------------------------              
              DE=C2*(-2.0*float(old_1)+2.0*float(old1)+2.0*alpha2*float(old0)+2.0)  
              new_1 = old_1 - 1               !  -1 to 0
              new0 = old0 + 1
              new1 = old1
!              print*, '-1 to 0   A'
            else   ! ------------------------------------------------------------------
              DE=C1*float(1+old1-old_1)
              new_1 = old_1 - 1                !  -1 to +1
              new0 = old0
              new1 = old1 + 1
            endif
          else if (agent <= old_1 + old0) then
            if (b == 0) then   ! ----------------------------------------------------
              DE=C2*(-2.0*float(old1)+2.0*float(old_1)-2.0*alpha2*float(old0)+2.0*alpha2)
              new_1 = old_1 + 1
              new0 = old0 - 1                ! 0 to -1
              new1 = old1
            else   ! -------------------------------------------------------------------
              DE=C2*(2.0*float(old1)-2.0*float(old_1)-2.0*alpha2*float(old0)+2.0*alpha2)
              new_1 = old_1
              new0 = old0 - 1 
              new1 = old1 + 1                 ! 0 to +1
            endif
          else
            if (b == 0) then ! ----------------------------------------------------
              DE=C2*(2.0*float(old_1)-2.0*float(old1)+2.0*alpha2*float(old0)+2.0)
              new_1 = old_1
              new0 = old0 + 1                 ! +1 to 0
              new1 = old1 - 1
            else  !----------------------------------------------------------------
              DE=C1*float(1+old_1-old1)
              new_1 = old_1 + 1
              new0 = old0
              new1 = old1 - 1             ! +1 to -1
            endif
          endif
        end subroutine system_evolves
! --------------------------------------------------------------------------
! --- ACCEPTANCE ---------------------------------------
! --- It evaluates whether the proposed change is accepted or not, maximizing the probability
! --  to go to unvisited configurations (landau-wang) and it delivers the energy, the 
! --  number of agents in each state and the position/index of the state in the energy array
! --  and dos array at the end of the current step ------------------------
          subroutine acceptance(prob_flip, Efinal, Einitial,&
                                old1, old_1, old0, p1, p_1, p0, &
                                old_E_bin, p_E_bin, old_ene, diff, & ! Input
!
                                n1, n_1, n0, ene, bin)  ! Output
!
          integer*16 :: old_E_bin, p_E_bin, bin   
          integer*16 :: p1, p_1, p0, n1, n_1, n0, old1, old_1, old0
          real*16 ene, old_ene, p_ene, diff, Efinal, Einitial, prob_flip
          real :: r1279 
          p_ene = old_ene + diff
          if ((p_ene <= Efinal).and.(p_ene >= Einitial)) then
            if (prob_flip > 0.0) then
              ene = p_ene      ! ACCEPT -----------------------------------
              n0=p0
              n1=p1
              n_1=p_1
              bin = p_E_bin
            else if (r1279() < exp(prob_flip)) then
              ene = p_ene       ! ACCEPT -----------------------------------
              n0=p0
              n1=p1
              n_1=p_1
              bin = p_E_bin
            else
              ene = old_ene        ! REJECT -----------------------------------
              n0=old0
              n1=old1
              n_1=old_1
              bin = old_E_bin
            endif
          else  
            ene = old_ene       ! REJECT -----------------------------------
            n0=old0
            n1=old1
            n_1=old_1
            bin = old_E_bin
          endif
          end subroutine acceptance
! -----------------------------------------------------------------------
      end module my_func
