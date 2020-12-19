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
!  It seeks for the first value in 'array' that matches the value 'val'
!  within a tolerance 'eps' and returns its index/position in the 'array'
        function find_index(length, array, val, ini, fin)
          integer*16 :: find_index, fin, i, length
          integer :: ini
          real*16 :: val, eps
          real*16, dimension (1:length) :: array
!
          eps = 0.000000001
          do i = ini, fin
            if (abs(array(i) - val) < eps) then
              find_index = i
              return
            else
              find_index = 0
            endif
          enddo
          return
        end function find_index
! -----------------------------------------------------------
!  INITIAL_T0 --------------------------------------
!  It generates an initial configurations of 'N' agents at consensus -1
!  which is a ground state for alpha<1
        subroutine initial_T0(N, &
                              array, n1, n_1, n0)
          integer*16 :: N, i , n1, n_1, n0
          integer*16, dimension(N) :: array
          do i = 1, N
            array(i) =  - 1    
          end do
          n_1 = N; n0 = 0; n1 = 0  ! Number of agents in each state
        end subroutine initial_T0
! -----------------------------------------------------------
!  INITIAL_inf --------------------------------------
!  It generates an initial configurations of 'N' agents randomly distributed
!  between the 3 states, which is compatible with infinite temperature
        subroutine initial_inf(N, &   ! Input
                              array, n1, n_1, n0)  ! Output
          integer*16 :: N, i, n1, n_1, n0
          integer*16, dimension(N) :: array
          do i = 1, N
            array(i) =  - 1
          end do
          n_1 = 0; n0 = 0; n1 = 0  ! Number of agents in each state
          do i = 1, N
            if (array(i) == -1) then
              n_1 = n_1 + 1
            else if (array(i) == 0) then
              n0 = n0 + 1
            else
              n1 = n1 + 1
            endif
          end do
        end subroutine initial_inf
!
! --------------------------------------------------------------------------
!  SYSTEM EVOLVES ----------------------------------------------------------
!  It proposes a single flip in one randomly chosen agent from 'array' and 
!  returns the new energy, and number of agents in each state associated 
!  to this change
        subroutine system_evolves(N, array, alpha2, old1, old_1, old0, & ! input
                                  new1, new_1, new0, DE, agent,new_spin, b) ! Output
          integer*16 :: N
          integer*16 :: agent, b, new_spin
          integer*16, dimension(N) :: array
          real*16 :: DE, C1, C2, alpha2
          integer*16 :: old1, old_1, old0, new1, new_1, new0
          real :: r1279
!
          C1 = -4.0/2.0
          C2 = -1.0/2.0
          agent = mod(int(N*r1279()), N) + 1  ! Pick a random spin
!   Change proposal and associated energy change 'DE' and new state populations ---------
            b = mod(int(2*r1279()), 2)   ! Change proposal
!            print*, agent, b
            if (array(agent) == -1) then
              if (b == 0) then
                DE=C2*(-2.0*float(old_1)+2.0*float(old1)+2.0*alpha2*float(old0)+2.0)
                new_spin = 0              !  -1 to 0
                new1=old1
                new0=old0+1
                new_1=old_1-1
              else
                DE=C1*float(1+old1-old_1)
                new_spin = 1              !  -1 to +1
                new1=old1+1
                new0=old0
                new_1=old_1-1
              endif
            else if (array(agent) == 0) then
              if (b == 0) then
                DE=C2*(-2.0*float(old1)+2.0*float(old_1)-2.0*alpha2*float(old0)+2.0*alpha2)
                new_spin = -1             ! 0 to -1
                new1=old1
                new0=old0-1
                new_1=old_1+1
              else
                DE=C2*(2.0*float(old1)-2.0*float(old_1)-2.0*alpha2*float(old0)+2.0*alpha2)
                new_spin = 1              !  0 to +1
                new0=old0-1
                new1=old1+1
                new_1=old_1
              endif
            else
              if (b == 0) then
                DE=C2*(2.0*float(old_1)-2.0*float(old1)+2.0*alpha2*float(old0)+2.0)
                new_spin = 0              ! +1 to 0
                new0=old0+1
                new1=old1-1
                new_1=old_1
              else
                DE=C1*float(1+old_1-old1)
                new_spin = -1             ! +1 to -1
                new1=old1-1
                new_1=old_1+1
                new0=old0
              endif
            endif
          end subroutine system_evolves
! --------------------------------------------------------------------------
! --- ACCEPTANCE ---------------------------------------
! --- It evaluates whether the proposed change is accepted or not, maximizing the probability
! --  to go to unvisited configurations (landau-wang) and it delivers the energy, the 
! --  number of agents in each state and the position/index of the state in the energy array
! --  and dos array at the end of the current step ------------------------
          subroutine acceptance(dos_dim,ln_dos,Efinal,Einitial,&
                                old1, old_1, old0, p1, p_1, p0, &
                                old_spin, p_spin, old_E_bin, p_E_bin, old_ene, diff, & ! Input
!
                                n1, n_1, n0, ene, spin, bin)  ! Output
!
          integer*16 :: old_E_bin, p_E_bin, bin, p_spin, old_spin, spin
          integer*16 :: dos_dim
          real*16, dimension(dos_dim) :: ln_dos 
          integer*16 :: p1, p_1, p0, n1, n_1, n0, old1, old_1, old0
          real*16 ene, old_ene, p_ene, diff, Efinal, Einitial, prob_flip
          real :: r1279 
          p_ene = old_ene + diff
          if ((p_ene <= Efinal).and.(p_ene >= Einitial)) then
            prob_flip = ln_dos(old_E_bin) - ln_dos(p_E_bin)
            if (prob_flip > 0.0) then
              ene = p_ene
              spin = p_spin    ! ACCEPT -----------------------------------
              n0=p0
              n1=p1
              n_1=p_1
              bin = p_E_bin
!              print*, 'A1'
            else if (r1279() < exp(prob_flip)) then
              ene = p_ene
              spin = p_spin    ! ACCEPT -------------------------------
              n0=p0
              n1=p1
              n_1=p_1
              bin = p_E_bin
!              print*, 'A2'
            else
              ene = old_ene
              spin = old_spin    ! REJECT -------------------------------
              n0=old0
              n1=old1
              n_1=old_1
              bin = old_E_bin
!              print*, 'R1'
            endif
          else  
            ene = old_ene
            spin = old_spin    ! REJECT -------------------------------
            n0=old0
            n1=old1
            n_1=old_1
            bin = old_E_bin
!            print*, 'R2'
          endif
          end subroutine
! -----------------------------------------------------------------------
      end module my_func
