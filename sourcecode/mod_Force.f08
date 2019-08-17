module mod_Force
  !
  ! This module should contain both the overlap computation and the force computation in seperate functions. The overlap computation
  ! may reveal that there is no actual overlap between the particles at all, so then the force computation will not be necessary. In
  ! that case, working with independent functions allows us to pass less data through the cache than if both area and force were
  ! computed in one function. The force computation should also include the adding of forces and torques to the total force and total
  ! torque of the respective particles. The conversion from forces and torques to linear and angular accelerations may also be
  ! implemented in this module; the accelerations (instead of the forces) can then be passed to the predictor-corrector module.
  !
  !
  ! NOMENCLATURE
  ! A_o       - overlap area [m^2]
  ! c_o       - centroid of overlap volume [m]
  ! f         - sum of all forces [N]
  ! f_b       - buoyancy force (just one component -> z-direction) [N]
  ! f_coh     - cohesive force component [N]
  ! f_d_n     - dissipative force component in normal direction [N]
  ! f_e       - elastic force component
  ! f_n       - normal force (of the overlap polyhedron)
  ! f_norm    - normal force (acting on the overlapping elements)
  ! f_t       - tangential force
  ! f_t_old_p - projected and tangential frictional force
  ! f_t_old_r - resized tangentional frictional force
  ! g         - gravitational constant
  ! i         - Index of first element
  ! j         - Index of second elements
  ! k         - current pair
  ! L_c       - characteristic length
  ! M_red     - reduced mass
  ! mu        - friction coefficient
  ! C_d       - drag coefficient
  ! n         - force direction, pointing from element i to element j
  ! r1        - vector connecting center of overlap volume and center of element i
  ! r2        - vector connecting center of overlap volume and center of element j
  ! rho_w     - density of water
  ! t1, t2, clock_rate, clock_max - parameters used for time measuring
  ! v_contact - contact velocity
  ! V_o       - overlap volume
  ! v_n       - normal component of contact velocity
  ! v_t       - tangential part of contact velocity
  ! Y         - young's modulus
  ! Y_coh     - a fraction of Y used for cohesive force
  ! Y_t       - young's modulus for tangential force computation (usally Y_t = nu * Y, where nu = poisson's ratio)
  !
  implicit none
contains
  subroutine ForceComputation (coh_coeff, gamma_n, dt, n_elem, n_ice, n_struct, n_wall, option)
    !changelog:
    !
    use omp_lib
    use mod_Global
    use mod_OverlapComputation
    use mod_Functions
    use mod_Neighbour

    implicit none

    double precision, intent (in) :: coh_coeff, gamma_n
    double precision, intent (in) :: dt
    integer,          intent (in) :: n_elem, n_ice, n_struct, n_wall, option

    integer :: i, j, k
    double precision :: M_red, L_c
    double precision :: f_n(3), f_coh, f_e, f_d_n
    double precision :: f_t(3), f_t_old_p(3), f_t_old_r(3), f_d_t(3)
    double precision :: f(3)
    double precision, allocatable, save :: forces(:,:), torques(:,:)
    double precision :: Y, Y_t
    double precision :: mu, C_d, Re
    double precision :: v_contact(3), v_n(3), v_t(3)
    double precision :: n(3), r1(3), r2(3)
    double precision :: A_o, V_o, c_o(3)
    double precision :: y_coh
    double precision :: f_b
    double precision :: G  = 9.81 !m/s**2
    double precision :: RHO_W = 1025.0 !saltwater

    if (.not. allocated (forces)) then
      allocate (forces(n_elem,3))
      allocate (torques(n_elem,3))
    end if

    call CalculateNeighbours

    !======== PARALLEL START ==========!
    !$OMP PARALLEL DEFAULT (NONE) &
    !$OMP PRIVATE (i, j, k, A_o, c_o, V_o, r1, r2, &
    !$OMP          L_c, M_red, n, v_contact, v_n, v_t, Y, Y_t, f_e, f_d_n, f_d_t, f_n, &
    !$OMP          mu, C_d, Re, f_t_old_p, f_t_old_r, f_t, y_coh, f_coh, f, f_b) &
    !$OMP SHARED  (overlapPairs, elements, G, RHO_W, nPairs, n_elem, n_ice, n_struct, n_wall, &
    !$OMP          forces, torques, gamma_n, dt, ice, option, coh_coeff)

    forces(:,:)  = 0.0
    torques(:,:) = 0.0

    !$OMP DO SCHEDULE (DYNAMIC) REDUCTION (+: forces, torques)
    do k = 1, nPairs

      i = overlapPairs(k)%owner(1)
      j = overlapPairs(k)%owner(2)

      !i and j could be boundaries, therefore one has to cycle current loop if two boundaries would be checked for overlap
      if (i > n_ice .and. j > n_ice) cycle

      call OverlapComputation (i, j, k, n, A_o, c_o, V_o)

      contact_forces: if (V_o > 0.0) then

        r1 = c_o - elements(i)%r
        r2 = c_o - elements(j)%r

        ! caluclation of reduced mass (i.e. the effective inertial mass), characteristic length and the damping force, UtDEM p.246
        L_c     = 4.0 * norm2 (r1) * norm2 (r2) / ( norm2 (r1) + norm2 (r2) )
        M_red   = elements(i)%mass * elements(j)%mass / (elements(i)%mass + elements(j)%mass)

        ! To ensure to calculate physically correct friction between structure and elements, since the friction between ice-ice,
        ! ice-wall and ice-consol.layer is unphysically large to take cohesive bond into account.
        mu = min (elements(i)%mu, elements(j)%mu)

        v_contact = (elements(i)%rtd1 + Cross (elements(i)%om, r1)) - (elements(j)%rtd1 + Cross (elements(j)%om, r2))
        v_n       = (n / norm2 (n)) * (dot_product (v_contact, n) / norm2 (n))! a*b = |a| |b| cos(alpha)      
        v_t       = v_contact - v_n

        ! calculation of stiffness: Remember that the time step is calculated using the ice's Young's Modulus.
        !                           Therefore a contact with the structure would be too stiff and the time
        !                           step would be too big to resolve this contact correctly (at least for high
        !                           impact velocities).
        Y   = min (elements(i)%youngs, elements(j)%youngs)
        Y_t = 0.33 * Y !poisson ratio of ice

        !===============!
        ! NORMAL FORCES !
        !===============!

        !---------------!
        ! elastic force !
        !---------------!
        f_e = Y * V_o / L_c

        !---------------!
        ! damping force !
        !---------------!
        f_d_n = gamma_n * sqrt (Y * M_red / L_c ** 3) * (V_o - overlapPairs(k)%overlap_volume_old) / dt

        overlapPairs(k)%overlap_volume_old = V_o

        ! To avoid the jump in the sum of elastic and damping forces at separation,
        ! one could set the total force to zero if the damping would reverse the sign
        ! of the force when it is added to the elastic part (Understanding the DEM, page 227)
        if (f_e + f_d_n < 0.0) f_d_n = - f_e

        f_n = (f_e + f_d_n) * n

        !===================!
        ! TANGENTIAL FORCES !
        !===================!

        ! Cundall-Strack friction from "Understanding the DEM"
        !------------------------------------------------------
        ! 1) Projection onto the new tangential plane. During the advance from time t-tau to time t, with new contact normal n^(t) and
        !    new tangential velocity v_t(t), we project the old tangential force f_t(t-tau) onto the new tangential plane
        f_t_old_p = overlapPairs(k)%f_t_old - (dot_product (overlapPairs(k)%f_t_old, n)) * n

        ! 2) Rescaling the old magnitude. We then rescale f_t(t-tau)^p to the magnitude of the previous tangential force ||f_t(t-tau)||
        f_t_old_r = norm2 (overlapPairs(k)%f_t_old) * (f_t_old_p / norm2 (f_t_old_p))

        ! if there is no tangential force from the time step before, point 2) would lead to
        ! NaN because of the division by 0. Therefore there is no old magnitude to rescale and
        ! it is set to zero
        if (norm2 (f_t_old_p) < 10.0d-14) f_t_old_r = 0.0

        ! 3) Vectorial addition of the new increment. The rescaled projection f_t(t-tau)^r is then incremented to the new tangential force
        f_t = f_t_old_r - Y_t * L_c * v_t * dt

        ! 4) Application of a cut-off, if necessary. Finally, a cut-off is applied if the result from the previous vector addition
        !    exceeds the maximal friction allowed (the dynamic friction).
        if (norm2 (f_t) > mu * norm2 (f_n)) then 
          f_t = (f_t / norm2 (f_t)) * mu * norm2 (f_n)
        end if

        overlapPairs(k)%f_t_old = f_t

        !----------!
        ! cohesion !
        !----------!
        ! Cohesion is calcualted after the friction force, because the maximum friction depends on f_n without cohesion
        ! The cohesion is calculated between rubble-rubble and rubble-consolidated layer.
        if (i <= n_ice+n_wall .and. j <= n_ice+n_wall) then
          Y_coh = coh_coeff * Y
          f_coh = Y_coh * A_o
          f_n = f_n - f_coh * n
        end if

        !=============!
        ! TOTAL FORCE !
        !=============!

        f = f_n + f_t

        !====================================!
        ! NORMAL PARTS ACTING ON C_1 AND C_2 !
        !====================================!

        forces(i,:)  = forces(i,:)  + f
        torques(i,:) = torques(i,:) + Cross ( r1, f )

        ! actio = reactio (f = -f)
        forces(j,:)  = forces(j,:)  - f
        torques(j,:) = torques(j,:) + Cross ( r2, -f )

      else !no contact
        overlapPairs(k)%f_t_old = 0.0 ! otherwise the calculation of the tangential force at a new contact between i and j
        !                               would start with the last tangential force from previous contact
      end if contact_forces

    end do
    !$OMP END DO

    !$OMP DO
    do i = 1, n_ice
      ice(i)%force  = forces(i,:)
      ice(i)%torque = torques(i,:)
    end do
    !$OMP END DO

    ! adding bouancy, gravitational and viscous forces for each rubble element
    !$OMP DO
    do i = 1, n_ice
      Re = norm2 (ice(i)%rtd1) * sqrt (ice(i)%A) / 0.0000017926
      C_d = 24.0 / Re + (2.6 * (Re / 5.0)) / (1 + (Re / 5.0) ** 1.52) +              & 
        & (0.411 * (Re / 263000.0) ** (-7.94)) / (1 + (Re / 263000.0) ** (-8.0)) +  (Re ** 0.8 / 461000.0)
      !C_d = 24.0 / Re + 2.0
      ice(i)%force  = ice(i)%force  - 0.5 * RHO_W * ice(i)%rtd1 * norm2 (ice(i)%rtd1) * C_d * ice(i)%A
      ice(i)%torque = ice(i)%torque - 0.02 * C_d * ice(i)%om ! just an approxiamtion for rotative damping
      if (ice(i)%r(3) > 0.0) then
        f_b = RHO_W * (ice(i)%length * ice(i)%width * ice(i)%thickness) * G
        ice(i)%force(3) = ice(i)%force(3) + ice(i)%mass * G - f_b
      else
        ice(i)%force(3) = ice(i)%force(3) + ice(i)%mass * G
      end if
    end do
    !$OMP END DO

    !$OMP END PARALLEL
    !=================!

    ! "structure" is a global variable which components "force" and "torque" are
    ! used for output. Therefore, one has to assign the corresponding values in
    ! forces and torques back to structure
    do i = 1, n_struct
      structure(i)%force  = forces(n_ice+n_wall+i,:)! maybe + buoyancy??
      structure(i)%torque = torques(n_ice+n_wall+i,:)
    end do
  end subroutine ForceComputation
end module mod_Force
