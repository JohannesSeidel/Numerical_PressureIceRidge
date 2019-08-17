module mod_PredCorr
  !
  !This contains the centers of mass, orientations and velocities of the particles (in SI units), which should be accessible from the
  !main program. Internally, it will have another set of positions and orientations, as well as their higher derivatives, which are
  !scaled with time in dimensionless units for the predictor and the corrector subroutines (; see the remarks in §2.6). The
  !gravitational constant g = 9.81 m/s^2 can also be assigned in this module, as it is an acceleration, not a force
  !
contains
  subroutine Predictor (dt, n_ice)
    use omp_lib
    use mod_Global

    implicit none

    double precision, intent (in) :: dt
    integer,          intent (in) :: n_ice

    integer          :: i
    double precision :: c1, c2

    c1 = dt
    c2 = 0.5 * dt ** 2

    !$OMP PARALLEL DEFAULT (NONE)     &
    !$OMP PRIVATE (i)                 &
    !$OMP SHARED (c1, c2, ice, n_ice)

    !$OMP DO
    do i = 1, n_ice
      ! rtd* stands for r-time-derivative *, i.e. rtd2 is the second time derivative of the vector r.
      ice(i)%r    = ice(i)%r    + c1*ice(i)%rtd1 + c2*ice(i)%rtd2
      ice(i)%rtd1 = ice(i)%rtd1 + c1*ice(i)%rtd2

      ice(i)%q    = ice(i)%q    + c1*ice(i)%qtd1 + c2*ice(i)%qtd2
      ice(i)%qtd1 = ice(i)%qtd1 + c1*ice(i)%qtd2
    end do
    !$OMP END DO

    !$OMP END PARALLEL
  end subroutine Predictor

  subroutine Corrector (dt, n_ice, n_struct, option)
    use omp_lib
    use mod_Functions
    use mod_Global

    implicit none

    double precision, intent (in)  :: dt   
    integer,          intent (in)  :: n_ice, n_struct, option

    integer          :: i
    double precision :: kin
    double precision, save :: hertz

    double precision :: coeff0, coeff1, coeff2

    double precision :: rotmat(3,3), rotmat_trans(3,3), momentum(3)
    double precision :: J(3,3), J_b_inv(3,3), J_inv(3,3)
    double precision :: omtd1(3)
    double precision :: raccel(3), rcorr(3)
    type (qtrn) :: qaccel, qcorr, q_conj
    logical     :: ok_flag

    coeff0 = 0.0
    coeff1 = 0.5 * dt
    coeff2 = 1.0

    !$OMP PARALLEL DEFAULT (NONE)                                             &
    !$OMP PRIVATE (i, raccel, rcorr, rotmat, rotmat_trans, J, J_b_inv, J_inv, &
    !$OMP          q_conj, momentum, omtd1, qaccel, qcorr, ok_flag)           &
    !$OMP SHARED (coeff0, coeff1, coeff2, n_ice, ice, kin)

    !$OMP DO REDUCTION (+: kin)
    do i = 1, n_ice
      !=======================!
      ! translat acceleration !
      !=======================!
      raccel = ice(i)%force / ice(i)%mass
      rcorr  = raccel - ice(i)%rtd2

      !=========================!
      ! quaternion acceleration !
      !=========================!
      !calculating the rotation matrix and its transpose from the quaternion
      rotmat       = Q2matrix (ice(i)%q)
      rotmat_trans = transpose (rotmat)

      !inertia tensor and its inverse in space-fixed coordinate system
      J = matmul (rotmat, matmul (ice(i)%inertia, rotmat_trans))
      call M33inv (ice(i)%inertia, J_b_inv, ok_flag)
      J_inv = matmul (rotmat, matmul (J_b_inv, rotmat_trans))

      q_conj = .conj. ice(i)%q

      ice(i)%om = 2.0 * Q2vec (ice(i)%qtd1 * q_conj)

      !momentum as function of inertia and omega
      momentum = matmul (J, ice(i)%om)

      !first time derivative of omega, a.k.a. Euler's equation of rotative motion
      omtd1 = matmul (J_inv, Cross (momentum, ice(i)%om) + ice(i)%torque)

      qaccel = 0.5d0 * (Vec2q (omtd1) * ice(i)%q + ice(i)%qtd1 * Vec2q (ice(i)%om))
      qcorr  = qaccel - ice(i)%qtd2

      !============!
      ! Correction !
      !============!
      ice(i)%r    = ice(i)%r    + coeff0 * rcorr
      ice(i)%rtd1 = ice(i)%rtd1 + coeff1 * rcorr
      ice(i)%rtd2 = ice(i)%rtd2 + coeff2 * rcorr

      ice(i)%q    = ice(i)%q    + coeff0 * qcorr
      ice(i)%qtd1 = ice(i)%qtd1 + coeff1 * qcorr
      ice(i)%qtd2 = ice(i)%qtd2 + coeff2 * qcorr

      q_conj = .conj. ice(i)%q
      ice(i)%om = 2.0 * Q2vec (ice(i)%qtd1 * q_conj)

      ! kinetic energy with updated omega and quaternion
      kin = kin + ice(i)%inertia(1,1) * ice(i)%om(1) ** 2   &
        + ice(i)%inertia(2,2) * ice(i)%om(2) ** 2   &
        + ice(i)%inertia(3,3) * ice(i)%om(3) ** 2   &
        + ice(i)%mass         * ice(i)%rtd1(1) ** 2 &
        + ice(i)%mass         * ice(i)%rtd1(2) ** 2 &
        + ice(i)%mass         * ice(i)%rtd1(3) ** 2
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    kin = 0.5 * kin

    !force output in 50Hz
    if (step == 1) hertz = 1.0/50.0
    if (step*dt > hertz) then
      call OutputForces (dt, kin, n_struct, option)
      hertz = hertz + 1.0/50.0
    end if
  end subroutine Corrector

  subroutine OutputForces (dt, kin, n_struct, option)
    use mod_Global
    implicit none
    double precision :: dt, kin, time
    integer          :: i, n_struct, option

    time = step * dt

    if (option /= 1) then
      write (30, *) time, ';' , sum ( [(structure(i)%force(1), i = 1, n_struct)] )
      write (31, *) time, ';' , sum ( [(structure(i)%force(2), i = 1, n_struct)] )
      write (32, *) time, ';' , sum ( [(structure(i)%force(3), i = 1, n_struct)] )
    end if

    write (26, *) time, ';', norm2 (ice(i)%force)
    write (33, *) time, ';', kin
  end subroutine OutputForces
end module mod_PredCorr
