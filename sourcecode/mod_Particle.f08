module mod_Particle
  !
  ! This module contains the initialization routines for the particle outlines, masses and moments of intertia. The update routines
  ! should also be kept here. For polygons, this would mean the positions of the corners; for ellipses, it would mean the
  ! coefficients of the corresponding polynomials. Additionals information such as the shapes of the particles at, e.g., the initial
  ! orientation can also be useful to keep here, as well as data structures for the sides. The bounding boxes should also be part of
  ! and computed in this module.
  !
  ! STARTED: 23.2.'15
  !

  ! the order of element types in the elements-vector is: 1) ice rubble 2) structure 3) consolidated layer and basin walls
contains
  subroutine InitializeParticles (density, length, mu, n_ice, n_struct, n_wall, thickness, width, youngs)
    use mod_Global
    use mod_Constants
    use mod_Functions
    use mod_Neighbour

    implicit none
    integer,          intent (in) :: n_ice, n_wall, n_struct
    double precision, intent (in) :: length, width, thickness, density, youngs, mu
    double precision :: diagonal
    double precision :: nx, ny, nz
    double precision :: v1(3), v2(3), v3(3)
    integer          :: i, j
    integer          :: n_elem
    double precision :: theta, phi, psi
    double precision :: varL, varT, varW !variation in length, thickness and width
    double precision, allocatable :: rand_num(:,:)

    allocate (rand_num(3,n_ice))

    varL = 0.5  !50%
    varT = 0.05 !5%
    varW = 0.1  !10%

    n_elem = n_ice + n_wall + n_struct 

    !============!
    ! ICE RUBBLE !
    !============!
    allocate (ice(n_ice))

    call Init_random_seed ()
    call random_number (rand_num)

    do i = 1, n_ice
      ice(i)%faces     = 6
      ice(i)%vertices  = 8
      ice(i)%youngs    = youngs
      ice(i)%density   = density ! kg/m^3
      ice(i)%mu        = mu

      ice(i)%length    = (1.0 - varL) * length    + 2.0 * rand_num(1,i) * varL * length
      ice(i)%width     = (1.0 - varW) * width     + 2.0 * rand_num(2,i) * varW * width
      ice(i)%thickness = (1.0 - varT) * thickness + 2.0 * rand_num(3,i) * varT * thickness

      ! mean reference area = sum of all face areas devided by number of faces
      ice(i)%A = (ice(i)%length*ice(i)%width + ice(i)%length*ice(i)%thickness + ice(i)%width*ice(i)%thickness) / 3.0

      ice(i)%mass          = ice(i)%length * ice(i)%width * ice(i)%thickness * ice(i)%density
      ice(i)%inertia       = 0.0
      ice(i)%inertia(1,1) = ice(i)%mass / 12.0 * (ice(i)%width**2  + ice(i)%thickness**2)
      ice(i)%inertia(2,2) = ice(i)%mass / 12.0 * (ice(i)%length**2 + ice(i)%thickness**2)
      ice(i)%inertia(3,3) = ice(i)%mass / 12.0 * (ice(i)%length**2 + ice(i)%width**2)
      ice(i)%vert_coord_body = reshape (&
        [ -0.5*ice(i)%width, -0.5*ice(i)%length,  0.5*ice(i)%thickness, &
        & -0.5*ice(i)%width,  0.5*ice(i)%length,  0.5*ice(i)%thickness, &
        &  0.5*ice(i)%width,  0.5*ice(i)%length,  0.5*ice(i)%thickness, &
        &  0.5*ice(i)%width, -0.5*ice(i)%length,  0.5*ice(i)%thickness, &
        & -0.5*ice(i)%width, -0.5*ice(i)%length, -0.5*ice(i)%thickness, &
        & -0.5*ice(i)%width,  0.5*ice(i)%length, -0.5*ice(i)%thickness, &
        &  0.5*ice(i)%width,  0.5*ice(i)%length, -0.5*ice(i)%thickness, &
        &  0.5*ice(i)%width, -0.5*ice(i)%length, -0.5*ice(i)%thickness], &
        shape (ice(i)%vert_coord_body) )
      ice(i)%face_vertex_table = reshape ([1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8, 1,2,3,4, 5,6,7,8],  &
        shape (ice(i)%face_vertex_table))
      ice(i)%vertex_face_table = reshape ([1,4,5, 1,2,5, 2,3,5, 3,4,5, 1,4,6, 1,2,6, 2,3,6, 3,4,6],&
        shape (ice(i)%vertex_face_table))
    end do

    diagonal = maxval ( [(sqrt (ice(i)%length**2 + ice(i)%width**2 + ice(i)%thickness**2), i = 1, n_ice)] ) 

    nx = 0
    ny = 1
    nz = 2

    do i = 1, n_ice
      !-------------------!
      ! initial positions !
      !-------------------!
      if (((nx + 1) * diagonal) <= rdge%width) then 
        nx = nx + 1
      else !end of ridge is reached and a new row is added
        nx = 1
        ny = ny + 1
        if (ny * diagonal >= rdge%length) then !a new layer z-1 is added
          ny = 1
          nz = nz + 1
        end if
      end if

      ice(i)%r(1) = (nx - 0.5) * diagonal
      ice(i)%r(2) = (ny - 0.5) * diagonal
      ice(i)%r(3) = (nz - 0.5) * diagonal

      ice(i)%rtd1(1) = (rand_num(1,i)-0.5)
      ice(i)%rtd1(2) = (rand_num(2,i)-0.5)
      ice(i)%rtd1(3) = (rand_num(3,i)-0.5)

      ice(i)%rtd2 = 0.0

      ice(i)%om   = 0.0
      ice(i)%qtd1%w = 0.0; ice(i)%qtd1%x = 0.0; ice(i)%qtd1%y = 0.0; ice(i)%qtd1%z = 0.0
      ice(i)%qtd2%w = 0.0; ice(i)%qtd2%x = 0.0; ice(i)%qtd2%y = 0.0; ice(i)%qtd2%z = 0.0

      theta = 2 * PI * rand_num(1,i)
      phi   = 2 * PI * rand_num(2,i)
      psi   = 2 * PI * rand_num(3,i)

      ice(i)%q%w = cos (0.5 * theta) * cos (0.5 * (phi + psi) )
      ice(i)%q%x = sin (0.5 * theta) * cos (0.5 * (phi - psi) )
      ice(i)%q%y = sin (0.5 * theta) * sin (0.5 * (phi - psi) )
      ice(i)%q%z = cos (0.5 * theta) * sin (0.5 * (phi + psi) )

      !Update vertex coordinates (the function 'Rotate_q' is declared in mod_Quaternion)
      do j = 1, ice(i)%vertices
        ice(i)%vert_coord_space(:, j) = ice(i)%r + Rotate_q (ice(i)%vert_coord_body(:, j), ice(i)%q)
      end do
    end do

    !============!
    ! BOUNDARIES !
    !============!

    allocate (boundaries(n_wall))

    !--------------------!
    ! CONSOLITATED LAYER !
    !--------------------!

    boundaries(1)%thickness = 0.05
    boundaries(1)%mass      = 10.0d+09
    boundaries(1)%mu        = mu
    boundaries(1)%r         = [0.5 * rdge%width, 0.5 * rdge%length, 0.5 * 0.9 * boundaries(1)%thickness]
    boundaries(1)%rtd1      = 0.0
    boundaries(1)%om        = 0.0
    boundaries(1)%youngs    = youngs
    boundaries(1)%force     = 0.0
    boundaries(1)%torque    = 0.0
    boundaries(1)%faces     = 6
    boundaries(1)%vertices  = 8
    boundaries(1)%vert_coord_space = reshape (&
      [ -1.0d+00,              0.0d+00,      0.045d+00 , &
      & -1.0d+00,              rdge%length,  0.045d+00 , &
      &  rdge%width + 3.0d+00, rdge%length,  0.045d+00 , &
      &  rdge%width + 3.0d+00, 0.0d+00,      0.045d+00 , &
      & -1.0d+00,              0.0d+00,     -0.005d+00 , &
      & -1.0d+00,              rdge%length, -0.005d+00 , &
      &  rdge%width + 3.0d+00, rdge%length, -0.005d+00 , &
      &  rdge%width + 3.0d+00, 0.0d+00,     -0.005d+00], &
      shape (boundaries(1)%vert_coord_space))

    boundaries(1)%face_vertex_table = reshape ([1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8, 1,2,3,4, 5,6,7,8],  &
      shape (boundaries(1)%face_vertex_table))
    boundaries(1)%vertex_face_table = reshape ([1,4,5, 1,2,5, 2,3,5, 3,4,5, 1,4,6, 1,2,6, 2,3,6, 3,4,6],&
      shape (boundaries(1)%vertex_face_table))

    ! face equations calculated with vertex coordinates
    do i = 1, boundaries(1)%faces
      v1 = boundaries(1)%vert_coord_space(:, boundaries(1)%face_vertex_table(1, i))
      v2 = boundaries(1)%vert_coord_space(:, boundaries(1)%face_vertex_table(2, i))
      v3 = boundaries(1)%vert_coord_space(:, boundaries(1)%face_vertex_table(3, i))

      boundaries(1)%face_equation(1:3, i) = Cross ( (v1 - v2), (v1 - v3) ) / norm2 (Cross ( (v1 - v2), (v1 - v3) ))
      if ( dot_product( boundaries(1)%face_equation(1:3, i), ( boundaries(1)%r - v1 ) ) > 0.0 ) then
        boundaries(1)%face_equation(1:3, i) = -1 * boundaries(1)%face_equation(1:3, i)
      end if

      boundaries(1)%face_equation(4, i) = dot_product( boundaries(1)%face_equation(1:3, i), v1)
    end do

    !-------------!
    ! basin walls !
    !-------------!
    boundaries(2)%width     = 0.5
    boundaries(2)%thickness = 0.05
    boundaries(2)%mass      = 10.0d+09
    boundaries(2)%mu        = mu
    boundaries(2)%r         = [0.5 * rdge%length, -0.5 * boundaries(2)%width, 1.25d+00]
    boundaries(2)%rtd1      = 0.0
    boundaries(2)%om        = 0.0
    boundaries(2)%youngs    = youngs
    boundaries(2)%force     = 0.0
    boundaries(2)%torque    = 0.0
    boundaries(2)%faces     = 6
    boundaries(2)%vertices  = 8
    boundaries(2)%vert_coord_space = reshape (&
      [ -1.0d+00,              -boundaries(2)%width,  10.0d+00, &
      & -1.0d+00,              -0.001d+00,            10.0d+00, &
      &  rdge%width + 3.0d+00, -0.001d+00,            10.0d+00, &
      &  rdge%width + 3.0d+00, -boundaries(2)%width,  10.0d+00, &
      & -1.0d+00,              -boundaries(2)%width, -0.5d+00,  &
      & -1.0d+00,              -0.001d+00,           -0.5d+00,  &
      &  rdge%width + 3.0d+00, -0.001d+00,           -0.5d+00,  &
      &  rdge%width + 3.0d+00, -boundaries(2)%width, -0.5d+00], &
      shape (boundaries(2)%vert_coord_space))

    boundaries(2)%face_vertex_table = reshape ([1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8, 1,2,3,4, 5,6,7,8],  &
      shape (boundaries(2)%face_vertex_table))
    boundaries(2)%vertex_face_table = reshape ([1,4,5, 1,2,5, 2,3,5, 3,4,5, 1,4,6, 1,2,6, 2,3,6, 3,4,6],&
      shape (boundaries(2)%vertex_face_table))

    ! face equations calculated with vertex coordinates
    do i = 1, boundaries(2)%faces
      v1 = boundaries(2)%vert_coord_space(:, boundaries(2)%face_vertex_table(1, i))
      v2 = boundaries(2)%vert_coord_space(:, boundaries(2)%face_vertex_table(2, i))
      v3 = boundaries(2)%vert_coord_space(:, boundaries(2)%face_vertex_table(3, i))

      boundaries(2)%face_equation(1:3, i) = Cross ( (v1 - v2), (v1 - v3) ) / norm2 (Cross ( (v1 - v2), (v1 - v3) ))
      if ( dot_product( boundaries(2)%face_equation(1:3, i), ( boundaries(2)%r - v1 ) ) > 0.0 ) then
        boundaries(2)%face_equation(1:3, i) = -1 * boundaries(2)%face_equation(1:3, i)
      end if

      boundaries(2)%face_equation(4, i) = dot_product( boundaries(2)%face_equation(1:3, i), v1)
    end do

    boundaries(3)%width     = 0.5
    boundaries(3)%thickness = 0.05
    boundaries(3)%mass      = 10.0d+09
    boundaries(3)%mu        = mu
    boundaries(3)%r         = [0.5 * rdge%width, rdge%length + 0.5 * boundaries(2)%width, 1.25d+00]
    boundaries(3)%rtd1      = 0.0
    boundaries(3)%om        = 0.0
    boundaries(3)%youngs    = youngs
    boundaries(3)%force     = 0.0
    boundaries(3)%torque    = 0.0
    boundaries(3)%faces     = 6
    boundaries(3)%vertices  = 8
    boundaries(3)%vert_coord_space = reshape (&
      [ -1.0d+00             , rdge%length + 0.001d+00          ,  10.0d+00, &
      & -1.0d+00             , rdge%length + boundaries(3)%width,  10.0d+00, &
      &  rdge%width + 1.0d+00, rdge%length + boundaries(3)%width,  10.0d+00, &
      &  rdge%width + 1.0d+00, rdge%length + 0.001d+00          ,  10.0d+00, &
      & -1.0d+00             , rdge%length + 0.001d+00          , -0.5d+00,  &
      & -1.0d+00             , rdge%length + boundaries(3)%width, -0.5d+00,  &
      &  rdge%width + 1.0d+00, rdge%length + boundaries(3)%width, -0.5d+00,  &
      &  rdge%width + 1.0d+00, rdge%length + 0.001d+00          , -0.5d+00], &
      shape (boundaries(3)%vert_coord_space))

    boundaries(3)%face_vertex_table = reshape ([1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8, 1,2,3,4, 5,6,7,8],  &
      shape (boundaries(3)%face_vertex_table))
    boundaries(3)%vertex_face_table = reshape ([1,4,5, 1,2,5, 2,3,5, 3,4,5, 1,4,6, 1,2,6, 2,3,6, 3,4,6],&
      shape (boundaries(3)%vertex_face_table))

    ! face equations calculated with vertex coordinates
    do i = 1, boundaries(3)%faces
      v1 = boundaries(3)%vert_coord_space(:, boundaries(3)%face_vertex_table(1, i))
      v2 = boundaries(3)%vert_coord_space(:, boundaries(3)%face_vertex_table(2, i))
      v3 = boundaries(3)%vert_coord_space(:, boundaries(3)%face_vertex_table(3, i))

      boundaries(3)%face_equation(1:3, i) = Cross ( (v1 - v2), (v1 - v3) ) / norm2 (Cross ( (v1 - v2), (v1 - v3) ))
      if ( dot_product( boundaries(3)%face_equation(1:3, i), ( boundaries(3)%r - v1 ) ) > 0.0 ) then
        boundaries(3)%face_equation(1:3, i) = -1 * boundaries(3)%face_equation(1:3, i)
      end if

      boundaries(3)%face_equation(4, i) = dot_product( boundaries(3)%face_equation(1:3, i), v1)
    end do

    allocate (structure(n_struct))

    !===================================!
    ! INITIALIZE ELEMENTS POINTER ARRAY !
    !===================================!
    allocate (elements(n_elem))

    !============!
    ! RUBBLE ICE !
    !============!
    do i = 1, n_ice
      elements(i)%mass      => ice(i)%mass
      elements(i)%mu        => ice(i)%mu
      elements(i)%rtd1      => ice(i)%rtd1
      elements(i)%r         => ice(i)%r
      elements(i)%om        => ice(i)%om
      elements(i)%youngs    => ice(i)%youngs
      elements(i)%force     => ice(i)%force
      elements(i)%torque    => ice(i)%torque
      elements(i)%faces     => ice(i)%faces
      elements(i)%vertices  => ice(i)%vertices
      elements(i)%face_equation     => ice(i)%face_equation
      elements(i)%face_vertex_table => ice(i)%face_vertex_table
      elements(i)%vert_coord_space  => ice(i)%vert_coord_space
      elements(i)%vertex_face_table => ice(i)%vertex_face_table
    end do

    !=======!
    ! WALLS !
    !=======!
    do i = 1, n_wall
      elements(n_ice + i)%mass     => boundaries(i)%mass 
      elements(n_ice + i)%mu       => boundaries(i)%mu 
      elements(n_ice + i)%rtd1     => boundaries(i)%rtd1
      elements(n_ice + i)%r        => boundaries(i)%r
      elements(n_ice + i)%om       => boundaries(i)%om
      elements(n_ice + i)%youngs   => boundaries(i)%youngs
      elements(n_ice + i)%force    => boundaries(i)%force
      elements(n_ice + i)%torque   => boundaries(i)%torque
      elements(n_ice + i)%faces    => boundaries(i)%faces
      elements(n_ice + i)%vertices => boundaries(i)%vertices
      elements(n_ice + i)%face_equation     => boundaries(i)%face_equation
      elements(n_ice + i)%face_vertex_table => boundaries(i)%face_vertex_table
      elements(n_ice + i)%vert_coord_space  => boundaries(i)%vert_coord_space
      elements(n_ice + i)%vertex_face_table => boundaries(i)%vertex_face_table
    end do

    !===========!
    ! STRUCTURE !
    !===========!
    do i = 1, n_struct
      elements(n_ice + n_wall + i)%mass     => structure(i)%mass 
      elements(n_ice + n_wall + i)%mu       => structure(i)%mu
      elements(n_ice + n_wall + i)%rtd1     => structure(i)%rtd1
      elements(n_ice + n_wall + i)%r        => structure(i)%r
      elements(n_ice + n_wall + i)%om       => structure(i)%om
      elements(n_ice + n_wall + i)%youngs   => structure(i)%youngs
      elements(n_ice + n_wall + i)%force    => structure(i)%force
      elements(n_ice + n_wall + i)%torque   => structure(i)%torque
      elements(n_ice + n_wall + i)%faces    => structure(i)%faces
      elements(n_ice + n_wall + i)%vertices => structure(i)%vertices
      elements(n_ice + n_wall + i)%face_equation     => structure(i)%face_equation
      elements(n_ice + n_wall + i)%face_vertex_table => structure(i)%face_vertex_table
      elements(n_ice + n_wall + i)%vert_coord_space  => structure(i)%vert_coord_space
      elements(n_ice + n_wall + i)%vertex_face_table => structure(i)%vertex_face_table
    end do

    !=====================!
    ! INIT BOUNDING BOXES !
    !=====================!

    !Values are assigned in subroutine UpdateParticles
    allocate (boundingBoxes(n_Elem))
    allocate (xEndPoints(2*n_elem)); allocate (yEndPoints(2*n_elem)); allocate (zEndPoints(2*n_Elem))
    allocate (overlapPairs(10*n_elem + n_struct*n_elem + n_wall*n_elem))

    !initialize just the array including ice and boundaries. The structure, if there is one, is added in "MakeStructure"
    j = 1
    do i = 1, 2*(n_ice+n_wall), 2
      xEndPoints(i)%mValue   =>  boundingBoxes(j)%mMin(1)
      xEndPoints(i)%mData    =  -j
      xEndPoints(i+1)%mValue =>  boundingBoxes(j)%mMax(1)
      xEndPoints(i+1)%mData  =   j

      yEndPoints(i)%mValue   =>  boundingBoxes(j)%mMin(2)
      yEndPoints(i)%mData    =  -j
      yEndPoints(i+1)%mValue =>  boundingBoxes(j)%mMax(2)
      yEndPoints(i+1)%mData  =   j

      zEndPoints(i)%mValue   =>  boundingBoxes(j)%mMin(3)
      zEndPoints(i)%mData    =  -j
      zEndPoints(i+1)%mValue =>  boundingBoxes(j)%mMax(3)
      zEndPoints(i+1)%mData  =   j

      j = j + 1
    end do
  end subroutine InitializeParticles

  subroutine MakeStructure (n_ice, n_struct, n_wall, option)
    use mod_Constants
    use mod_Global
    use mod_Functions
    implicit none
    integer, intent (in) :: n_ice, n_struct, n_wall, option
    integer          :: i, id, j
    double precision :: v1(3), v2(3), v3(3)
    double precision :: angle, axis(3)
    type (qtrn)      :: q

    select case (option)
    case (1) !ridge creation bars
      angle = atan ((0.5d+00 * (rdge%width - rdge%keel_width)) / rdge%keel_height)
      axis  = [0.0d+00, 1.0d+00, 0.0d+00]

      do i = 1, n_struct
        if (i == 1) then
          structure(i)%r(1) = minval ( [(ice(j)%vert_coord_space(1,:), j = 1, n_ice)] ) - 0.5d+00
          structure(i)%rtd1 = [0.2d+00, 0.0d+00, 0.0d+00]
        else if (i == 2) then
          structure(i)%r(1) = maxval ( [(ice(j)%vert_coord_space(1,:), j = 1, n_ice)] ) + 0.5d+00
          structure(i)%rtd1 = [-0.2d+00, 0.0d+00, 0.0d+00]
        end if
        structure(i)%r(2) = 0.5d+00 * rdge%length
        structure(i)%r(3) = maxval ( [(ice(j)%vert_coord_space(3,:), j = 1, n_ice)] ) + 0.5d+00

        structure(i)%reverse   = .false. !reverse the velocity
        structure(i)%mass      = 10.0d+09
        structure(i)%mu        = 0.0d+00
        structure(i)%om        = 0.0d+00
        structure(i)%youngs    = 200.0d+09
        structure(i)%force     = 0.0d+00
        structure(i)%torque    = 0.0d+00
        structure(i)%faces     = 6
        structure(i)%vertices  = 8
        structure(i)%length    = 0.1d+00
        structure(i)%width     = rdge%length
        structure(i)%thickness = 2.5d0 * structure(i)%r(3)
        structure(i)%vert_coord_body = reshape (&
          [ -0.5d0 * structure(i)%length, -0.5d0 * structure(i)%width,  0.01d0 * structure(i)%thickness,  &
          & -0.5d0 * structure(i)%length,  0.5d0 * structure(i)%width,  0.01d0 * structure(i)%thickness,  &
          &  0.5d0 * structure(i)%length,  0.5d0 * structure(i)%width,  0.01d0 * structure(i)%thickness,  &
          &  0.5d0 * structure(i)%length, -0.5d0 * structure(i)%width,  0.01d0 * structure(i)%thickness,  &
          & -0.5d0 * structure(i)%length, -0.5d0 * structure(i)%width, -0.99d0 * structure(i)%thickness,  &
          & -0.5d0 * structure(i)%length,  0.5d0 * structure(i)%width, -0.99d0 * structure(i)%thickness,  &
          &  0.5d0 * structure(i)%length,  0.5d0 * structure(i)%width, -0.99d0 * structure(i)%thickness,  &
          &  0.5d0 * structure(i)%length, -0.5d0 * structure(i)%width, -0.99d0 * structure(i)%thickness], &
          shape (structure(i)%vert_coord_space))

        if (i == 2) axis = -1.0d+00 * axis

        q%w = cos (0.5d0 * angle)
        q%x = axis(1) * sin (0.5d0 * angle)
        q%y = axis(2) * sin (0.5d0 * angle)
        q%z = axis(3) * sin (0.5d0 * angle)

        do j = 1, structure(i)%vertices
          structure(i)%vert_coord_space(:,j) = structure(i)%r + Rotate_q (structure(i)%vert_coord_body(:,j), q)
        end do

        structure(i)%face_vertex_table = reshape ([1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8, 1,2,3,4, 5,6,7,8],  &
          shape (structure(i)%face_vertex_table))
        structure(i)%vertex_face_table = reshape ([1,4,5, 1,2,5, 2,3,5, 3,4,5, 1,4,6, 1,2,6, 2,3,6, 3,4,6],&
          shape (structure(i)%vertex_face_table))

        ! face equations calculated with vertex coordinates
        do j = 1, structure(i)%faces
          v1 = structure(i)%vert_coord_space(:, structure(i)%face_vertex_table(1, j))
          v2 = structure(i)%vert_coord_space(:, structure(i)%face_vertex_table(2, j))
          v3 = structure(i)%vert_coord_space(:, structure(i)%face_vertex_table(3, j))

          structure(i)%face_equation(1:3, j) = Cross ( (v2 - v1), (v3 - v1) ) / norm2 (Cross ( (v2 - v1), (v3 - v1) ))
          if ( dot_product( structure(i)%face_equation(1:3, j), ( structure(i)%r - v1 ) ) > 0.0d0 ) then
            structure(i)%face_equation(1:3, j) = -1 * structure(i)%face_equation(1:3, i)
          end if

          structure(i)%face_equation(4, j) = dot_product( structure(i)%face_equation(1:3, j), v1)
        end do
      end do
    case (2) !punch device
      angle = 0.25d0 * PI
      axis  = [0.0d0, 0.0d0, 1.0d0]

      q%w = cos (0.5d0 * angle)
      q%x = axis(1) * sin (0.5d0 * angle)
      q%y = axis(2) * sin (0.5d0 * angle)
      q%z = axis(3) * sin (0.5d0 * angle)

      structure(1)%length    = sqrt (PI * 0.125 ** 2) !diameter of punch device is 0.25m
      structure(1)%width     = sqrt (PI * 0.125 ** 2)
      structure(1)%thickness = 1.0d0
      structure(1)%mass      = 122 !kg
      structure(1)%mu        = 0.1d0
      structure(1)%r         = [0.5d+00*rdge%width, 0.5d0*rdge%length, -0.5d+00*structure(1)%thickness]
      structure(1)%rtd1      = [0.00d+00, 0.0d+00, 0.0067d+00]
      structure(1)%om        = 0.0
      structure(1)%youngs    = 200.0d+09 !steel
      structure(1)%force     = 0.0
      structure(1)%torque    = 0.0
      structure(1)%faces     = 6
      structure(1)%vertices  = 8
      structure(1)%vert_coord_body = reshape (&
        [ -0.5d0 * structure(1)%length, -0.5d0 * structure(1)%width,  0.5d0 * structure(1)%thickness,  &
        & -0.5d0 * structure(1)%length,  0.5d0 * structure(1)%width,  0.5d0 * structure(1)%thickness,  &
        &  0.5d0 * structure(1)%length,  0.5d0 * structure(1)%width,  0.5d0 * structure(1)%thickness,  &
        &  0.5d0 * structure(1)%length, -0.5d0 * structure(1)%width,  0.5d0 * structure(1)%thickness,  &
        & -0.5d0 * structure(1)%length, -0.5d0 * structure(1)%width, -0.5d0 * structure(1)%thickness,  &
        & -0.5d0 * structure(1)%length,  0.5d0 * structure(1)%width, -0.5d0 * structure(1)%thickness,  &
        &  0.5d0 * structure(1)%length,  0.5d0 * structure(1)%width, -0.5d0 * structure(1)%thickness,  &
        &  0.5d0 * structure(1)%length, -0.5d0 * structure(1)%width, -0.5d0 * structure(1)%thickness], &
        shape (structure(1)%vert_coord_space))

      do j = 1, structure(1)%vertices
        structure(1)%vert_coord_space(:,j) = structure(1)%r + Rotate_q (structure(1)%vert_coord_body(:,j), q)
      end do

      structure(1)%face_vertex_table = reshape ([1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8, 1,2,3,4, 5,6,7,8],  &
        shape (structure(1)%face_vertex_table))
      structure(1)%vertex_face_table = reshape ([1,4,5, 1,2,5, 2,3,5, 3,4,5, 1,4,6, 1,2,6, 2,3,6, 3,4,6],&
        shape (structure(1)%vertex_face_table))

      ! face equations calculated with vertex coordinates
      do j = 1, structure(1)%faces
        v1 = structure(1)%vert_coord_space(:, structure(1)%face_vertex_table(1, j))
        v2 = structure(1)%vert_coord_space(:, structure(1)%face_vertex_table(2, j))
        v3 = structure(1)%vert_coord_space(:, structure(1)%face_vertex_table(3, j))

        structure(1)%face_equation(1:3, j) = Cross ( (v1 - v2), (v1 - v3) ) / norm2 (Cross ( (v1 - v2), (v1 - v3) ))
        if ( dot_product( structure(1)%face_equation(1:3, j), ( structure(1)%r - v1 ) ) > 0.0d0 ) then
          structure(1)%face_equation(1:3, j) = -1 * structure(1)%face_equation(1:3, j)
        end if

        structure(1)%face_equation(4, j) = dot_product( structure(1)%face_equation(1:3, j), v1)
      end do
    case (3) !test model
      do i = 1, n_struct
        angle = 0.25 * PI
        angle = 0.00
        axis  = [0.0, 0.0, 1.0]

        q%w = cos (0.5 * angle)
        q%x = axis(1) * sin (0.5 * angle)
        q%y = axis(2) * sin (0.5 * angle)
        q%z = axis(3) * sin (0.5 * angle)

        structure(i)%length    = 0.7
        structure(i)%width     = 0.7
        structure(i)%thickness = 0.7
        structure(i)%mass      = 200
        structure(i)%mu        = 0.11
        structure(i)%r         = [-structure(i)%length, 0.5*rdge%length, 0.22 + 0.5*structure(i)%thickness]
        structure(i)%rtd1      = [0.045d+00, 0.0d+00, 0.0d+00]
        structure(i)%om        = 0.0
        structure(i)%youngs    = 200.0d+09
        structure(i)%force     = 0.0
        structure(i)%torque    = 0.0
        structure(i)%faces     = 6
        structure(i)%vertices  = 8
        structure(i)%vert_coord_body = reshape (&
          [ -0.5d0 * structure(i)%length, -0.5d0 * structure(i)%width,  0.5d0 * structure(i)%thickness,  &
          & -0.5d0 * structure(i)%length,  0.5d0 * structure(i)%width,  0.5d0 * structure(i)%thickness,  &
          &  0.5d0 * structure(i)%length,  0.5d0 * structure(i)%width,  0.5d0 * structure(i)%thickness,  &
          &  0.5d0 * structure(i)%length, -0.5d0 * structure(i)%width,  0.5d0 * structure(i)%thickness,  &
          & -0.5d0 * structure(i)%length, -0.5d0 * structure(i)%width, -0.5d0 * structure(i)%thickness,  &
          & -0.5d0 * structure(i)%length,  0.5d0 * structure(i)%width, -0.5d0 * structure(i)%thickness,  &
          &  0.5d0 * structure(i)%length,  0.5d0 * structure(i)%width, -0.5d0 * structure(i)%thickness,  &
          &  0.5d0 * structure(i)%length, -0.5d0 * structure(i)%width, -0.5d0 * structure(i)%thickness], &
          shape (structure(i)%vert_coord_space))

        do j = 1, structure(i)%vertices
          structure(i)%vert_coord_space(:,j) = structure(i)%r + Rotate_q (structure(i)%vert_coord_body(:,j), q)
        end do

        structure(i)%face_vertex_table = reshape ([1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8, 1,2,3,4, 5,6,7,8],  &
          shape (structure(i)%face_vertex_table))
        structure(i)%vertex_face_table = reshape ([1,4,5, 1,2,5, 2,3,5, 3,4,5, 1,4,6, 1,2,6, 2,3,6, 3,4,6],&
          shape (structure(i)%vertex_face_table))

        ! face equations calculated with vertex coordinates
        do j = 1, structure(i)%faces
          v1 = structure(i)%vert_coord_space(:, structure(i)%face_vertex_table(1, j))
          v2 = structure(i)%vert_coord_space(:, structure(i)%face_vertex_table(2, j))
          v3 = structure(i)%vert_coord_space(:, structure(i)%face_vertex_table(3, j))

          structure(i)%face_equation(1:3, j) = Cross ( (v1 - v2), (v1 - v3) ) / norm2 (Cross ( (v1 - v2), (v1 - v3) ))
          if ( dot_product( structure(i)%face_equation(1:3, j), ( structure(i)%r - v1 ) ) > 0.0d0 ) then
            structure(i)%face_equation(1:3, j) = -1 * structure(i)%face_equation(1:3, i)
          end if

          structure(i)%face_equation(4, j) = dot_product( structure(i)%face_equation(1:3, j), v1)
        end do
      end do
    end select

    !add structure to existing x-,y- and zEndPoints
    do i = 1, n_struct
      id = n_ice + n_wall + i
      xEndPoints( 2*id-1 )%mData = -id
      yEndPoints( 2*id-1 )%mData = -id
      zEndPoints( 2*id-1 )%mData = -id

      xEndPoints( 2*id-1 )%mValue => boundingBoxes(id)%mMin(1)
      yEndPoints( 2*id-1 )%mValue => boundingBoxes(id)%mMin(2)
      zEndPoints( 2*id-1 )%mValue => boundingBoxes(id)%mMin(3)

      xEndPoints( 2*id )%mData = id
      yEndPoints( 2*id )%mData = id
      zEndPoints( 2*id )%mData = id

      xEndPoints( 2*id )%mValue => boundingBoxes(id)%mMax(1)
      yEndPoints( 2*id )%mValue => boundingBoxes(id)%mMax(2)
      zEndPoints( 2*id )%mValue => boundingBoxes(id)%mMax(3)
    end do
  end subroutine MakeStructure

  subroutine UpdateParticles (dt, n_ice, n_elem, n_struct, option)
    use omp_lib
    use mod_Global
    use mod_Functions
    implicit none

    integer, intent (in) :: n_ice, n_elem, n_struct
    integer, intent (in) :: option
    double precision, intent (in) :: dt
    integer :: i, j
    double precision :: v1(3), v2(3), v3(3)
    double precision :: va(3), vb(3), vc(3), vd(3), n(3), p1(3), p2(3)
    double precision :: d

    !$OMP PARALLEL DEFAULT (NONE)    &
    !$OMP PRIVATE (i, j, v1, v2, v3) &
    !$OMP SHARED (ice, n_ice)

    !$OMP DO
    !==== UPDATE ICE ====!
    do i = 1, n_ice
      ! normalization of new quaternion, .norm. is declared in type_quaternion.f08 which itself is included in mod_Global
      call Norm_quaternion (ice(i)%q)

      !Update vertex coordinates (the function 'Rotate_q' is declared in mod_Quaternion)
      do j = 1, ice(i)%vertices
        ice(i)%vert_coord_space(:, j) = ice(i)%r + Rotate_q (ice(i)%vert_coord_body(:, j), ice(i)%q)
      end do

      !update face equations
      do j = 1, ice(i)%faces
        v1 = ice(i)%vert_coord_space(:, ice(i)%face_vertex_table(1, j))
        v2 = ice(i)%vert_coord_space(:, ice(i)%face_vertex_table(2, j))
        v3 = ice(i)%vert_coord_space(:, ice(i)%face_vertex_table(3, j))

        ice(i)%face_equation(1:3, j) = Cross ( (v1 - v2), (v1 - v3) ) / norm2 (Cross ( (v1 - v2), (v1 - v3) ))
        if ( dot_product( ice(i)%face_equation(1:3, j), ( ice(i)%r - v1 ) ) > 0.0d0 ) then
          ice(i)%face_equation(1:3, j) = -1 * ice(i)%face_equation(1:3, j)
        end if

        ice(i)%face_equation(4, j) = dot_product( ice(i)%face_equation(1:3, j), v1)
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    !==== UPDATE STRUCTURE ====!
    if (option == 1) then
      !plane equation of the water surface
      n = [0.0d+00, 0.0d+00, 1.0d+00]
      d = 0.0d+00
      !2 vertices each for the inner face of the slide bars, which touch the ice
      !elements
      va = structure(1)%vert_coord_space(:,4)
      vb = structure(1)%vert_coord_space(:,8)
      vc = structure(2)%vert_coord_space(:,1)
      vd = structure(2)%vert_coord_space(:,5)
      !the intersection points of the two edges of the faces and the water plane
      p1 = va + (d - dot_product (n, va) / dot_product (n, (vb - va))) * (vb - va)
      p2 = vc + (d - dot_product (n, vc) / dot_product (n, (vd - vc))) * (vd - vc)
      !reverse the velocity if the bars have reached the ridge's width
      if ((p2(1) - p1(1) <= rdge%width) .and. (structure(1)%reverse .eqv. .false.)) then
        do i = 1, n_struct
          structure(i)%reverse = .true.
          structure(i)%rtd1    = -1.0 * structure(i)%rtd1
        end do
      end if

      do i = 1, n_struct
        structure(i)%r = structure(i)%r + structure(i)%rtd1 * dt

        do j = 1, structure(i)%vertices
          structure(i)%vert_coord_space(:,j) = structure(i)%vert_coord_space(:,j) + structure(i)%rtd1 * dt
        end do

        do j = 1, structure(i)%faces
          v1 = structure(i)%vert_coord_space(:, structure(i)%face_vertex_table(1, j))
          v2 = structure(i)%vert_coord_space(:, structure(i)%face_vertex_table(2, j))
          v3 = structure(i)%vert_coord_space(:, structure(i)%face_vertex_table(3, j))

          structure(i)%face_equation(1:3, j) = Cross ( (v1 - v2), (v1 - v3) ) / norm2 (Cross ( (v1 - v2), (v1 - v3) ))
          if ( dot_product( structure(i)%face_equation(1:3, j), ( structure(i)%r - v1 ) ) > 0.0d0 ) then
            structure(i)%face_equation(1:3, j) = -1 * structure(i)%face_equation(1:3, j)
          end if

          structure(i)%face_equation(4, j) = dot_product( structure(i)%face_equation(1:3, j), v1)
        end do
      end do

    else if (option == 2) then
      structure(1)%r = structure(1)%r + structure(1)%rtd1 * dt

      do j = 1, structure(1)%vertices
        structure(1)%vert_coord_space(:,j) = structure(1)%vert_coord_space(:,j) + structure(1)%rtd1 * dt
      end do

      do j = 1, structure(1)%faces
        v1 = structure(1)%vert_coord_space(:, structure(1)%face_vertex_table(1, j))
        v2 = structure(1)%vert_coord_space(:, structure(1)%face_vertex_table(2, j))
        v3 = structure(1)%vert_coord_space(:, structure(1)%face_vertex_table(3, j))

        structure(1)%face_equation(1:3, j) = Cross ( (v1 - v2), (v1 - v3) ) / norm2 (Cross ( (v1 - v2), (v1 - v3) ))
        if ( dot_product( structure(1)%face_equation(1:3, j), ( structure(1)%r - v1 ) ) > 0.0 ) then
          structure(1)%face_equation(1:3, j) = -1 * structure(1)%face_equation(1:3, j)
        end if

        structure(1)%face_equation(4, j) = dot_product( structure(1)%face_equation(1:3, j), v1)
      end do

    else if (option == 3) then
      do i = 1, n_struct
        structure(i)%r = structure(i)%r + structure(i)%rtd1 * dt

        do j = 1, structure(i)%vertices
          structure(i)%vert_coord_space(:,j) = structure(i)%vert_coord_space(:,j) + structure(i)%rtd1 * dt
        end do

        do j = 1, structure(i)%faces
          v1 = structure(i)%vert_coord_space(:, structure(i)%face_vertex_table(1, j))
          v2 = structure(i)%vert_coord_space(:, structure(i)%face_vertex_table(2, j))
          v3 = structure(i)%vert_coord_space(:, structure(i)%face_vertex_table(3, j))

          structure(i)%face_equation(1:3, j) = Cross ( (v1 - v2), (v1 - v3) ) / norm2 (Cross ( (v1 - v2), (v1 - v3) ))
          if ( dot_product( structure(i)%face_equation(1:3, j), ( structure(i)%r - v1 ) ) > 0.0d0 ) then
            structure(i)%face_equation(1:3, j) = -1 * structure(i)%face_equation(1:3, j)
          end if

          structure(i)%face_equation(4, j) = dot_product( structure(i)%face_equation(1:3, j), v1)
        end do
      end do
    end if

    !==== UPDATE BOUNDINGBOXES ====!
    call UpdateBoundingBoxes (n_elem)
  end subroutine UpdateParticles

  subroutine UpdateBoundingBoxes (n_elem)
    use mod_Global
    implicit none
    integer, intent (in) :: n_elem
    integer :: i
    do i = 1, n_elem
      boundingBoxes(i)%mMax(1) = maxval (elements(i)%vert_coord_space(1,:))
      boundingBoxes(i)%mMax(2) = maxval (elements(i)%vert_coord_space(2,:))
      boundingBoxes(i)%mMax(3) = maxval (elements(i)%vert_coord_space(3,:))
      boundingBoxes(i)%mMin(1) = minval (elements(i)%vert_coord_space(1,:))
      boundingBoxes(i)%mMin(2) = minval (elements(i)%vert_coord_space(2,:))
      boundingBoxes(i)%mMin(3) = minval (elements(i)%vert_coord_space(3,:))
    end do
  end subroutine

  subroutine TerminationCheck (dt, n_elem, option, termination)
    ! Termination checks
    ! Ridge creation: After the push bar's velocity is reversed, the main loop
    !                 is terminated after 20 seconds
    ! Punch test: The main loop is terminated after the the bottom of the punch
    !             device penetrated the ridge plus additional 20cm  
    ! Model test: Main loop is terminated after the x-min coordinate of the model
    !             is greater than the ridge width plus additional 50cm  
    !
    use mod_Global
    implicit none
    double precision, intent (in)  :: dt
    integer,          intent (in)  :: option, n_elem
    logical,          intent (out) :: termination
    double precision, save         :: timer

    if (step == 1) timer = 0.0

    select case (option)
    case (1)
      if (structure(1)%reverse) then
        timer = timer + dt
        if (timer > 20.0) then
          write (*,*) "Ridge creation complete"
          termination = .true.
        end if
      end if
    case (2)
      if (boundingBoxes(n_elem)%mMax(3) > rdge%keel_height + 0.2) then
        write (*,*) "Punch test complete"
        termination = .true.
      end if
    case (3)
      !For structures containing more than 1 element this condition has to be
      !adjusted
      if (boundingBoxes(n_elem)%mMin(1) > rdge%width + 0.5) then
        write (*,*) "Model test complete"
        termination = .true.
      end if
    end select
  end subroutine
end module mod_Particle
