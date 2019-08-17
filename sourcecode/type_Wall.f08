module type_Wall
  type, public :: wall
    logical          :: reverse
    double precision :: mu
    double precision :: mass
    double precision :: rtd1(3)
    double precision :: r(3)
    double precision :: om(3)
    double precision :: youngs
    double precision :: force(3)
    double precision :: torque(3)
    double precision :: length, width, thickness
    integer :: faces, vertices
    double precision :: face_equation(4, 6)
    integer :: face_vertex_table(4, 6)
    integer :: vertex_face_table(3, 6)
    double precision :: vert_coord_body(3, 8)
    double precision :: vert_coord_space(3, 8)
  end type
end module type_Wall
