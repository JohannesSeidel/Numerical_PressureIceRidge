module type_Iceblock
  use type_Quaternion
  type, public :: iceblock
    double precision :: width, length, thickness
    double precision :: density
    double precision :: mu
    integer          :: faces
    integer          :: vertices
    double precision :: mass
    double precision :: inertia(3,3)
    double precision :: youngs
    double precision :: vert_coord_space(3,8)
    double precision :: vert_coord_body(3,8)
    double precision :: face_equation(4,6)
    integer          :: face_vertex_table(4, 6)
    integer          :: vertex_face_table(3, 8)
    double precision :: r(3), rtd1(3), rtd2(3)
    double precision :: om(3)
    type (qtrn)      :: q, qtd1, qtd2, qconj !Quaternion
    double precision :: force(3), torque(3)
    double precision :: A !mean reference area for calculating the drag forces
  end type iceblock
end module type_Iceblock
