module type_Element
  type, public :: element
    double precision, pointer :: force(:)
    double precision, pointer :: mass
    double precision, pointer :: mu
    double precision, pointer :: r(:)
    double precision, pointer :: rtd1(:)
    double precision, pointer :: om(:)
    double precision, pointer :: youngs
    double precision, pointer :: torque(:)
    double precision, pointer :: length, width, thickness
    integer,          pointer :: faces, vertices
    double precision, pointer :: face_equation(:, :)
    integer,          pointer :: face_vertex_table(:, :)
    integer,          pointer :: vertex_face_table(:, :)
    double precision, pointer :: vert_coord_space(:, :)
  end type
end module type_Element
