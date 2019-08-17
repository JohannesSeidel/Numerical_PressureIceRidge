module type_Quaternion
  !
  ! A module for basic quaternion operations
  ! and 3D spatial rotations using quaternion representation
  ! Quaternions are a 4-component analogue to complex numbers
  !
  ! .. math::
  !
  ! \mathbf{q} = w + x \mathbf{i} + y \mathbf{j} + z \mathbf{k} = [w,x,y,z]
  !
  ! where the three imaginary components obey :math:`\mathbf{ijk} = \mathbf{i}^2 = \mathbf{j}^2 = \mathbf{k}^2 = -1`.
  ! This leads to a structure similar to that of complex numbers, except that the quaternion product defined
  ! according to the above rules is non-commutative :math:`\mathbf{q}_1\mathbf{q}_2 \ne \mathbf{q}_2\mathbf{q}_1`.
  !
  ! It turns out that unit quaternions :math:`||\mathbf{q}|| = \sqrt{w^2+x^2+y^2+z^2} = 1` represent the space
  ! of 3D rotations so that the rotation by angle :math:`\alpha` around unit axis :math:`\mathbf{u} = [x,y,z]`
  ! is represented by the quaternion
  !
  ! .. math::
  !
  ! \mathbf{q} = [\cos \frac{\alpha}{2}, x \sin \frac{\alpha}{2}, y \sin \frac{\alpha}{2}, z \sin \frac{\alpha}{2}]
  !
  ! and joining rotations is represented as a a quaternion product
  ! (rotating first by :math:`\mathbf{q}_1`, then by :math:`\mathbf{q}_2` yields the combined rotation
  ! of :math:`\mathbf{q}_{12} = \mathbf{q}_2 \mathbf{q}_1`).
  implicit none
  ! *norm_tolerance the threshold value for the norm for treating vectors as zero vectors
  double precision, parameter :: norm_tolerance = 1.0e-8

  ! The quarternion type. It only contains four
  ! double precision components, but the main advantage for
  ! defining it as a custom type is the possibility
  ! to write routines and operators for quaternion
  ! algebra.
  ! *w the "double precision" component of the quaternion
  ! *x an "imaginary" component of the quaternion
  ! *y an "imaginary" component of the quaternion
  ! *z an "imaginary" component of the quaternion
  type qtrn
    double precision :: w
    double precision :: x
    double precision :: y
    double precision :: z
  end type qtrn

  ! Quarternion product operator
  interface operator(.qp.)
    module procedure qprod
  end interface

  ! Quarternion product operator
  interface operator(*)
    module procedure qprod
  end interface

  ! Quarternion norm operator
  interface operator(.norm.)
    module procedure qnorm
  end interface

  ! Quarternion conjugate operator
  interface operator(.c.)
    module procedure qconj
  end interface

  ! Quarternion conjugate operator
  interface operator(.conj.)
    module procedure qconj
  end interface

  ! Quarternion inverse operator
  interface operator(.inv.)
    module procedure qinv
  end interface

  ! Quarternion to angle operator
  interface operator(.angle.)
    module procedure q2angle
  end interface

  ! Quarternion to axis operator
  interface operator(.axis.)
    module procedure q2axis
  end interface

  ! Quarternion to matrix operator
  interface operator(.mat.)
    module procedure q2matrix
  end interface

  ! Quarternion + scalar operator
  interface operator(+)
    module procedure qplus, qplusq
  end interface

  ! Quarternion - scalar operator
  interface operator(-)
    module procedure qminus, qminusq
  end interface

  ! Scalar * quarternion operator
  interface operator(*)
    module procedure qtimes
  end interface

  ! Quarternion / scalar operator
  interface operator(/)
    module procedure qdiv
  end interface

  ! Vector dot product
  interface operator(.o.)
    module procedure dot
  end interface

  ! Vector norm
  interface operator(.norm.)
    module procedure vnorm
  end interface

  ! Overloading the "rotate" subroutine
  interface rotate
    module procedure rotate_q, rotate_au, rotate_a
  end interface
contains
  function qplus (q, r) result (qn)
    ! Returns the quarternion :math:`\mathbf{q}` added by scalar :math:`r`
    ! component-wise
    ! *q a quaternion
    ! *r a double precision scalar
    ! *qn q+r
    implicit none
    double precision, intent (in) :: r
    type (qtrn),      intent (in) :: q
    type (qtrn)                   :: qn
    qn%w = q%w + r
    qn%x = q%x + r
    qn%y = q%y + r
    qn%z = q%z + r
  end function qplus

  function qplusq (q1, q2) result (qpq)
    ! Returns the sum of two quaternions
    ! 
    implicit none
    type (qtrn), intent (in) :: q1, q2
    type (qtrn)              :: qpq
    qpq%w = q1%w + q2%w
    qpq%x = q1%x + q2%x
    qpq%y = q1%y + q2%y
    qpq%z = q1%z + q2%z
  end function qplusq

  function qminus (q, r) result (qn)
    ! Returns the quarternion :math:`\mathbf{q}` subtracted by scalar :math:`r`
    ! component-wise
    ! *q a quaternion
    ! *r a double precision scalar
    ! *qn q-r
    implicit none
    double precision,        intent (in) :: r
    type (qtrn), intent (in) :: q
    type (qtrn)              :: qn
    qn%w = q%w - r
    qn%x = q%x - r
    qn%y = q%y - r
    qn%z = q%z - r
  end function qminus

  function qminusq (q1, q2) result (qmq)
    ! Returns the quaternion subtraction
    !
    implicit none
    type (qtrn), intent (in) :: q1, q2
    type (qtrn)              :: qmq
    qmq%w = q1%w - q2%w
    qmq%x = q1%x - q2%x
    qmq%y = q1%y - q2%y
    qmq%z = q1%z - q2%z
  end function qminusq

  function qtimes (r, q) result (qn)
    ! Returns the quaternion :math:`\mathbf{q}` multiplied by scalar :math:`r`
    ! component-wise
    ! *q a quaternion
    ! *r a double precision scalar
    ! *qn r*q
    implicit none
    double precision,        intent (in) :: r
    type (qtrn), intent (in) :: q
    type (qtrn)              :: qn
    qn%w = q%w * r
    qn%x = q%x * r
    qn%y = q%y * r
    qn%z = q%z * r
  end function qtimes

  function qdiv(q,r) result(qn)
    ! Returns the quarternion :math:`\mathbf{q}` divided by scalar :math:`r`
    ! component-wise
    ! *q a quaternion
    ! *r a double precision scalar
    ! *qn q/r
    implicit none
    double precision,       intent (in) :: r
    type(qtrn), intent (in) :: q
    type(qtrn)              :: qn
    qn%w = q%w / r
    qn%x = q%x / r
    qn%y = q%y / r
    qn%z = q%z / r 
  end function qdiv

  function qprod (q1,q2) result (qn)
    ! Returns the quarternion product :math:`\mathbf{q}_1\mathbf{q}_2`
    ! Note that the product is non-commutative: :math:`\mathbf{q}_1\mathbf{q}_2 \ne \mathbf{q}_2\mathbf{q}_1`
    ! *q1 a quaternion
    ! *q2 a quaternion
    ! *qn q1*q2
    implicit none
    type (qtrn), intent (in) :: q1, q2
    type (qtrn)              :: qn
    qn%w = q1%w*q2%w - q1%x*q2%x - q1%y*q2%y - q1%z*q2%z
    qn%x = q1%x*q2%w + q1%w*q2%x - q1%z*q2%y + q1%y*q2%z
    qn%y = q1%y*q2%w + q1%z*q2%x + q1%w*q2%y - q1%x*q2%z
    qn%z = q1%z*q2%w - q1%y*q2%x + q1%x*q2%y + q1%w*q2%z
  end function qprod

  function qconj (q) result (cq)
    ! Returns the quarternion conjugate of :math:`\mathbf{q}`: :math:`\mathbf{q}^* = w+x\mathbf{i}+y\mathbf{j}+z\mathbf{k} \to w-x\mathbf{i}-y\mathbf{j}-z\mathbf{k}`
    ! *q a quaternion
    ! *cq conjugate of q
    implicit none
    type(qtrn), intent (in) :: q
    type(qtrn)              :: cq
    cq%w =  q%w
    cq%x = -q%x
    cq%y = -q%y
    cq%z = -q%z
  end function qconj

  function qnorm (q) result (norm)
    ! Returns the quarternion norm of :math:`\mathbf{q}`: :math:`||\mathbf{q}|| = \sqrt{w^2+x^2+y^2+z^2}`
    ! *q a quaternion
    ! *norm the norm of q
    implicit none
    type (qtrn), intent (in) :: q
    double precision         :: norm
    norm = sqrt (q%w*q%w + q%x*q%x + q%y*q%y + q%z*q%z)
  end function qnorm

  function qinv (q) result (iq)
    ! Returns the quarternion inverse of :math:`\mathbf{q}`: :math:`\mathbf{q}^*/||\mathbf{q}||`
    ! *q a quaternion
    ! *iq inverse of q
    implicit none
    type(qtrn), intent (in) :: q
    type(qtrn)              :: iq
    iq = qconj (q) / qnorm (q)
  end function qinv

  function rot2q(a,u) result(q)
    ! Returns the quarternion representing a rotation
    ! around axis :math:`\mathbf{u}` by angle :math:`\alpha`
    ! *a angle in radians
    ! *u 3D vector, defining an axis of rotation
    ! *q representing the rotation
    implicit none
    double precision, intent (in) :: a, u(3)
    type (qtrn)                   :: q
    double precision              :: normer
    q%w    = cos (a / 2)
    normer = sin (a / 2) / sqrt (u(1)**2 + u(2)**2 + u(3)**2)
    q%x    = u(1) * normer
    q%y    = u(2) * normer
    q%z    = u(3) * normer
  end function rot2q

  function Vec2q (v) result (q)
    ! Returns the quarternion representing a rotation
    ! around axis :math:`\mathbf{v}` by angle :math:`||\mathbf{v}||`. If :math:`\mathbf{v} = 0`,
    ! the quaternion :math:`\mathbf{q} = [1 0 0 0]` will be returned.
    ! *v 3D vector, defining both the angle and axis of rotation
    ! *q representing the rotation
    implicit none
    double precision, intent (in) :: v(3)
    type (qtrn)                   :: q
    q%w = 0.0d0
    q%x = v(1)
    q%y = v(2)
    q%z = v(3)
    !double precision              :: normer, a
    !a = .norm.v
    !if (abs(a) < norm_tolerance) then ! zero rotation
    !  q%w = 1.0
    !  q%x = 0.0
    !  q%y = 0.0
    !  q%z = 0.0
    !  return
    !end if
    !q%w    = cos (a / 2)
    !normer = sin (a / 2) / a
    !q%x    = v(1) * normer
    !q%y    = v(2) * normer
    !q%z    = v(3) * normer
  end function vec2q
  
  function Q2Vec (q) result (v)
    ! Returns the vector part of a quaternion
    implicit none
    type (qtrn) :: q
    double precision :: v(3)
    v(1) = q%x
    v(2) = q%y
    v(3) = q%z
  end function Q2Vec

  function q2angle(q) result(a)
    ! Returns the angle of rotation described by the
    ! UNIT quarternion :math:`\mathbf{q}`. Note that the unity of :math:`\mathbf{q}`
    ! is not checked (as it would be time consuming to
    ! calculate the norm all the time if we know
    ! the quaternions used have unit length).
    ! *q a quaternion representation of rotation
    ! *a angle of rotation
    implicit none
    double precision         :: a
    type (qtrn), intent (in) :: q
    a = 2.0 * acos (q%w)
  end function q2angle

  function Euler2q (theta, phi, psi) result (q)
    !Returns the quaternion for given euler angles
    !ref: Understanding the Discrete Element Method, p.27
    implicit none
    double precision, intent (in) :: theta, phi, psi
    type (qtrn)                   :: q
    q%w = cos (theta / 2.0) * cos ((phi + psi) / 2.0)
    q%x = sin (theta / 2.0) * cos ((phi - psi) / 2.0)
    q%y = sin (theta / 2.0) * sin ((phi - psi) / 2.0)
    q%z = cos (theta / 2.0) * sin ((phi + psi) / 2.0)
  end function Euler2q

  function q2axis(q) result(u)
    ! Returns the axis of rotation described by the
    ! UNIT quarternion :math:`\mathbf{q}`. Note that the unity of :math:`\mathbf{q}`
    ! is not checked (as it would be time consuming to
    ! calculate the norm all the time if we know
    ! the quaternions used have unit length).
    ! *q a quaternion representation of rotation
    ! *u axis of rotation
    implicit none
    double precision        :: u(3), normer
    type(qtrn), intent (in) :: q
    normer = 1.0 / sin(q2angle(q))
    u(1)   = q%x * normer
    u(2)   = q%y * normer
    u(3)   = q%z * normer
  end function q2axis

  function Q2matrix (q) result (mat)
    ! Returns the rotation matrix described by the
    ! UNIT quarternion :math:`\mathbf{q}`. Note that the unity of :math:`\mathbf{q}`
    ! is not checked (as it would be time consuming to
    ! calculate the norm all the time if we know
    ! the quaternions used have unit length).
    ! *q a quaternion representation of rotation
    ! *mat rotation matrix
    !ref: Understanding the Dicscrete Element Method, p.27
    ! anm.: just works with the convention of a quaternion used in UtDEM
    implicit none
    type (qtrn), intent (in) :: q
    double precision         :: mat(3,3) 
    double precision         :: w, x, y, z
    w = q%w
    x = q%x
    y = q%y
    z = q%z
    mat(1, 1) = w ** 2 + x ** 2 - y ** 2 + z ** 2
    mat(1, 2) = 2.0d0 * (x * y + w * z)
    mat(1, 3) = 2.0d0 * (x * z - w * y)
    mat(2, 1) = 2.0d0 * (x * y - w * z)
    mat(2, 2) = w ** 2 - x ** 2 + y ** 2 - z ** 2
    mat(2, 3) = 2.0d0 * (y * z + w * x)
    mat(3, 1) = 2.0d0 * (x * z + w * y)
    mat(3, 2) = 2.0d0 * (y * z - w * x)
    mat(3, 3) = w ** 2 - x ** 2 - y ** 2 + z ** 2
  end function Q2matrix

  function Rotate_q (vec,q) result (v)
    ! Returns the 3D vector rotated according to
    ! the UNIT quarternion :math:`\mathbf{q}`. Note that the unity of :math:`\mathbf{q}`
    ! is not checked (as it would be time consuming to
    ! calculate the norm all the time if we know
    ! the quaternions used have unit length).
    ! *q a quaternion representation of rotation
    ! *vec vector to be rotated
    ! *v the rotated vector
    implicit none
    double precision, intent (in) :: vec(3)
    type(qtrn),       intent (in) :: q
    double precision              :: v(3)
    type(qtrn)                    :: temp, newq
    temp%w = 0.0d0
    temp%x = vec(1)
    temp%y = vec(2)
    temp%z = vec(3)
    newq = q .qp. temp .qp. (.c.q)
    v(1) = newq%x
    v(2) = newq%y
    v(3) = newq%z
  end function rotate_q

  function rotate_au(vec,a,u) result(v)
    ! Returns the vector rotated according to
    ! the axis :math:`\mathbf{u}` and angle :math:`\alpha`.
    ! *u axis of rotation
    ! *a angle of rotation
    ! *vec vector to be rotated
    ! *v the rotated vector
    implicit none
    double precision, intent (in) :: vec(3),a,u(3)
    double precision              :: v(3)
    v = rotate_q(vec,rot2q(a,u))
  end function rotate_au

  function rotate_a(vec,da) result(v)
    ! Returns the vector rotated according to
    ! the vector :math:`\mathbf{d}`. The axis of rotation is given by
    ! the direction of :math:`\mathbf{d}` and the angle by :math:`||\mathbf{d}||`.
    ! *da rotation vector (e.g., angular velocity x time :math:`\mathbf{\omega} t`)
    ! *vec vector to be rotated
    ! *v the rotated vector
    implicit none
    double precision, intent (in) :: vec(3),da(3)
    double precision :: v(3)
    v = rotate_q(vec,rot2q(sqrt(da(1)**2 + da(2)**2 + da(3)**2),da) )
  end function rotate_a

  function dot(v,u) result(product)
    ! Normal dot product of vectors :math:`\mathbf{v}\cdot\mathbf{u}` (Note: for 3-vectors only!)
    ! *v vector
    ! *u vector
    ! *product v . u
    implicit none
    double precision, intent (in) :: v(3), u(3)
    double precision              :: product
    product = v(1)*u(1) + v(2)*u(2) + v(3)*u(3)
  end function dot

  function vnorm(v) result(normi)
    ! Norm of a vector, :math:`||\mathbf{v}||`
    ! *v vector
    ! *normi norm of v
    implicit none
    double precision, intent (in) :: v(3)
    double precision              :: normi
    normi = sqrt(v.o.v)
  end function vnorm

  subroutine norm_quaternion(qq)
    ! norms the given quaternion
    ! *qq quaternion to be normed to unity
    implicit none
    type(qtrn) :: qq
    qq = qq / (.norm. qq)
    !if (abs(.norm.qq) > norm_tolerance) then
    !  qq = qq/.norm.qq
    !else
    !  qq%w = 1.0
    !  qq%x = 0.0
    !  qq%y = 0.0
    !  qq%z = 0.0
    !end if
    return
  end subroutine norm_quaternion
end module type_Quaternion
