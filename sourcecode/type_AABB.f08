module type_AABB
  type aabb
    double precision :: mMin(3), mMax(3)
  end type

  type endPoint
    integer                   :: mData
    double precision, pointer :: mValue
  end type
end module type_AABB
