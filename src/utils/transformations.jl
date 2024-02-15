transform_vector!(vy, dy_dx, vx) = mul!(vy, dy_dx, vx)

to_static(u::AbstractVector) = SVector{length(u), eltype(u)}(u)