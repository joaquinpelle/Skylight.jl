@inline function det4x4(A)
    @inbounds return (A[13] * A[10] * A[7] * A[4] - A[9] * A[14] * A[7] * A[4] -
                      A[13] * A[6] * A[11] * A[4] + A[5] * A[14] * A[11] * A[4] +
                      A[9] * A[6] * A[15] * A[4] - A[5] * A[10] * A[15] * A[4] -
                      A[13] * A[10] * A[3] * A[8] + A[9] * A[14] * A[3] * A[8] +
                      A[13] * A[2] * A[11] * A[8] - A[1] * A[14] * A[11] * A[8] -
                      A[9] * A[2] * A[15] * A[8] + A[1] * A[10] * A[15] * A[8] +
                      A[13] * A[6] * A[3] * A[12] - A[5] * A[14] * A[3] * A[12] -
                      A[13] * A[2] * A[7] * A[12] + A[1] * A[14] * A[7] * A[12] +
                      A[5] * A[2] * A[15] * A[12] - A[1] * A[6] * A[15] * A[12] -
                      A[9] * A[6] * A[3] * A[16] + A[5] * A[10] * A[3] * A[16] +
                      A[9] * A[2] * A[7] * A[16] - A[1] * A[10] * A[7] * A[16] -
                      A[5] * A[2] * A[11] * A[16] + A[1] * A[6] * A[11] * A[16])
end
@inline function inv4x4(A)
    ideterminant = 1 / det4x4(A)
    @inbounds B = @SMatrix [(A[2, 3] * A[3, 4] * A[4, 2] - A[2, 4] * A[3, 3] * A[4, 2] +
                             A[2, 4] * A[3, 2] * A[4, 3] - A[2, 2] * A[3, 4] * A[4, 3] -
                             A[2, 3] * A[3, 2] * A[4, 4] + A[2, 2] * A[3, 3] * A[4, 4]) *
                            ideterminant
        (A[2, 4] * A[3, 3] * A[4, 1] - A[2, 3] * A[3, 4] * A[4, 1] -
         A[2, 4] * A[3, 1] * A[4, 3] + A[2, 1] * A[3, 4] * A[4, 3] +
         A[2, 3] * A[3, 1] * A[4, 4] - A[2, 1] * A[3, 3] * A[4, 4]) * ideterminant
        (A[2, 2] * A[3, 4] * A[4, 1] - A[2, 4] * A[3, 2] * A[4, 1] +
         A[2, 4] * A[3, 1] * A[4, 2] - A[2, 1] * A[3, 4] * A[4, 2] -
         A[2, 2] * A[3, 1] * A[4, 4] + A[2, 1] * A[3, 2] * A[4, 4]) * ideterminant
        (A[2, 3] * A[3, 2] * A[4, 1] - A[2, 2] * A[3, 3] * A[4, 1] -
         A[2, 3] * A[3, 1] * A[4, 2] + A[2, 1] * A[3, 3] * A[4, 2] +
         A[2, 2] * A[3, 1] * A[4, 3] - A[2, 1] * A[3, 2] * A[4, 3]) * ideterminant
        (A[1, 4] * A[3, 3] * A[4, 2] - A[1, 3] * A[3, 4] * A[4, 2] -
         A[1, 4] * A[3, 2] * A[4, 3] + A[1, 2] * A[3, 4] * A[4, 3] +
         A[1, 3] * A[3, 2] * A[4, 4] - A[1, 2] * A[3, 3] * A[4, 4]) * ideterminant
        (A[1, 3] * A[3, 4] * A[4, 1] - A[1, 4] * A[3, 3] * A[4, 1] +
         A[1, 4] * A[3, 1] * A[4, 3] - A[1, 1] * A[3, 4] * A[4, 3] -
         A[1, 3] * A[3, 1] * A[4, 4] + A[1, 1] * A[3, 3] * A[4, 4]) * ideterminant
        (A[1, 4] * A[3, 2] * A[4, 1] - A[1, 2] * A[3, 4] * A[4, 1] -
         A[1, 4] * A[3, 1] * A[4, 2] + A[1, 1] * A[3, 4] * A[4, 2] +
         A[1, 2] * A[3, 1] * A[4, 4] - A[1, 1] * A[3, 2] * A[4, 4]) * ideterminant
        (A[1, 2] * A[3, 3] * A[4, 1] - A[1, 3] * A[3, 2] * A[4, 1] +
         A[1, 3] * A[3, 1] * A[4, 2] - A[1, 1] * A[3, 3] * A[4, 2] -
         A[1, 2] * A[3, 1] * A[4, 3] + A[1, 1] * A[3, 2] * A[4, 3]) * ideterminant
        (A[1, 3] * A[2, 4] * A[4, 2] - A[1, 4] * A[2, 3] * A[4, 2] +
         A[1, 4] * A[2, 2] * A[4, 3] - A[1, 2] * A[2, 4] * A[4, 3] -
         A[1, 3] * A[2, 2] * A[4, 4] + A[1, 2] * A[2, 3] * A[4, 4]) * ideterminant
        (A[1, 4] * A[2, 3] * A[4, 1] - A[1, 3] * A[2, 4] * A[4, 1] -
         A[1, 4] * A[2, 1] * A[4, 3] + A[1, 1] * A[2, 4] * A[4, 3] +
         A[1, 3] * A[2, 1] * A[4, 4] - A[1, 1] * A[2, 3] * A[4, 4]) * ideterminant
        (A[1, 2] * A[2, 4] * A[4, 1] - A[1, 4] * A[2, 2] * A[4, 1] +
         A[1, 4] * A[2, 1] * A[4, 2] - A[1, 1] * A[2, 4] * A[4, 2] -
         A[1, 2] * A[2, 1] * A[4, 4] + A[1, 1] * A[2, 2] * A[4, 4]) * ideterminant
        (A[1, 3] * A[2, 2] * A[4, 1] - A[1, 2] * A[2, 3] * A[4, 1] -
         A[1, 3] * A[2, 1] * A[4, 2] + A[1, 1] * A[2, 3] * A[4, 2] +
         A[1, 2] * A[2, 1] * A[4, 3] - A[1, 1] * A[2, 2] * A[4, 3]) * ideterminant
        (A[1, 4] * A[2, 3] * A[3, 2] - A[1, 3] * A[2, 4] * A[3, 2] -
         A[1, 4] * A[2, 2] * A[3, 3] + A[1, 2] * A[2, 4] * A[3, 3] +
         A[1, 3] * A[2, 2] * A[3, 4] - A[1, 2] * A[2, 3] * A[3, 4]) * ideterminant
        (A[1, 3] * A[2, 4] * A[3, 1] - A[1, 4] * A[2, 3] * A[3, 1] +
         A[1, 4] * A[2, 1] * A[3, 3] - A[1, 1] * A[2, 4] * A[3, 3] -
         A[1, 3] * A[2, 1] * A[3, 4] + A[1, 1] * A[2, 3] * A[3, 4]) * ideterminant
        (A[1, 4] * A[2, 2] * A[3, 1] - A[1, 2] * A[2, 4] * A[3, 1] -
         A[1, 4] * A[2, 1] * A[3, 2] + A[1, 1] * A[2, 4] * A[3, 2] +
         A[1, 2] * A[2, 1] * A[3, 4] - A[1, 1] * A[2, 2] * A[3, 4]) * ideterminant
        (A[1, 2] * A[2, 3] * A[3, 1] - A[1, 3] * A[2, 2] * A[3, 1] +
         A[1, 3] * A[2, 1] * A[3, 2] - A[1, 1] * A[2, 3] * A[3, 2] -
         A[1, 2] * A[2, 1] * A[3, 3] + A[1, 1] * A[2, 2] * A[3, 3]) * ideterminant]
    SMatrix{4, 4, eltype(B), 16}(B)
end
"""
Computes the inverse of a 4x4 symmetric matrix using the closed-form solution for the inverse.

Parameters:
- inv_matrix: mutable array of size (4,4) to store the resulting inverse matrix.
- matrix: array of size (4,4) containing the input symmetric matrix.

Returns: nothing.
"""
function inv4x4sym!(inv_matrix::Matrix{Float64}, matrix::Matrix{Float64})
    @inbounds begin
        a, b, c, d = matrix[1, 1], matrix[1, 2], matrix[1, 3], matrix[1, 4]
        e, f, g = matrix[2, 2], matrix[2, 3], matrix[2, 4]
        h, i = matrix[3, 3], matrix[3, 4]
        j = matrix[4, 4]

        A = e * h * j - e * i^2 - f^2 * j + 2 * f * g * i - g^2 * h
        B = b * h * j - b * i^2 - f * c * j + f * i * d + g * c * i - g * h * d
        C = b * f * j - b * g * i - c * e * j + c * g^2 + d * e * i - g * d * f
        D = b * f * i - b * g * h - c * e * i + c * g * f + d * e * h - d * f^2
        det = a * A - b * B + c * C - d * D

        det != 0 || error("The matrix is singular and cannot be inverted.")

        inv_det = 1.0 / det

        inv_matrix[1, 1] = A * inv_det
        inv_matrix[1, 2] = inv_matrix[2, 1] = -B * inv_det
        inv_matrix[1, 3] = inv_matrix[3, 1] = C * inv_det
        inv_matrix[1, 4] = inv_matrix[4, 1] = -D * inv_det

        E = a * h * j - a * i^2 - c^2 * j + 2 * c * d * i - d^2 * h
        F = a * f * j - a * g * i - b * c * j + b * d * i + d * c * g - d^2 * f
        G = a * f * i - a * g * h - b * c * i + b * d * h + c^2 * g - c * d * f

        inv_matrix[2, 2] = E * inv_det
        inv_matrix[2, 3] = inv_matrix[3, 2] = -F * inv_det
        inv_matrix[2, 4] = inv_matrix[4, 2] = G * inv_det

        H = a * e * j - a * g^2 - b^2 * j + 2 * b * g * d - d^2 * e
        I = a * e * i - a * g * f - b^2 * i + +b * f * d + b * c * g - c * e * d

        inv_matrix[3, 3] = H * inv_det
        inv_matrix[3, 4] = inv_matrix[4, 3] = -I * inv_det

        J = a * e * h - a * f^2 - b^2 * h + 2 * b * c * f - c^2 * e

        inv_matrix[4, 4] = J * inv_det
    end

    return nothing
end

"""
Computes the determinant of a 4x4 symmetric matrix using a closed-form formula.

Parameters:
- matrix: array of size (4,4) containing the input symmetric matrix.

Returns: the determinant of the input matrix.
"""
function det4x4sym(matrix)
    @inbounds begin
        a, b, c, d = matrix[1, 1], matrix[1, 2], matrix[1, 3], matrix[1, 4]
        e, f, g = matrix[2, 2], matrix[2, 3], matrix[2, 4]
        h, i = matrix[3, 3], matrix[3, 4]
        j = matrix[4, 4]
    end

    A = e * h * j - e * i^2 - f^2 * j + 2 * f * g * i - g^2 * h
    B = b * h * j - b * i^2 - f * c * j + f * i * d + g * c * i - g * h * d
    C = b * f * j - b * g * i - c * e * j + c * g^2 + d * e * i - g * d * f
    D = b * f * i - b * g * h - c * e * i + c * g * f + d * e * h - d * f^2
    det = a * A - b * B + c * C - d * D

    return det
end
