"""
Computes the inverse of a 4x4 symmetric matrix using the closed-form solution for the inverse.

Parameters:
- inv_matrix: mutable array of size (4,4) to store the resulting inverse matrix.
- matrix: array of size (4,4) containing the input symmetric matrix.

Returns: nothing.
"""
function inverse_4x4_symmetric!(inv_matrix::Matrix{Float64},matrix::Matrix{Float64})
    
    @inbounds begin
        a, b, c, d = matrix[1, 1], matrix[1, 2], matrix[1, 3], matrix[1, 4]
        e, f, g = matrix[2, 2], matrix[2, 3], matrix[2, 4]
        h, i = matrix[3, 3], matrix[3, 4]
        j = matrix[4, 4]

        A = e * h * j - e * i^2 - f^2 * j + 2 * f * g * i - g^2 * h
        B = b * h * j - b * i^2 - f * c * j + f * i * d + g * c * i - g * h * d
        C = b * f * j - b * g * i - c * e * j + c * g^2 + d * e * i - g * d * f
        D = b * f * i - b * g * h - c * e * i + c * g * f  + d * e * h - d * f^2
        det = a * A - b * B + c * C - d * D

        @assert det != 0 "The matrix is singular and cannot be inverted."

        inv_det = 1.0 / det

        inv_matrix[1, 1] = A * inv_det
        inv_matrix[1, 2] = inv_matrix[2, 1] = -B * inv_det
        inv_matrix[1, 3] = inv_matrix[3, 1] =  C * inv_det
        inv_matrix[1, 4] = inv_matrix[4, 1] = -D * inv_det

        E = a * h * j - a * i^2 - c^2 * j + 2 * c * d * i - d^2 * h
        F = a * f * j - a * g * i - b * c * j + b * d * i + d * c * g - d^2 * f
        G = a * f * i - a * g * h - b * c * i + b * d * h + c^2* g - c * d * f

        inv_matrix[2, 2] = E * inv_det
        inv_matrix[2, 3] = inv_matrix[3, 2] = -F * inv_det
        inv_matrix[2, 4] = inv_matrix[4, 2] = G * inv_det

        H = a * e * j - a * g^2 - b^2 * j + 2 * b * g * d - d^2 * e
        I = a * e * i - a * g * f - b^2 * i + + b * f * d + b * c * g - c * e * d

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
function determinant_4x4_symmetric(matrix)

    a, b, c, d = matrix[1, 1], matrix[1, 2], matrix[1, 3], matrix[1, 4]
    e, f, g = matrix[2, 2], matrix[2, 3], matrix[2, 4]
    h, i = matrix[3, 3], matrix[3, 4]
    j = matrix[4, 4]

    A = e * h * j - e * i^2 - f^2 * j + 2 * f * g * i - g^2 * h
    B = b * h * j - b * i^2 - f * c * j + f * i * d + g * c * i - g * h * d
    C = b * f * j - b * g * i - c * e * j + c * g^2 + d * e * i - g * d * f
    D = b * f * i - b * g * h - c * e * i + c * g * f  + d * e * h - d * f^2
    det = a * A - b * B + c * C - d * D

    return det

end