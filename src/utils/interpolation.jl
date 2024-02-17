function build_interpolator(file; delim = '\t', col_x = 1, col_y = 2, kind = "linear", kwargs...)
    data = readdlm(file, delim, Float64, '\n')
    return my_interpolation(data[:, col_y], data[:, col_x], kind = kind, kwargs...)
end

function my_interpolation(ydata, xdata; kind, kwargs...)
    if kind == "linear"
        return LinearInterpolation(ydata, xdata; kwargs...)
    elseif kind == "cubicspline"
        return CubicSpline(ydata, xdata; kwargs...)
    elseif kind == "loglinear"
        return LinearInterpolation(log10.(ydata), log10.(xdata); kwargs...)
    elseif kind == "logcubicspline"
        return CubicSpline(log10.(ydata), log10.(xdata); kwargs...)
    else
        error("Kind $kind not defined.")
    end
end

function type_of_interpolator(kind)
    typeof(my_interpolation(ones(2), [0.0, 1.0]; kind = kind))
end
