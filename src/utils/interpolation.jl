function build_interpolator(file; delim='\t', col_x=1, col_y=2, kind="linear")
    data = readdlm(file, delim, Float64, '\n')
    return my_interpolation(data[:,col_y], data[:,col_x], kind=kind)
end


function my_interpolation(ydata, xdata; kind=nothing)
    
    if kind == "linear"
        return LinearInterpolation(ydata,xdata)
    elseif kind == "cubicspline" 
        return CubicSpline(ydata,xdata)
    elseif kind == "loglinear"
        return LinearInterpolation(log10.(ydata),log10.(xdata))
    elseif kind == "logcubicspline" 
        return CubicSpline(log10.(ydata),log10.(xdata))
    else error("Kind not defined.")        
    end

end

function type_of_interpolator(kind)
    typeof(my_interpolation(ones(2),[0.0, 1.0]; kind=kind))
end