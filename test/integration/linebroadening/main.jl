using Skylight
using CairoMakie
using DelimitedFiles

set_theme!(; fonts = (; regular = "Times New Roman"))

include("linebroadening.jl")

line_broadening_test(ProgradeRotation(), 
    npixels = [2000],
    tols = [1e-15],
    num_bins = [30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200])

line_broadening_test(RetrogradeRotation(), 
    npixels = [2000],
    tols = [1e-15],
    num_bins = [30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200])