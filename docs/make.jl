using Documenter
using Skylight

makedocs(
    sitename = "Skylight",
    format = Documenter.HTML(),
    modules = [Skylight],
    checkdocs = :none, #Change this to :all in the future, once all docstrings are supposed to be included.
    pages = [
        "Home" => "index.md",
        "Getting started" => "gettingstarted.md",
        "Tutorials" => [
            "Novikov-Thorne disk around a Kerr black hole" => "examples/disk_kerr.md",
            "Lamppost corona above a Kerr black hole" => "examples/corona_kerr.md",
            "Ion torus around a Kerr black hole" => "examples/torus_kerr.md",
            "f(R)-Kerr black hole" => "examples/frkerr.md",
            "Johannsen black hole" => "examples/johannsen.md",
            "Boson star" => "examples/bosonstar.md",
            "Star across charged wormhole" => "examples/chargedwormhole.md",
            "RAR galactic core" => "examples/rar.md",
            "Skylight logo" => "examples/skylightlogo.md",
        ],
        "Spacetimes" => "spacetimes.md",
        "Radiative models" => "radiativemodels.md",
        "Camera" => "camera.md",
        "Configurations" => "configurations/configurations.md",
        "Initial data" => "initialdata.md",
        "Radiative transfer" => "radiativetransfer/radiativetransfer.md",
        "Postprocess" => [
            "Images" => "postprocess/images.md",
            "Spectra" => "postprocess/spectra.md",
            "Line emission" => "postprocess/lineemission.md",
            "Light curves" => "postprocess/lightcurves.md",
        ],
        "Utils" => [
            "Geometry and linear algebra" => "utils/geometry.md",
            "Vector fields" => "utils/vectorfields.md",
            "Constants" => "utils/constants.md",
            "Units and dimensions" => "utils/units.md",
        ],
        "Developers" => [
            "Adding a spacetime" => "developers/spacetimes.md",
            "Adding a radiative model" => "developers/radiativemodels.md"
        ],
        "api.md",
        "publications.md",
        "citing.md",
        "contributing.md",
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/joaquinpelle/Skylight.jl.git",
    branch = "gh-pages",
    target = "build",
)
