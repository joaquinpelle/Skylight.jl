using Documenter
using Skylight

makedocs(
    sitename = "Skylight",
    format = Documenter.HTML(),
    modules = [Skylight],
    pages = [
        "Home" => "index.md",
        "Getting started" => "gettingstarted.md",
        "Spacetimes" => [
            "Catalogue" => "spacetimes/catalogue.md",
            "Geometric functions" => "spacetimes/functions.md",
            "Abstract types and traits" => "spacetimes/othertypes.md"
            ],
        "Radiative models" => [
            "Catalogue" => "radiativemodels/catalogue.md",
            "Radiative processes" => "radiativemodels/radiativeprocesses.md",
            "Functions" => "radiativemodels/functions.md",
        ],
        "Configurations" => "configurations/configurations.md",
        "Initial data" => [
            "Pinhole camera" => "initialization/pinholecamera.md",
            "Initialization" => "initialization/initialization.md"],
        "Radiative transfer" => "radiativetransfer/radiativetransfer.md",
        "Postprocess" => [
            "Images" => "postprocess/images.md",
            "Spectra" => "postprocess/spectra.md",
            "Line emission" => "postprocess/lineemission.md",
            "Light curves" => "postprocess/lightcurves.md",
        ],
        "Examples" => [
            "Thin disk around a Kerr black hole" => "examples/disk_kerr.md",
            "Lamppost corona above a Kerr black hole" => "examples/corona_kerr.md",
            "Ion torus around a Kerr black hole" => "examples/torus_kerr.md",
            "f(R)-Kerr black hole" => "examples/frkerr.md",
            "Johannsen black hole" => "examples/johannsen.md",
            "Boson star" => "examples/bosonstar.md",
            "Star across charged wormhole" => "examples/chargedwormhole.md",
            "RAR galactic core" => "examples/rar.md",
            "Skylight's logo" => "examples/skylightlogo.md",
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
