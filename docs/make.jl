using Documenter
using Skylight

makedocs(
    sitename = "Skylight",
    format = Documenter.HTML(),
    modules = [Skylight]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
