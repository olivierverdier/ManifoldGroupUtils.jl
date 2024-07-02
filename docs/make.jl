push!(LOAD_PATH, "../src/")

using Documenter
using GroupTools


makedocs(
    sitename = "GroupTools",
    format = Documenter.HTML(),
    modules = [GroupTools]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
