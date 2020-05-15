using Documenter, CamparyBindings

makedocs(;
    modules=[CamparyBindings],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/orkolorko/CamparyBindings.jl/blob/{commit}{path}#L{line}",
    sitename="CamparyBindings.jl",
    authors="Isaia Nisoli",
    assets=String[],
)

deploydocs(;
    repo="github.com/orkolorko/CamparyBindings.jl",
)
