# Code for: *Quantum jumps are inescapable*

## Code structure
The code in repository is structured such that the bulk of the code lives inside `/src/main.jl` along useful utility code in `/src/utils.jl` both of which are cojointly called by the helper file `/src/src.jl` which simply invoques them. This is what the command `include("../src/src.jl")`is doing in all the rest of the repository. Finally, you'll find all the code for generating the figures inside `/data` and the resulting figures in `/figs`. 

We hope the code is well documented enough for you to follow along. Otherwise, dont hesitate to contribute or post questions in this repository.

## Enviroment
Here's some relevant info on the versions of `Julia` and the packages used:
- `julia_version = "1.10.1"`
- `DataFrames = "1.6.1"`
- `LaTeXStrings = "1.3.1"`
- `LinearAlgebra = see Julia version`
- `Plots = "1.40.4"`
- `ProgressBars = "1.5.1"`
- `StatsBase = "0.34.3"`

of course all the info is inside `/spin_idea_env/Manifest.toml`.

## Contributing 
If you spot any mistakes, you want to improve the code or simply want to discuss either submit an issue to GitHub or send us an email.


## Enjoy!

#
#

# To do:
- Actually document.
- Split data generation and data plotting