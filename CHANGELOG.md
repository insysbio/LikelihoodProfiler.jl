# Change Log

<<<<<<< HEAD
## 0.4.0 - support for CB

- changes in get_right_endpoint for :CICO_ONE_PASS to use function scan_loss_func(x::Vector{Float64})
- add get_endpoint for functions (works for CICO only)


## 0.3.2 - experimental

An experimental version includes a draft version of `get_endpoint()` for calculating confidence bands (CB).

## 0.3.1

- handling errors in loss function
- add :LOSS_ERROR_STOP
- fix bug in constraint for :CICO
- repository structure: add manifest file
- documentation updates

## 0.3.0
=======
## 0.3.2

- box constraints in CICO
- move Jupiter notebook to another repos insysbio/likelihoodprofiler-cases
- update docs and readme

## 0.3.1
>>>>>>> e5d7fed7750e5a0face3309d580270df3da06a1c

- add case: PK model with saturation
- handling errors in loss function
- add :LOSS_ERROR_STOP
- fix bug in constraint for :CICO
- repository structure: add manifest file
- add new tests
- documentation updates

## 0.2.1

- added tests for plot and adapted grid
- precompilation added
- multiple bug fixes

## 0.2.0

- code refactoring
- plot parameter profiles

## 0.1.0 - First Release

- :ONE_PATH method for estimation of confidence intervals. Based on CICO approach, see documentation
- :D2D_PLE method for estimation of confidence intervals. Uses rewritten derivation free method from d2d
- optional transformation to :log space for particular parameters for :ONE_PATH method
unit tests for test cases series
- params_plot function for preparation smooth curve based on modified adapted_grid from https://github.com/JuliaPlots/PlotUtils.jl
