# Change Log

## 0.5.1

- add method `get_optimal`
- update gradient-based interface, `loss_grad`, `scan_grad` options
- support of `loss_grad`, `scan_grad` as explicit functions

## 0.5.0

- support of scales in `profile`
- use `:lin` as a synonym of `:direct`
- julia 1.7 compatibility
- progress info
- `autodiff` argument in `get_right_endpoint`
- `show` method for `Endpoint`, `ParamInterval`
- add support of :FAILURE return code

## 0.4.0

- supporting gradient methods of fitting: :LD_MMA :LD_SLSQP :LD_CCSAQ
- add `get_endpoint()`, `get_interval()` methods for `scan_func::Function` to calculate Confidence Bands (CB)
- Multiple auto-test updates
- Migration from Travis to GH Actions

## 0.3.3

- NLOpt update to v0.6
- minor compatibility and CI/CD updates

## 0.3.2

- box constraints in CICO
- move Jupiter notebook to another repos insysbio/likelihoodprofiler-cases
- update docs and readme

## 0.3.1

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
