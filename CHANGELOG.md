# Change Log

## 0.3.0

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
