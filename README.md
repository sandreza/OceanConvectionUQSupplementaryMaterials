# LocalOceanUQSupplementaryMaterials
Supplementary code and information

## Folders
* [LES](#les)
* [figure_scripts](#figure_scripts)
* [mcmc_data](#mcmc_data)
* [scripts](#scripts)
* [scrc](#scr)
* [tests](#tests)

## LES
JLD2 files of horizontally averaged statistics from the Boussinesq equations.

## figure_scripts
Scripts for generating figures.

## mcmc_data
JLD2 files that contain MCMC runs for various cases.

## scripts
Scripts for generating data for figures. This includes MCMC data and optimization data.

## src
Contains the "Core Functionality" module. This contains several convenient functions such as
* Converting [Oceananigans](https://github.com/climate-machine/Oceananigans.jl) JLD2 files into a struct that makes interfacing with [OceanTurb](https://github.com/glwagner/OceanTurb.jl) easy.
* Implementation of a RWMCMC method
* Implementation of loss functions and forward map concepts for [OceanTurb](https://github.com/glwagner/OceanTurb.jl).

## Tests
Functions for testing functionality in src.
