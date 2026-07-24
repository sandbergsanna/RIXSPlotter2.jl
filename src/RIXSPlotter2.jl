################################################################################
#
#   MODULE FOR PLOTTING RIXS SPECTRA
#   AUTHOR Jan Attig
#          Anna Sandberg (2.0)
#   JULIA v.1+
#
################################################################################

module RIXSPlotter2

    # using
    using PyPlot
    using LinearAlgebra
    using RIXSCalculator2
    using Base.Threads # for parallelising
    using PyCall # for animation

    # include files
    include("eigensystem_plotting.jl")
    include("lab_system_plotting.jl")
    include("spectrum_plotting.jl")
    include("spectrum_heatmap_plotting.jl")
    include("momentum_transfer_plotting.jl")
    include("mean_square_error_plotting.jl")

end
