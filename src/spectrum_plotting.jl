################################################################################
#
#   Spectrum plotting functions
#
################################################################################


# plot a spectrum that has been calculated before
function plot_spectrum(
        spectrum    :: S,
        energies    :: Vector{<:Real};
        new_figure  :: Bool = true,
        show_figure :: Bool = true,
        plot_label  :: String = "",
        plot_color  :: Any = "b",
        plot_transitions :: Bool = false
    ) where {T,S<:AbstractSpectrum{T}}

    # configure the plot
    if new_figure
        figure()
        xlabel("energy [meV]")
        ylabel("intensity [a.u.]")
    end

    # plot the spectrum
    I=[intensity(spectrum, omega) for omega in energies]
    plot(
        energies, I,
        label=plot_label, color=plot_color
    )

    # add a legend
    if plot_label != ""
        legend()
    end

    # maybe plot transitions
    if plot_transitions
        transitions = spectrum.transitions
        intensities = weight.(transitions)
        intensities ./= sum(intensities)
        intensities ./= maximum(intensities)
        for t in 1:length(transitions)
            if intensities[t] > 0.01
                axvline(frequency(transitions[t]), color="k", alpha=intensities[t])
            end
        end
    end

    # set the x limits
    xlim(energies[1], energies[end])
    ylim(0, ylim()[2])

    # tighten the layout
    tight_layout()

    # show the plot
    if show_figure
        show()
    end
    # return RIXS intensities
    return I
end

# export functions
export plot_spectrum

# function save spectrum
function save_spectrum_theta(
    lab :: LabSystem,
    energies    :: Vector{<:Real},
    intensities :: Vector{<:Real},
    theta:: Real, 
    twotheta :: Real, 
    dQ :: Real,
    filename :: String
)
    # open file
    f = open(filename*".txt", "w")
    # write header line hamiltonian
    lines = split(string(lab.hamiltonian), "\n")
    for l in lines
        print(f,"# ",l, "\n")
    end
    # write header with multiplet energies
    print(f, "# RIXS intensity (arb. units) as a function of energy loss \n")
    print(f, "# theta=$(theta), twotheta=$(twotheta), dQ= $(dQ) \n#\n")
    # write header line
    hl = "# E \t I"
    print(f, hl, "\n")
    # write body
    for i in 1:length(energies)
        l = "$(energies[i])\t$(intensities[i])"
        print(f, l, "\n")
    end
    # close file
    close(f)  
end
# export functions
export save_spectrum_theta