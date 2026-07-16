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

# calc and plot spectrum
function plot_spectrum(
        lab :: LabSystem,
        energies    :: Vector{<:Real},
        dq_values   :: Vector{<:Real},
        q_beam      :: Real,
        linewidth      :: Real;
        new_figure  :: Bool = true,
        show_figure :: Bool = true,
        plot_labels  :: Vector{<:Any} =[],
        plot_colors  :: Vector{<:Any} = [],
        plot_transitions :: Bool = false,
        parallel :: Bool = false
    )

    # configure the plot
    if new_figure
        figure()
        xlabel("energy [meV]")
        ylabel("intensity [a.u.]")
    end

    # prepare vector for intensities
    I = zeros(length(energies),length(dq_values))
    spectrums=Vector{AbstractSpectrum}(undef, length(dq_values))
    # maybe parallellize
    if parallel
        # check how many threads available
        if nthreads()==1
            println("Note: Only 1 thread available!")
        end
        n_blas = BLAS.get_num_threads() # save the current BLAS thread count so we can restore it after
        BLAS.set_num_threads(1) # set the number of threads for LiearAlgebra, to avoid oversubscription
        # create one independent copy of op per thread -- set_parameter! and recalculate! mutate the operator in place
        lab_copies = [deepcopy(lab) for _ in 1:nthreads()]
        # Calculate spectrum and intensities
        @threads for i in 1:length(dq_values)
            lab_t = lab_copies[threadid()] # temporary labsystem
            set_dQ!(lab_t, dq_values[i], q_beam)
            recalculate_dipole_operators!(lab_t)
            spectrums[i]=get_spectrum(lab_t;linewidth = linewidth)
            I[:,i]=[intensity(spectrums[i], omega) for omega in energies]
        end
        # restore BLAS threads so we don't affect other code outside this function
        BLAS.set_num_threads(n_blas)
    else
        for i in 1:length(dq_values)
            set_dQ!(lab, dq_values[i], q_beam)
            recalculate_dipole_operators!(lab)
            spectrums[i]=get_spectrum(lab;linewidth = linewidth)
            I[:,i]=[intensity(spectrums[i], omega) for omega in energies]
        end
    end
    # plot intensities
    if plot_colors ==[]
        if plot_labels == []
            for i in eachindex(dq_values)
                plot(energies, I[:,i])
            end
        else
            for i in eachindex(dq_values)
                plot(energies, I[:,i], label=plot_labels[i])
            end
        end
    else
        if plot_labels == []
            for i in eachindex(dq_values)
                plot(energies, I[:,i], color=plot_colors[i])
            end
        else
            for i in eachindex(dq_values)
                plot(energies, I[:,i], label=plot_labels[i], color=plot_colors[i])
            end
        end
    end

    # add a legend
    if plot_labels != []
        legend()
    end

    # maybe plot transitions
    if plot_transitions
        for i in eachindex(dq_values)
            transitions = spectrums[i].transitions
            intensities = weight.(transitions)
            intensities ./= sum(intensities)
            intensities ./= maximum(intensities)
            for t in 1:length(transitions)
                if intensities[t] > 0.01
                    if plot_colors !=[]
                        axvline(frequency(transitions[t]), color=plot_colors[i], alpha=intensities[t])
                    else
                        axvline(frequency(transitions[t]), color="k", alpha=intensities[t])
                    end
                end
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