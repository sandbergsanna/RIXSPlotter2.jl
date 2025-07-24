# Spectrum 2D heatmap plotting function, energy vs momentum transfer
function plot_spectrum_heatmap(
        lab         :: LabSystem,
        dq_values   :: Vector{<:Real},
        q_beam      :: Real,
        energies    :: Vector{<:Real},
        lwidth      :: Real;
        new_figure  :: Bool = true,
        show_figure :: Bool = true,
        dumpfile :: String=""
    )

    # configure the plot
    if new_figure
        figure()
        ylabel("Energy (meV)")
        xlabel("dq (pi)")
    end
    
    # Calculate spectrum and intensities
    intensities = zeros(length(energies),length(dq_values))
    for i in 1:length(dq_values)
        set_dQ!(lab, dq_values[i], q_beam)
        recalculate!(lab)
        spectrum=get_spectrum(lab;linewidth = lwidth)
        intensities[:,i]=[intensity(spectrum, omega) for omega in energies]
    end

    # plot the spectrum
    heatmap=imshow(intensities,cmap="Spectral_r",origin="lower")
    cbar=colorbar(heatmap)
    cbar.set_ticks(ticks=[minimum(intensities),maximum(intensities)/2,maximum(intensities)], labels=[0,0.5,1])
    cbar.ax.tick_params(labelsize=15,direction="in")

    # set the x and y ticks
    yticks(collect(range(0, stop=length(energies)-1, length=5)),[round(E,digits=1) for E in collect(range(energies[1],stop=energies[end],length=5))])
    xticks(collect(range(0, stop=length(dq_values)-1, length=5)),[round(dq, digits=2) for dq in collect(range(dq_values[1]/pi,stop=dq_values[end]/pi,length=5))])

    # tighten the layout
    #tight_layout()

    # show the plot
    if show_figure
        show()
    end
    # saving
    if dumpfile != ""
        # open file 
        f = open(dumpfile, "w")
        # write header line
        lines = split(string(lab1.hamiltonian), "\n")
        for l in lines
            print(f,"# ",l, "\n")
        end
        # write out q_beam and intensities
        print(f, "# RIXS intensity heatmap \n")
        print(f,"# RIXS intensities have been rescaled by dividing by the maximum intensity: $(maximum(intensities)) \n")
        print(f, "# linewidth: $(lwidth)\n#\n")
        print(f, "# - First column: Energy (meV) \n")
        print(f, "# - First row: dq (pi) \n")
        print(f,"# - Remaining matrix: RIXS intensity (arb. units) \n#\n")
        # write first line
        l = "Energy/dq"
        for i in 1:length(dq_values)
            l=l*"\t$(dq_values[i]/pi)"
        end
        print(f,l,"\n")
        # write body
        for i in 1:length(energies)
            l="$(energies[i])"
            for j in 1:length(dq_values)
                l=l*"\t$(intensities[i,j]/maximum(intensities))"
            end
            print(f, l,"\n")
        end
        # close file
        close(f)
    end
end

# Spectrum 2D heatmap plotting function, energy vs momentum transfer, unit l 
function plot_spectrum_heatmap(
        lab         :: LabSystem,
        dq_values   :: Vector{<:Real},
        q_beam      :: Real,
        energies    :: Vector{<:Real},
        lwidth      :: Real,
        unit_l      :: Real;
        new_figure  :: Bool = true,
        show_figure :: Bool = true,
        dumpfile :: String=""
    )

    # configure the plot
    if new_figure
        figure()
        ylabel("Energy (meV)",fontsize=15)
        xlabel("dq (L)",fontsize=15)
    end
    
    # Calculate spectrum and intensities
    intensities = zeros(length(energies),length(dq_values))
    for i in 1:length(dq_values)
        set_dQ!(lab, dq_values[i], q_beam)
        recalculate!(lab)
        spectrum=get_spectrum(lab;linewidth = lwidth)
        intensities[:,i]=[intensity(spectrum, omega) for omega in energies]
    end

    # plot the spectrum
    heatmap=imshow(intensities,cmap="Spectral_r",origin="lower")
    cbar=colorbar(heatmap)
    cbar.set_ticks(ticks=[minimum(intensities),maximum(intensities)/2,maximum(intensities)], labels=[0,0.5,1])
    cbar.ax.tick_params(labelsize=15,direction="in")

    # set the x and y ticks
    yticks(collect(range(0, stop=length(energies)-1, length=5)),[round(E,digits=1) for E in collect(range(energies[1],stop=energies[end],length=5))])
    xticks(collect(range(0, stop=length(dq_values)-1, length=5)),[round(dq, digits=1) for dq in collect(range(dq_values[1]*unit_l/pi,stop=dq_values[end]*unit_l/pi,length=5))])

    # show the plot
    if show_figure
        show()
    end
    # saving
    if dumpfile != ""
        # open file 
        f = open(dumpfile, "w")
        # write header line
        lines = split(string(lab1.hamiltonian), "\n")
        for l in lines
            print(f,"# ",l, "\n")
        end
        # write out q_beam and intensities
        print(f, "# RIXS intensity heatmap \n")
        print(f,"# RIXS intensities have been rescaled by dividing by the maximum intensity: $(maximum(intensities)) \n")
        print(f, "# linewidth: $(lwidth)\n#\n")
        print(f, "# - First column: Energy (meV) \n")
        print(f, "# - First row: (0 0 l) (r.l.u.) \n")
        print(f,"# - Remaining matrix: RIXS intensity (arb. units) \n#\n")
        # write first line
        l = "Energy/l"
        for i in 1:length(dq_values)
            l=l*"\t$(dq_values[i]*unit_l/pi)"
        end
        print(f,l,"\n")
        # write body
        for i in 1:length(energies)
            l="$(energies[i])"
            for j in 1:length(dq_values)
                l=l*"\t$(intensities[i,j]/maximum(intensities))"
            end
            print(f, l,"\n")
        end
        # close file
        close(f)
    end
end

export plot_spectrum_heatmap
