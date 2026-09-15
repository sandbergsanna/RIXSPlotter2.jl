# function that plots ane energy panel plot together with the spectrum
function plot_energy_and_spectrum(
        lab         :: LabSystem,
        parameters  :: Vector{Symbol},
        param_values:: Vector{<:Vector{<:Number}},
        energies    :: Vector{<:Real},
        dq_values   :: Vector{<:Real},
        q_beam      :: Real,
        linewidth   :: Real;
        color            :: String = "k",
        plot_labels      :: Vector{<:Any} = [],
        plot_colors      :: Vector{<:Any} = [],
        plot_transitions :: Bool = false,
        parallel         :: Bool = false,
        figsize          = (12, 5)
    )
    # create figure
    fig, (ax1, ax2) = subplots(1, 2, figsize=figsize)

    # left panel: energy evolution
    show_energy_evolution_panelplot(lab.hamiltonian, parameters, param_values;subtract_GS=true, new_figure=true, color=color, parallel=parallel, ax=ax1)

    # right panel: spectrum
    # set parameters 
    for (i,param) in enumerate(parameters)
        set_parameter!(lab, param, param_values[i][end];  recalculate=true, site=:all)
    end
    I = plot_spectrum(lab, energies, dq_values, q_beam, linewidth;new_figure=true, show_figure=false,plot_labels=plot_labels, plot_colors=plot_colors,plot_transitions=plot_transitions, parallel=parallel, ax=ax2)

    tight_layout()
    show()

    return fig, I
end
# export 
export plot_energy_and_spectrum
