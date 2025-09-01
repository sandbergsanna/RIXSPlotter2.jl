#######################################
#  PLOTTING against momentum transfer #
#######################################



function plot_dq_dependence(ls :: LabSystem, arg, dq_file :: AbstractString, energy_value :: Real, q_beam :: Real=15; plot_label  :: String = "", factor :: Real = 1.0, kwargs...)

    # open file and extract values
    f = open(dq_file, "r")
    lines = strip.(readlines(f))
    lines = [l for l in lines if !startswith(l, "#") && isdigit(l[1])]
    vals = [
        get_dq_eps(l, factor = factor) for l in lines
    ]
    close(f)

    # make a list of intensities
    intensities = zeros(length(vals))

    # iterate over all dq values
    for q in 1:length(vals)
        setup_dQ!(ls, vals[q][1], vals[q][2], q_beam)
        recalculate_dipole_operators!(ls)
        spectrum = get_spectrum(ls, arg; kwargs...)
        intensities[q] = intensity(spectrum, energy_value)
    end

    # plot the data
    plot(1:length(vals), intensities, label="E = $(energy_value)" * (plot_label=="" ? "" : " "*plot_label))
    xlabel("index")
    ylabel("intensity [a.u.]")
    legend()
end
function plot_dq_dependence(ls :: LabSystem, dq_file :: AbstractString, energy_value :: Real, q_beam :: Real=15; plot_label  :: String = "", factor :: Real = 1.0, kwargs...)

    # open file and extract values
    f = open(dq_file, "r")
    lines = strip.(readlines(f))
    lines = [l for l in lines if !startswith(l, "#") && isdigit(l[1])]
    vals = [
        get_dq_eps(l, factor = factor) for l in lines
    ]
    close(f)

    # make a list of intensities
    intensities = zeros(length(vals))

    # iterate over all dq values
    for q in 1:length(vals)
        setup_dQ!(ls, vals[q][1], vals[q][2], q_beam)
        recalculate_dipole_operators!(ls)
        spectrum = get_spectrum(ls; kwargs...)
        intensities[q] = intensity(spectrum, energy_value)
    end

    # plot the data
    plot(1:length(vals), intensities, label="E = $(energy_value)" * (plot_label=="" ? "" : " "*plot_label))
    xlabel("index")
    ylabel("intensity [a.u.]")
    legend()
end

function plot_dq_dependence(ls :: LabSystem, arg, dq_file :: AbstractString, energy_values :: Vector{<:Real}, q_beam :: Real=100; plot_label  :: String = "", factor :: Real = 1.0, kwargs...)

    # open file and extract values
    f = open(dq_file, "r")
    lines = strip.(readlines(f))
    lines = [l for l in lines if !startswith(l, "#") && isdigit(l[1])]
    vals = [
        get_dq_eps(l, factor = factor) for l in lines
    ]
    close(f)

    # make a list of intensities
    intensities = zeros(length(vals), length(energy_values))

    # iterate over all dq values
    for q in 1:length(vals)
        setup_dQ!(ls, vals[q][1], vals[q][2], q_beam)
        recalculate_dipole_operators!(ls)
        spectrum = get_spectrum(ls, arg; kwargs...)
        for i in 1:length(energy_values)
            intensities[q,i] = intensity(spectrum, energy_values[i])
        end
    end

    # plot the data
    for i in 1:length(energy_values)
        plot(1:length(vals), intensities[:,i], label="E = $(energy_values[i])" * (plot_label=="" ? "" : " "*plot_label))
    end
    xlabel("index")
    ylabel("intensity [a.u.]")
    legend()
end
function plot_dq_dependence(ls :: LabSystem, dq_file :: AbstractString, energy_values :: Vector{<:Real}, q_beam :: Real=100; plot_label  :: String = "", factor :: Real = 1.0, kwargs...)

    # open file and extract values
    f = open(dq_file, "r")
    lines = strip.(readlines(f))
    lines = [l for l in lines if !startswith(l, "#") && isdigit(l[1])]
    vals = [
        get_dq_eps(l, factor = factor) for l in lines
    ]
    close(f)

    # make a list of intensities
    intensities = zeros(length(vals), length(energy_values))

    # iterate over all dq values
    for q in 1:length(vals)
        setup_dQ!(ls, vals[q][1], vals[q][2], q_beam)
        recalculate_dipole_operators!(ls)
        spectrum = get_spectrum(ls; kwargs...)
        for i in 1:length(energy_values)
            intensities[q,i] = intensity(spectrum, energy_values[i])
        end
    end

    # plot the data
    for i in 1:length(energy_values)
        plot(1:length(vals), intensities[:,i], label="E = $(energy_values[i])" * (plot_label=="" ? "" : " "*plot_label))
    end
    xlabel("index")
    ylabel("intensity [a.u.]")
    legend()
end


function plot_dq_dependence(ls :: LabSystem, dq_values, energy_value :: Real, q_beam :: Real=100; plot_label  :: String = "", kwargs...)
    # make a list of intensities
    intensities = zeros(length(dq_values))

    # iterate over all dq values
    for q in 1:length(dq_values)
        set_dQ!(ls, dq_values[q], q_beam)
        recalculate_dipole_operators!(ls)
        spectrum = get_spectrum(ls; kwargs...)
        intensities[q] = intensity(spectrum, energy_value)
    end

    # plot the data
    plot(dq_values ./ pi, intensities, label="E = $(energy_value)" * (plot_label=="" ? "" : " "*plot_label))
    xlabel("dq [pi]")
    ylabel("intensity [a.u.]")
    legend()
end

function plot_dq_dependence(ls :: LabSystem, arg, dq_values, energy_value :: Real, q_beam :: Real=100; plot_label  :: String = "", kwargs...)
    # make a list of intensities
    intensities = zeros(length(dq_values))

    # iterate over all dq values
    for q in 1:length(dq_values)
        set_dQ!(ls, dq_values[q], q_beam)
        recalculate_dipole_operators!(ls)
        spectrum = get_spectrum(ls, arg; kwargs...)
        intensities[q] = intensity(spectrum, energy_value)
    end

    # plot the data
    plot(dq_values ./ pi, intensities, label="E = $(energy_value)" * (plot_label=="" ? "" : " "*plot_label))
    xlabel("dq [pi]")
    ylabel("intensity [a.u.]")
    legend()
end

function plot_dq_dependence(ls :: LabSystem, dq_values, i_from :: Int64, i_to :: Int64, q_beam :: Real=100)
    # pass on to other function
    plot_dq_dependence(ls, dq_values, ls.eigensys[:vectors][i_from], ls.eigensys[:vectors][i_to], q_beam, annotate=false)
    # format the plot
    title("intensities for <$(i_to)|D|$(i_from)>")
    xlabel("dq [pi]")
    ylabel("intensity [a.u.]")
    legend()
end


function plot_dq_dependence(ls :: LabSystem, dq_values, from :: Vector{<:Number}, to :: Vector{<:Number}, q_beam :: Real=100; annotate::Bool=true)

    # print the states
    println("Calculating dq dependence of the following MP states:\n")
    print("Initial "); printMPState(from,basis(ls)); println("")
    print("Final "); printMPState(to,basis(ls)); println("")


    # make a list of intensities
    amplitudes_h  = zeros(Complex{Float64}, length(dq_values))
    amplitudes_v  = zeros(Complex{Float64}, length(dq_values))
    intensities_h = zeros(length(dq_values))
    intensities_v = zeros(length(dq_values))
    intensities   = zeros(length(dq_values))

    # iterate over all dq values
    for q in 1:length(dq_values)
        set_dQ!(ls, dq_values[q], q_beam)
        recalculate_dipole_operators!(ls)
        # calculate amplitudes
        amplitudes_h[q]  = get_amplitude(ls.dipole_hor, from,to)
        amplitudes_v[q]  = get_amplitude(ls.dipole_ver, from,to)
        intensities_h[q] = abs.(amplitudes_h[q])^2
        intensities_v[q] = abs.(amplitudes_v[q])^2
        intensities[q]   = intensities_v[q] + intensities_h[q]
    end

    # plot the data
    plot(dq_values ./ pi, real.(amplitudes_h), label="real(A_h)", color=(0,0.2,0.5), alpha=0.4, linestyle="--")
    plot(dq_values ./ pi, imag.(amplitudes_h), label="imag(A_h)", color=(0,0.2,0.5), alpha=0.4, linestyle=":")
    plot(dq_values ./ pi,      intensities_h , label="I_h",       color=(0,0.2,0.5), alpha=1.0, linestyle="-")
    plot(dq_values ./ pi, real.(amplitudes_v), label="real(A_v)", color=(0.7,0,0.1), alpha=0.4, linestyle="--")
    plot(dq_values ./ pi, imag.(amplitudes_v), label="imag(A_v)", color=(0.7,0,0.1), alpha=0.4, linestyle=":")
    plot(dq_values ./ pi,      intensities_v , label="I_v",       color=(0.7,0,0.1), alpha=1.0, linestyle="-")
    plot(dq_values ./ pi,      intensities   , label="I",         color="k"        , alpha=1.0, linestyle="-")
    # format the plot
    if annotate
        title_string = "intensities for <TO|FROM>"
        title(title_string)
        xlabel("dq [pi]")
        ylabel("intensity [a.u.]")
        legend()
    end
end

function plot_dq_dependence(ls :: LabSystem, dq_values, energy_values :: Vector{<:Real}, q_beam :: Real=100; plot_label  :: String = "", kwargs...)
    # make a list of intensities
    intensities = zeros(length(dq_values), length(energy_values))

    # iterate over all dq values
    for q in 1:length(dq_values)
        set_dQ!(ls, dq_values[q], q_beam)
        recalculate_dipole_operators!(ls)
        spectrum = get_spectrum(ls; kwargs...)
        for i in 1:length(energy_values)
            intensities[q,i] = intensity(spectrum, energy_values[i])
        end
    end

    # plot the data
    for i in 1:length(energy_values)
        plot(dq_values ./ pi, intensities[:,i], label="E = $(energy_values[i])" * (plot_label=="" ? "" : " "*plot_label))
    end
    xlabel("dq [pi]")
    ylabel("intensity [a.u.]")
    legend()
end

function plot_dq_dependence(ls :: LabSystem, arg, dq_values, energy_values :: Vector{<:Real}, q_beam :: Real=100; plot_label  :: String = "", kwargs...)
    # make a list of intensities
    intensities = zeros(length(dq_values), length(energy_values))

    # iterate over all dq values
    for q in 1:length(dq_values)
        set_dQ!(ls, dq_values[q], q_beam)
        recalculate_dipole_operators!(ls)
        spectrum = get_spectrum(ls, arg; kwargs...)
        for i in 1:length(energy_values)
            intensities[q,i] = intensity(spectrum, energy_values[i])
        end
    end

    # plot the data
    for i in 1:length(energy_values)
        plot(dq_values ./ pi, intensities[:,i], label="E = $(energy_values[i])" * (plot_label=="" ? "" : " "*plot_label))
    end
    xlabel("dq [pi]")
    ylabel("intensity [a.u.]")
    legend()
end

# export the function
export plot_dq_dependence

# Function that plots the theta dependence of given multiplets
function plot_theta_dependence_multiplets(lab:: LabSystem, theta_values::Vector{<:Real}, twotheta :: Real, dQ :: Real, to_multiplets::Vector{<:Integer}; new_figure:: Bool=true,dumpfile::String="")
    # unique energies
    unique_energies,indices=multiplets(lab.eigensys)
    # calc I vs theta for multiplets
    I=[theta_dependence_multiplet(lab,theta_values,twotheta,dQ,i) for i in to_multiplets]
    # if new figure
    if new_figure
        figure()
        # labels
        ylabel("RIXS intensity (arb. units)",fontsize=15)
        xlabel(L"\Theta (deg)",fontsize=15)
        # ticks
        minorticks_on()
        tick_params(axis="both",which="both",labelsize=15,direction="in")
        tick_params(axis="both", which="major",length=5)
        tick_params(axis="both", which="minor",length=3)
    end
    # plot
    for i in eachindex(to_multiplets)
        plot(theta_values,I[i],label="E=$(round(unique_energies[to_multiplets[i]]-unique_energies[1]))")
    end
    # legend
    legend()
    # saving
    if dumpfile != ""
        # open file
        f = open(dumpfile*".txt", "w")
        # write header line hamiltonian
        lines = split(string(lab.hamiltonian), "\n")
        for l in lines
            print(f,"# ",l, "\n")
        end
        # write header with multiplet energies
        print(f, "# RIXS intensity (arb. units) for multiplets as a function of theta \n")
        print(f, "# dQ= $(dQ) \n#\n")
        el="# Multiplet energies: \n"
        for i in to_multiplets
            el=el*" \t$(round(unique_energies[i]-unique_energies[1]))"
        end
        print(f,el,"\n")
        # write header line
        hl = "# theta \t"
        for i in to_multiplets
            hl = hl * "\t I(E_$(i))"
        end
        print(f, hl, "\n")
        # write body
        for i in 1:length(theta_values)
            l = "$(theta_values[i])"
            for j in eachindex(to_multiplets)
                l = l * "\t$(I[j][i])"
            end
            print(f, l, "\n")
        end
        # close file
        close(f)
    end
end
export plot_theta_dependence_multiplets

# Function that plots the dq dependence of given multiplets
function plot_dq_dependence_multiplets(lab:: LabSystem, dq_values::Vector{<:Real}, q_beam :: Real, to_multiplets::Vector{<:Integer}; new_figure:: Bool=true,dumpfile::String="")
    # unique energies
    unique_energies,indices=multiplets(lab.eigensys)
    # calc I vs theta for multiplets
    I=[dq_dependence_multiplet(lab,dq_values,q_beam,i) for i in to_multiplets]
    # if new figure
    if new_figure
        figure()
        # labels
        ylabel("RIXS intensity (arb. units)",fontsize=15)
        xlabel("dq [pi]",fontsize=15)
        # ticks
        minorticks_on()
        tick_params(axis="both",which="both",labelsize=15,direction="in")
        tick_params(axis="both", which="major",length=5)
        tick_params(axis="both", which="minor",length=3)
    end
    # plot
    for i in eachindex(to_multiplets)
        plot(dq_values./pi,I[i],label="E=$(round(unique_energies[to_multiplets[i]]-unique_energies[1]))")
    end
    # legend
    legend()
    # saving
    if dumpfile != ""
        # open file
        f = open(dumpfile*".txt", "w")
        # write header line hamiltonian
        lines = split(string(lab.hamiltonian), "\n")
        for l in lines
            print(f,"# ",l, "\n")
        end
        # write header with multiplet energies
        print(f, "# RIXS intensity (arb. units) for multiplets as a function of theta \n")
        print(f, "# dQ= $(dQ) \n#\n")
        el="# Multiplet energies: \n"
        for i in to_multiplets
            el=el*" \t$(round(unique_energies[i]-unique_energies[1]))"
        end
        print(f,el,"\n")
        # write header line
        hl = "# theta \t"
        for i in to_multiplets
            hl = hl * "\t I(E_$(i))"
        end
        print(f, hl, "\n")
        # write body
        for i in 1:length(theta_values)
            l = "$(theta_values[i])"
            for j in eachindex(to_multiplets)
                l = l * "\t$(I[j][i])"
            end
            print(f, l, "\n")
        end
        # close file
        close(f)
    end
end
export plot_dq_dependence_multiplets
