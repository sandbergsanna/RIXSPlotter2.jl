# calculates the mean square error given model data and experimental data
function calc_mse(
    model_data :: Vector{<:Real},
    exp_data  :: Vector{<:Real}
    )
    return sum((model_data-exp_data).^2)/length(model_data)
end
# export 
export calc_mse

# help function that computes the mean square error for a given set of parameters
function _compute_mse(
    lab :: LabSystem,
    param1 :: Symbol,
    param2 :: Symbol,
    param1_value :: Real,
    param2_value :: Real,
    theta_values :: Vector{<:Real},
    dQ :: Real,
    twotheta :: Real,
    lwidth :: Real,
    energy_loss_exp :: Vector{<:Vector{<:Real}},
    intensities_exp_norm :: Vector{<:Vector{<:Real}}
    )
     # set parameters
    set_parameter!(lab,param1,param1_value,site=:all,recalculate=true)
    set_parameter!(lab,param2,param2_value,site=:all,recalculate=true)
    recalculate!(lab)
    # prepare theory intensities
    intensities_theory=[similar(energy_loss_exp[i]) for i in eachindex(theta_values)]
    for i in eachindex(theta_values)
        # set scattering angles
        set_scattering_angles_deg!(lab, theta_values[i],twotheta, dQ)
        recalculate_dipole_operators!(lab)
        # get spectrum
        spectrum=get_spectrum(lab,;linewidth = lwidth)
        for j in eachindex(energy_loss_exp[i])
            intensities_theory[i][j]=intensity(spectrum,1000*energy_loss_exp[i][j])
        end
    end
    # maximum value
    max_theory=maximum([maximum(intensities_theory[i]) for i in eachindex(intensities_theory)])
    intensities_theory_norm=intensities_theory./max_theory
    # calc mean square error
    mean_square_error=sum([calc_mse(intensities_theory_norm[i],intensities_exp_norm[i]) for i in eachindex(theta_values)])/length(theta_values)
    #return
    return mean_square_error
end

# compute mean square error, option the parallelise
function compute_mse(
    lab :: LabSystem,
    param1 :: Symbol,
    param2 :: Symbol,
    param1_vec :: Vector{<:Real},
    param2_vec :: Vector{<:Real},
    theta_values :: Vector{<:Real},
    dQ :: Real,
    twotheta :: Real,
    lwidth :: Real,
    energy_loss_exp :: Vector{<:Vector{<:Real}},
    intensities_exp_norm :: Vector{<:Vector{<:Real}};
    parallel :: Bool = false
    )
    # Set up mean square error matrix
    mean_square_error=zeros(length(param2_vec),length(param1_vec))
    # maybe parallellize
    if parallel
        # check how many threads available
        if nthreads()==1
            println("Note: Only 1 thread available!")
        end
        n_blas = BLAS.get_num_threads() # save the current BLAS thread count so we can restore it after
        BLAS.set_num_threads(1) # set the number of threads for LinearAlgebra, to avoid oversubscription
        # flatten (k,m) pairs into a single task list
        tasks = [(k,m) for k in eachindex(param1_vec) for m in eachindex(param2_vec)]
        # sweep parameters
       @threads for t in eachindex(tasks)
            k, m = tasks[t]
            lab_t = deepcopy(lab)  # create a copy of lab
            # compute mean square errors
            mean_square_error[m,k]=_compute_mse(lab_t,param1,param2,param1_vec[k],param2_vec[m],theta_values,dQ,twotheta,lwidth,energy_loss_exp,intensities_exp_norm)
        end
        # restore BLAS threads so we don't affect other code outside this function
        BLAS.set_num_threads(n_blas)  
    else
        # save original parameters
        param1_init=get_parameter(lab.hamiltonian,param1,site=:all)
        param2_init=get_parameter(lab.hamiltonian,param2,site=:all)
        # sweep parameters
        for k in eachindex(param1_vec)
            for m in eachindex(param2_vec)
                # compute mean square errors
                mean_square_error[m,k]=_compute_mse(lab,param1,param2,param1_vec[k],param2_vec[m],theta_values,dQ,twotheta,lwidth,energy_loss_exp,intensities_exp_norm)
            end
        end
        # reset parameters
        set_parameter!(lab, param1,  param1_init, site=:all)
        set_parameter!(lab, param2,  param2_init, site=:all)
    end
    # return mean square errors
    return mean_square_error
end
# export
export compute_mse

# plotting help function for plotting mean square error
function _plot_mse(
        mean_square_error :: Matrix{<:Real}, 
        parameter1_vec :: Vector{<:Real},
        parameter2_vec :: Vector{<:Real};
        digits1 :: Int64 = 0,
        digits2 :: Int64 = 0,
        vmin = nothing,
        vmax = nothing
    )
        heatmap=imshow(mean_square_error,cmap="Purples",origin="lower",vmin=vmin,vmax=vmax)
        colorbar(heatmap)

        # set the x and y ticks
        xticks(collect(range(0, stop=length(parameter1_vec)-1, step=2)),[parameter1_vec[i] for i in range(1, stop=length(parameter1_vec), step=2)])
        yticks(collect(range(0, stop=length(parameter2_vec)-1, step=2)),[parameter2_vec[i] for i in range(1, stop=length(parameter2_vec), step=2)])
    
end

# calculating and plotting mean square error
function plot_mse(
    lab :: LabSystem,
    param1 :: Symbol,
    param2 :: Symbol,
    param1_vec :: Vector{<:Real},
    param2_vec :: Vector{<:Real},
    theta_values :: Vector{<:Real},
    dQ :: Real,
    twotheta :: Real,
    lwidth :: Real,
    energy_loss_exp :: Vector{<:Vector{<:Real}},
    intensities_exp_norm :: Vector{<:Vector{<:Real}};
    parallel :: Bool = false,
    vmin = nothing,
    vmax = nothing
    )
    # compute mean square errors
    mean_square_error=compute_mse(lab,param1,param2,param1_vec,param2_vec,theta_values,dQ,twotheta,lwidth,energy_loss_exp,intensities_exp_norm,parallel=parallel)
    # plot mean square error
    _plot_mse(mean_square_error,param1_vec,param2_vec,vmin=vmin,vmax=vmax)
end
# export 
export plot_mse

function plot_mse_animation(
    lab :: LabSystem,
    param1 :: Symbol,
    param2 :: Symbol,
    param3 :: Symbol,
    param1_vec :: Vector{<:Real},
    param2_vec :: Vector{<:Real},
    param3_vec :: Vector{<:Real},
    theta_values :: Vector{<:Real},
    dQ :: Real,
    twotheta :: Real,
    lwidth :: Real,
    energy_loss_exp :: Vector{<:Vector{<:Real}},
    intensities_exp_norm :: Vector{<:Vector{<:Real}};
    parallel :: Bool = false,
    vmin = nothing,
    vmax = nothing,
    output_file :: String = "mse_animation.gif",  # ".gif", does not support mp4 at the moment
    fps :: Int = 1
    )
    # save original param3
    param3_init = get_parameter(lab.hamiltonian, param3, site=:all)
    # pre-compute all MSE matrices first 
    println("Computing MSE matrices...")
    mse_frames = Vector{Matrix{Float64}}(undef, length(param3_vec))
    for (n, p3) in enumerate(param3_vec)
        println("  Frame $n / $(length(param3_vec)): $param3 = $p3")
        set_parameter!(lab, param3, p3; site=:all, recalculate=true)
        # compute MSE matrix directly (without plotting)
        mse_frames[n] = compute_mse(
            lab, param1, param2,
            param1_vec, param2_vec,
            theta_values, dQ, twotheta, lwidth,
            energy_loss_exp, intensities_exp_norm;
            parallel=parallel
        )
    end

    # restore original param3
    set_parameter!(lab, param3, param3_init; site=:all, recalculate=true)

    # determine color limits across all frames if not specified
    global_vmin = isnothing(vmin) ? minimum(minimum.(mse_frames)) : vmin
    global_vmax = isnothing(vmax) ? maximum(maximum.(mse_frames)) : vmax

    println("Building animation...")
    # set up figure
    fig, ax = subplots(figsize=(6,5))
    # initial frame
    im = ax.imshow(
        mse_frames[1],
        origin="lower",
        aspect="auto",
        extent=[param1_vec[1], param1_vec[end], param2_vec[1], param2_vec[end]],
        vmin=global_vmin, vmax=global_vmax,
        cmap="Purples"
    )
    colorbar(im, ax=ax)
    ax.set_xlabel(String(param1))
    ax.set_ylabel(String(param2))
    title_obj = ax.set_title("$param3 = $(round(param3_vec[1], digits=3))")
    tight_layout()

    # update function called for each frame
    function update(n)
        im.set_data(mse_frames[n+1])  # +1 because Python is 0-indexed
        title_obj.set_text("$param3 = $(round(param3_vec[n+1], digits=3))")
        return (im, title_obj)
    end

    # import animation
    animation = pyimport("matplotlib.animation")
    # create animation
    anim = animation.FuncAnimation(
        fig,
        update,
        frames=length(param3_vec),
        interval=1000÷fps,  # milliseconds between frames
        blit=true
    )
    # save
    if endswith(output_file, ".gif")
        anim.save(output_file, writer="pillow", fps=fps)
    else
        println("Only .gif format supported.")
    end
    println("Animation saved to: $output_file")
    close(fig)
    return anim
end
# export 
export plot_mse_animation