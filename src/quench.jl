function perturbed_energies(sys::System) 

    fname = "energies/energies_" * file_name(sys)

    if isfile(fname)
        energies = vec(load_data(fname))
    else
        basis = perturbed_basis(sys) 
        energies = basis.energies
    end

    return energies
end

function perturbed_states(sys::System; save_states::Bool = false) 

    fname = "states/states_" * file_name(sys)

    if isfile(fname)
        states = load_data(fname)
    else
        basis = perturbed_basis(sys; save_states = save_states) 
        states = basis.states
    end

    return states
end

function perturbed_coefficients(sys::System) 

    fname = "coefficients/coefficients_" * file_name(sys)

    if isfile(fname)
        coefficients = vec(load_data(fname))
    else
        basis = perturbed_basis(sys) 
        coefficients = overlap(sys,basis.states)
    end

    return coefficients
end

function perturbed_basis(sys::System; save_states::Bool = false) 

    fname1 = "energies/energies_" * file_name(sys)
    fname2 = "states/states_" * file_name(sys)

    if isfile(fname1) && isfile(fname2)
        energies = vec(load_data(fname1))
        states = load_data(fname2)

        return Basis(energies,states)
    else
        x,dx = position_grid(sys)
        H = eigen(Hamiltonian(sys))
        energies = H.values
        states = (1/dx)*H.vectors

        basis = truncate_basis(sys,Basis(energies,states))
        coefficients = overlap(sys,basis.states)

        save_data(fname1,basis.energies)
        save_data("coefficients/coefficients_" * file_name(sys),coefficients)
        if save_states
            save_data(fname2,basis.states)
        end
    end

    return basis
end

function overlap(sys::System, basis::Matrix{Float64}; state::Union{Nothing,Vector{Float64}} = nothing)

    if typeof(state) == Nothing
        E,state = initial_state(sys; analytic = true)
    end

    state = vec(reshape(state,(sys.N^2,1)))

    x,dx = position_grid(sys)
    a = (state'*basis)'*dx^2

    return a
end