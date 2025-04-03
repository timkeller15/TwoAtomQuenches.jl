function loschmidt_echo(a::Vector{Float64}, energies::Vector{Float64}, t::Float64; E0::Float64 = 0.)

    phases = @. exp(1im*(E0 - energies)*t)
    L = abs2(abs2.(a)'*phases)

    return L
end

mean_loschmidt_echo(a::Vector{Float64}) = sum(abs.(a).^4)

function evolved_state(a::Vector{Float64}, states::Matrix{Float64}, energies::Vector{Float64}, t::Float64)

    factors = @. a*exp(-1im*energies*t)
    state = states*factors

    return state
end

evolved_state(a::Vector{Float64}, basis::Basis, t::Float64) = evolved_state(a,basis.states,basis.energies,t)

twodim(state::Union{Vector{Float64},Vector{ComplexF64}}, N::Int64) = length(state) == N^2 ? reshape(state,(N,N)) : nothing 

function compute_rspdm(state::Union{Vector{Float64},Vector{ComplexF64}}, N::Int64, dx::Float64)

    if length(state) == N^2
        rspdm = (twodim(state,N)*twodim(state,N)')*dx
    else
        rspdm = nothing
    end

    return rspdm
end 

function single_particle_density(state::Union{Vector{Float64},Vector{ComplexF64}}, N::Int64, dx::Float64; rspdm::Bool = false)

    if rspdm 
        density = real.(diag(compute_rspdm(state,N,dx)))
    else
        density = sum(abs2.(twodim(state,N)); dims = 2)*dx
    end

    return density
end

function split_step(sys::System, state::Union{Matrix{Float64},Matrix{ComplexF64}}, dt::Float64, timesteps::Int64)

    # kinetic energy operator
    p,dp = momentum_grid(sys.N,sys.xrange)
    K = [0.5*(p^2 + q^2) for p in p, q in p]
    Ekin = @. exp(-dt*1im*K)

    # potential operator
    x,dx = position_grid(sys.N,sys.xrange)
    V = [0.5*(x^2 + y^2) for x in x, y in x]
    mid = Int(ceil(sys.N/2))
    V[mid,:] .+= sys.κ/dx
    V[:,mid] .+= sys.κ/dx
    V .+= (sys.g/dx)*I(sys.N)
    Vpot = @. exp(-0.5*dt*1im*V)
 
    for i = 1:timesteps
        # Step 1 - half-step in position space
        state = Vpot.*state

        # Step 2 - full step in momentum space
        state = ifft(Ekin.*fft(state))

        # Step 3 - half-step in position space
        state = Vpot.*state
    end

    return state
end

function spectral_function(sys::System; discrete::Bool = true, scaled::Bool = true, range::Union{Nothing,Tuple} = nothing)

    energies = perturbed_energies(sys)
    a = perturbed_coefficients(sys)

    if scaled
        factor = 1.
    else
        factor = 2*π
    end

    if discrete
        ω = energies
        A = factor*abs2.(a)
    else
        E0 = 0.
        dt = 1e-3
        steps = 2^22
        tspan = dt*(0:(steps-1))

        ν = [sum(abs2.(a).*exp.(1im*(E0 .- energies)*t)) for t in tspan]

        ω = 2*π*fftshift(fftfreq(steps,1/dt))
        A = fftshift(factor*real.(ifft(ν)))

        if typeof(range) <: Tuple{Real,Real} && range[1] < range[2]
            inds = @. range[1] < ω < range[2]
            ω = ω[inds]
            A = A[inds]
        end
    end

    return [ω,A]
end