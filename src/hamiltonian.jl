struct System
    N::Int64
    xrange::Float64
    g::Float64
    κ::Float64

    System(; N = 65, xrange = 3., g = 0., κ = 0.) = new(N,xrange,g,κ)
end

mutable struct Basis
    size::Int64
    energies::Vector{Float64}
    states::Matrix{Float64}

    Basis(energies,states) = new(length(energies),energies,states)
end

function truncate_basis(sys::System, basis::Basis; threshold::Float64 = 1e-4) 

    if threshold > 0.
        coefficients = overlap(sys,basis.states)

        inds = @. abs(coefficients) > threshold
        
        return Basis(basis.energies[inds], basis.states[:,inds])
    else
        return basis
    end
end

function position_grid(N::Int64, xrange::Float64)

    # making sure that the barrier position x = 0 is contained in the grid 
    if isodd(N)
        x = collect(LinRange(-xrange,xrange,N))
    else
        x = xrange*(2*collect(1:N)/N .- 1)
    end

    dx = abs(x[2] - x[1])
    
    return [x,dx]
end

position_grid(sys::System) = position_grid(sys.N,sys.xrange)

function momentum_grid(N::Int64, xrange::Float64)

    # momentum ordering for FFT 
    dp = π/xrange
    p = dp*vcat(0:N/2-1, -N/2:-1)
    
    return [p,dp]
end

function Hamiltonian(N::Int64, xrange::Float64, g::Float64, κ::Float64)

    x,dx = position_grid(N,xrange)

    V = diagm(0.5*x.^2)
    mid = Int(ceil(N/2))
    V[mid,mid] += κ/dx

    T = [i == j ? π^2*(1 - 1/N^2)/6 : (-1)^(i-j)*(π/N)^2*cos(π*(i-j)/N)/sin(π*(i-j)/N)^2 for i = 1:N, j = 1:N]

    inds = [(l-1)*N + l for l = 1:N]
    interaction = zeros(N^2)
    interaction[inds] .= g/dx

    H = (1/dx)^2*kron(T,I(N)) .+ (1/dx)^2*kron(I(N),T) .+ kron(V,I(N)) .+ kron(I(N),V) .+ diagm(interaction)

    return H
end

Hamiltonian(sys::System) = Hamiltonian(sys.N,sys.xrange,sys.g,sys.κ)