# eigenvalue equation for two interacting atoms
f(x,p) = -2*sqrt(2)*gamma(-0.5*x + 0.75)/gamma(-0.5*x + 0.25) - p 

kummerU(a,b,z) = abs(z) < 25 ? gamma(1-b)*pFq((1,a),(1,b),z)/gamma(a + 1 - b) + gamma(b-1)*z^(1-b)*pFq((1,a+1-b),(1,2-b),z)/gamma(a) : 0.

function initial_state(N::Int64, xrange::Float64, g::Float64; analytic::Bool = true)

    x,dx = position_grid(N,xrange)
    
    if analytic
        uspan = (-0.5*g^2, 1.5)
        energy = getfield(solve(IntervalNonlinearProblem(f,uspan,g)),:u) 

        state = [exp(-0.5*(x^2 + y^2))*kummerU(0.25-0.5*energy,0.5,0.5*x^2 - x*y + 0.5*y^2) for x in x, y in x]
        state /= sqrt(sum(abs2.(state))*dx^2)

        energy += 0.5
    else
        sys = System(N = N, xrange = xrange, g = g, κ = 0.)
        basis = perturbed_basis(sys)
        energy = basis.energies[1]
        state = basis.states[:,1]
    end

    return [energy,state]
end

initial_state(sys::System; analytic::Bool = true) = initial_state(sys.N,sys.xrange,sys.g; analytic = analytic)

function Tonks_Girardeau(sys::System; levels::Int64 = 50, N::Int64 = 2049)

    energies = zeros(levels)
    coefficients = zeros(levels)

    x,dx = position_grid(N,sys.xrange)
    ϕ = (1/sqrt(sqrt(π)))*[exp(-0.5*x^2) for x in x]

    # only even states are affected by barrier and have non-zero coefficients
    for (i,n) in enumerate(2*(0:(levels-1)))

        if n == 0
            uspan = (-sys.κ^2, 1.49)
        else
            uspan = n .+ (-0.49, 1.49)
        end

        energy = getfield(solve(IntervalNonlinearProblem(f,uspan,sqrt(2)*sys.κ)),:u) 
        state = [exp(-0.5*x^2)*kummerU(0.25-0.5*energy,0.5,x^2) for x in x]
        state /= sqrt(sum(abs2.(state))*dx)

        energies[i] = energy
        coefficients[i] = (ϕ'*state)*dx
    end

    return [energies,coefficients]
end