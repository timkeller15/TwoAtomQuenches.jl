module TwoAtomQuenches

  # for formatting strings
   using Printf 

  # for saving / loading data
  using DelimitedFiles

  # for calculating perturbed and analytic states 
  using LinearAlgebra, NonlinearSolve, SpecialFunctions, HypergeometricFunctions, FFTW

  include("hamiltonian.jl")
  export System, Basis, truncate_basis, position_grid, momentum_grid, Hamiltonian

  include("io.jl")
  export save_data, load_data, file_name

  include("analytics.jl")
  export kummerU, initial_state, Tonks_Girardeau 

  include("quench.jl")
  export perturbed_energies, perturbed_states, perturbed_coefficients, perturbed_basis, overlap 

  include("evaluation.jl")
  export loschmidt_echo, mean_loschmidt_echo, evolved_state, twodim, compute_rspdm, single_particle_density, split_step, spectral_function
end