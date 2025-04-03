function save_data(fname,data)
    open(fname, "w") do file
        writedlm(file, data)
    end
end

load_data(fname) = readdlm(fname)

file_name(sys::System) = @sprintf("N%d_xrange%2.2f_g%2.2f_k%2.2f.dat",sys.N,sys.xrange,sys.g,sys.κ)