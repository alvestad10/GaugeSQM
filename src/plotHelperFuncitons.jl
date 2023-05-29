


function get_nr_bins(_PLoop)    
    dx = 0.05
    PLoop = copy(_PLoop)
    PLoop[isnan.(PLoop) .|| abs.(PLoop) .> 1000] .= 100.
    return floor(Int64,(maximum(PLoop) - minimum(PLoop))/dx)
end