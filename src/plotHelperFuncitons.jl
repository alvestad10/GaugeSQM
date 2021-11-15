


function get_nr_bins(_PLoop)    
    dx = 0.05
    return floor(Int64,(maximum(_PLoop) - minimum(_PLoop))/dx)
end