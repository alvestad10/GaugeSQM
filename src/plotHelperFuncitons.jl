


function get_nr_bins(_PLoop)    
    dx = 0.01
    return floor(Int64,(maximum(_PLoop) - minimum(_PLoop))/dx)
end