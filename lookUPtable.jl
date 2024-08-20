function lookUPtable()

global f_minus = zeros(Int,SIZE)
global f_plus = zeros(Int,SIZE)

    for i=1:SIZE
        global f_minus[i] = round(Int,SCOPE*abs(log(1 - exp(-i/SCOPE))))
        global f_plus[i] = round(Int,SCOPE*abs(log(1 + exp(-i/SCOPE))))
    end

end

function get_f_minus(value)

    idx = ceil(Int, value)
    return f_minus[idx]

end

function get_f_plus(value)

    idx = ceil(Int, value)
    return f_plus[idx]

end