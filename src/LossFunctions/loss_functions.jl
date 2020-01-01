# Take output of the forward map and the les data and create loss functions


"""
closure_T_nll(𝒢, les; weight = 1, subsample = 1,
                       series=false, power = 2, f1 = mean, f2 = maximum )

# Description
- Generate a loss function from the forward map and les data

# Arguments
- `𝒢` : (function), Forward map output
- `les`: (OceananigansData), les data for comparing

# Keyword Arguments
- `weight`: (scalar), multiplies the loss function by a weight
- `subsample`: (UnitRange), if time was subsampled then this should be passed in
- `series`:(boolean), whether or not to output timeseries information
- `power`:(Real), exponent in difference between temperature fields
- `f1`:(function), function to reduce spatial data at a given point in time to a number
- `f2`: (function), function to reduce temporal data from f1 into a  scalar

"""
function closure_T_nll(𝒢, les; weight = 1, subsample = 1,
                       series=false, power = 2, f1 = mean, f2 = maximum )
    if subsample != 1
        time_index = subsample
    else
        time_index = 1:length(les.t)
    end
    function nll(𝑪)
        T = 𝒢(𝑪)
        Nt = length(time_index)
        loss_value = zeros(Nt)
        for i in 1:length(time_index)
            ti = time_index[i]
            Tᵖ = T[:, i]
            Tˢ = avg( les.T[:, ti], length(Tᵖ) )
            ΔT = abs.((Tˢ - Tᵖ)).^power
            loss_value[i] = weight * f1(ΔT)
        end
        if series
            return loss_value
        else
            return f2(loss_value)
        end
    end
    return nll
end
