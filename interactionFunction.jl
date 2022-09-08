Hv(x) = float(x > 0.0)

function getforce(r)
    U0 = 2650
    U1 = 30
    U2 = 2
    U3 = 1
    A0 = 8
    A1 = 2
    A2 = 25
    A3 = 26
    force = 0
    force += U0 * r * exp(-(((r / A0)^2)))
    force += U2 * exp(-r / A2)
    force -= U3 * (r - A3)^2 * Hv(r - A3)
    force += U1 * (r - A1) * Hv(r - A1)
    return force
end


function U(r)
    U0 = 2400
    a0 = 8
    U1 = 2
    a1 = 35
    return U0 * exp(-(r / a0)^2) + U1 * (r - a1)^2 * Hv(r - a1)
end

function lenjones(r, sig=100, eps=10)
    f = 4 * eps
    f *= (12 * sig^12) * r^(-13) -
         6 * (sig^6) * r^(-7)
    return f
end

function gravity(r)
    return 1e4 / (r ^ 2)
end 

"""Border force"""
function fborder(H)
    Hmax = 1/20
    F0 = 0.
    Fmax = 1250
    if (H>0)
        return F0
    elseif (0 ≥ H ≥ -Hmax)
        return (Fmax/Hmax)*H 
    elseif (H < -Hmax)
        return Fmax 
    else 
        error("finthebug")
    end
end

