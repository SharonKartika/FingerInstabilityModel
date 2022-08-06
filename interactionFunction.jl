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
