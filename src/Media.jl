module Media

using Convex

"""
    MSub(ϵ_r, h, t, σ, tan_δ)

Represents a microstrip substrate with relative dielectric ϵ_r, substrate height
h (in meters), conductor thickness t (in meters), metal conductivity σ, and
loss tangent tan_δ.
"""
struct MSub
    ϵ_r::Real
    h::Real
    t::Real
    σ::Real
    tan_δ::Real
end

"""
Utility function for microstrip calculations.
"""
function Z_01(u::Real, substrate::MSub)
    η_0 = 376.73
    f = 6 + (2*π - 6)*exp(-(30.666/u)^0.7528)
    η_0/(2*π)*log(f/u + sqrt(1 + 4/u^2))
end

"""
    ϵ_eff(f, w, substrate::MSub)

Calculates the effective relative permittivity of a microstrip line of width w
(in meters) at frequency f (in Hz) on a given substrate.
"""
function ϵ_eff(f::Real, w::Real, substrate::MSub)
    u = w / substrate.h
    t = substrate.t / substrate.h
    η_0 = 376.73
    fg = f * 1e-9
    bs = substrate.h * 1e3

    Δu1 = t / π * log(1 + 4 * ℯ / t * tanh(sqrt(6.517 * u))^2)
    Δur = (1/2) * (1 + 1/cosh(sqrt(substrate.ϵ_r - 1))) * Δu1
    u1 = u + Δu1
    ur = u + Δur

    a = 1 + (1/49) * log((ur^4 + ur^2/2704) / (ur^4 + 0.432)) + (1/18.7) * log(1 + (ur/18.1)^3)
    b = 0.564 * ((substrate.ϵ_r - 0.9) / (substrate.ϵ_r + 3))^0.053
    y = (substrate.ϵ_r + 1) / 2 + (substrate.ϵ_r - 1) / 2 * (1 + 10/ur)^(-a * b)

    ϵ_eff = y * (Z_01(u1, substrate) / Z_01(ur, substrate))^2

    P_1 = 0.27488 + u*(0.6315 + 0.525*(1 + 0.0157*fg*bs)^(-20)) - 0.065683*exp(-8.7513*u)
    P_2 = 0.33622*(1 - exp(-0.03442 * substrate.ϵ_r))
    P_3 = 0.0363*exp(-4.6*u) * (1 - exp(-(fg*bs/38.7)^4.97))
    P_4 = 1 + 2.751 * (1 - exp(-(substrate.ϵ_r/15.916)^8))
    P = P_1 * P_2 * (fg*bs * (0.1844 + P_3 * P_4))^1.5763

    substrate.ϵ_r - (substrate.ϵ_r - ϵ_eff)/(1 + P)
end

"""
    Z_0(f, w, substrate::MSub)

Calculates the characterisitc impedance of a microstrip line of width w (in
meters) at frequency f (in Hz) on a given substrate.
"""
function Z_0(f::Real, w::Real, substrate::MSub)
    u = w / substrate.h
    t = substrate.t / substrate.h
    η_0 = 376.73
    fg = f * 1e-9
    bs = substrate.h * 1e3

    Δu1 = t / π * log(1 + 4 * exp(1) / t * tanh(sqrt(6.517 * u))^2)
    Δur = (1/2) * (1 + 1/cosh(sqrt(substrate.ϵ_r - 1))) * Δu1
    u1 = u + Δu1
    ur = u + Δur

    a = 1 + (1/49) * log((ur^4 + ur^2/2704) / (ur^4 + 0.432)) + (1/18.7) * log(1 + (ur/18.1)^3)
    b = 0.564 * ((substrate.ϵ_r - 0.9) / (substrate.ϵ_r + 3))^0.053
    y = (substrate.ϵ_r + 1) / 2 + (substrate.ϵ_r - 1) / 2 * (1 + 10/ur)^(-a * b)

    ϵ_eff = y * (Z_01(u1, substrate) / Z_01(ur, substrate))^2

    P_1 = 0.27488 + u*(0.6315 + 0.525*(0.0157*fg*bs + 1)^(-20)) - 0.065683*exp(-8.7513*u)
    P_2 = 0.33622*(1 - exp(-0.03442 * substrate.ϵ_r))
    P_3 = 0.0363*exp(-4.6*u) * (1 - exp(-(fg*bs/38.7)^4.97))
    P_4 = 1 + 2.751 * (1 - exp(-(substrate.ϵ_r/15.916)^8))
    P = P_1 * P_2 * (fg*bs * (0.1844 + P_3 * P_4))^1.5763

    ϵ_eff_f = substrate.ϵ_r - (substrate.ϵ_r - ϵ_eff)/(1 + P)

    Z_01(ur, substrate) / sqrt(y) * sqrt(ϵ_eff / ϵ_eff_f) * (ϵ_eff_f - 1) / (ϵ_eff - 1)
end

"""
    mline_taper_props(f, w, substrate::MSub)

Calculates the characterisitc impedance of a microstrip line of width w (in
meters) at frequency f (in Hz) on a given substrate.
"""
function mline_taper_props(f::Real, w::IndexAtom, l::Real, substrate::MSub)
    u = w / substrate.h
    t = substrate.t / substrate.h
    η_0 = 376.73
    fg = f * 1e-9
    bs = substrate.h * 1e3

    tanh_arg = sqrt(6.517 * u)
    tanh_eval = (exp(2*tanh_arg) - 1) * invpos(exp(2*tanh_arg) + 1)
    cosh_arg = sqrt(substrate.ϵ_r - 1)
    cosh_eval = (exp(cosh_arg) + exp(-cosh_arg)) / 2
    Δu1 = t / π * log(1 + 4 * ℯ / t * square(tanh_eval))
    Δur = (1/2) * (1 + 1/cosh_eval) * Δu1
    u1 = u + Δu1
    ur = u + Δur

    a = 1 + (1/49) * log((ur^4 + ur^2/2704) / (ur^4 + 0.432)) + (1/18.7) * log(1 + (ur/18.1)^3)
    b = 0.564 * ((substrate.ϵ_r - 0.9) / (substrate.ϵ_r + 3))^0.053
    y = (substrate.ϵ_r + 1) / 2 + (substrate.ϵ_r - 1) / 2 * (1 + 10/ur)^(-a * b)

    ϵ_eff = y * (Z_01(u1, substrate) / Z_01(ur, substrate))^2

    P_1 = 0.27488 + u*(0.6315 + 0.525*(0.0157*fg*bs + 1)^(-20)) - 0.065683*exp(-8.7513*u)
    P_2 = 0.33622*(1 - exp(-0.03442 * substrate.ϵ_r))
    P_3 = 0.0363*exp(-4.6*u) * (1 - exp(-(fg*bs/38.7)^4.97))
    P_4 = 1 + 2.751 * (1 - exp(-(substrate.ϵ_r/15.916)^8))
    P = P_1 * P_2 * (fg*bs * (0.1844 + P_3 * P_4))^1.5763

    ϵ_eff_f = substrate.ϵ_r - (substrate.ϵ_r - ϵ_eff)/(1 + P)

    z_0 = Z_01(ur, substrate) / sqrt(y) * sqrt(ϵ_eff / ϵ_eff_f) * (ϵ_eff_f - 1) / (ϵ_eff - 1)

    t_d = l * sqrt(ϵ_eff_f) ./ 2.998e8

    return (ϵ_eff_f, z_0, t_d)
end

end # module
