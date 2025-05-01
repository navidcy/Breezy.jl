struct IdealGas{FT}
    molar_mass :: FT
    heat_capacity :: FT # specific heat capacity at constant pressure
end

function IdealGas(FT = Oceananigans.defaults.FloatType;
                  molar_mass = 0.02897,
                  heat_capacity = 1005)

    return IdealGas{FT}(convert(FT, molar_mass),
                        convert(FT, heat_capacity))
end

struct Condensation{FT}
    reference_vaporization_enthalpy :: FT
    energy_reference_temperature :: FT
    triple_point_temperature :: FT
    triple_point_pressure :: FT
    liquid_heat_capacity :: FT
end

function Condensation(FT = Oceananigans.defaults.FloatType;
                      reference_vaporization_enthalpy = 2500800,
                      energy_reference_temperature = 273.16,
                      liquid_heat_capacity = 4181,
                      triple_point_temperature = 273.16,
                      triple_point_pressure = 611.657)

    return Condensation{FT}(convert(FT, reference_vaporization_enthalpy),
                            convert(FT, energy_reference_temperature),
                            convert(FT, triple_point_temperature),
                            convert(FT, triple_point_pressure),
                            convert(FT, liquid_heat_capacity))
end

struct AtmosphereThermodynamics{FT, C}
    molar_gas_constant :: FT
    gravitational_acceleration :: FT
    dry_air :: IdealGas{FT}
    vapor :: IdealGas{FT}
    condensation :: C
end

"""
    AtmosphereThermodynamics(FT = Oceananigans.defaults.FloatType;
                             gravitational_acceleration = 9.81,
                             molar_gas_constant = 8.314462618,
                             dry_air_molar_mass = 0.02897,
                             dry_air_isentropic_exponent = 2/7,
                             vapor_molar_mass = 0.018015,
                             vapor_isentropic_exponent = 4.03,
                             condensation = nothing)

Create `AtmosphereThermodynamics` with parameters that correpsond to the composition of dry air
in Earth's atmosphere and water vapor. The default `isnothing(condensation)` implies thermodynamics
appropriate for unsaturated and therefore non-condensing air.
"""
function AtmosphereThermodynamics(FT = Oceananigans.defaults.FloatType;
                                  gravitational_acceleration = 9.81,
                                  molar_gas_constant = 8.314462618,
                                  dry_air_molar_mass = 0.02897,
                                  dry_air_heat_capacity = 1005,
                                  vapor_molar_mass = 0.018015,
                                  vapor_heat_capacity = 1850,
                                  condensation = Condensation(FT))

    dry_air = IdealGas(FT; molar_mass = dry_air_molar_mass,
                           heat_capacity = dry_air_heat_capacity)

    vapor = IdealGas(FT; molar_mass = vapor_molar_mass,
                         heat_capacity = vapor_heat_capacity)

    return AtmosphereThermodynamics(convert(FT, molar_gas_constant),
                                    convert(FT, gravitational_acceleration),
                                    dry_air, vapor, condensation)
end

#####
##### Computations
#####

const AT = AtmosphereThermodynamics
const IG = IdealGas

@inline vapor_gas_constant(thermo::AT)    = thermo.molar_gas_constant / thermo.vapor.molar_mass
@inline dry_air_gas_constant(thermo::AT)  = thermo.molar_gas_constant / thermo.dry_air.molar_mass

const NonCondensingAtmosphereThermodynamics{FT} = AtmosphereThermodynamics{FT, Nothing} 

"""
    saturation_vapor_pressure(T, thermo)

Compute the saturation vapor pressure over a liquid surface by integrating
the Clausius-Clapeyron relation,

```math
dp/dT = ℒᵛ / (Rᵛ T^2)
```
"""
@inline function saturation_vapor_pressure(T, thermo)
    ℒ₀  = thermo.condensation.reference_vaporization_enthalpy
    T₀  = thermo.condensation.energy_reference_temperature
    Tᵗʳ = thermo.condensation.triple_point_temperature
    pᵗʳ = thermo.condensation.triple_point_pressure
    cᵖˡ = thermo.condensation.liquid_heat_capacity
    cᵖᵛ = thermo.vapor.heat_capacity
    Rᵛ  = vapor_gas_constant(thermo)
    Δcᵖ = cᵖˡ - cᵖᵛ
    Δϰ = Δcᵖ / Rᵛ
    return pᵗʳ * (T / Tᵗʳ)^Δϰ * exp((ℒ₀ - Δcᵖ * T₀) * (1/Tᵗʳ - 1/T) / Rᵛ)
end

# Over a liquid surface
@inline function saturation_specific_humidity(T, ρ, thermo)
    p★ = saturation_vapor_pressure(T, thermo)
    Rᵛ = vapor_gas_constant(thermo)
    return p★ / (ρ * Rᵛ * T)
end

@inline function mixture_gas_number(q, thermo::AT)
    Rᵈ = dry_air_gas_constant(thermo)   
    Rᵛ = vapor_gas_constant(thermo)   
    return Rᵈ * (1 - q) + Rᵛ * q
end

"""
    mixture_heat_capacity(q, thermo)

Compute the heat capacity of state air given the total specific humidity q
and assuming that condensate mass ratio qℓ ≪ q, where qℓ is the mass ratio of
liquid condensate.
"""
@inline function mixture_heat_capacity(q, thermo::AT)
    cᵖᵈ = thermo.dry_air.heat_capacity
    cᵖᵛ = thermo.vapor.heat_capacity
    return cᵖᵈ * (1 - q) + cᵖᵛ * q
end

#=
For example, if we would like to account for the mass of condensate consistently:
struct WarmCondensate{FT}
    vapor :: FT
    liquid :: FT
    total :: FT
end

# Then more correctly:
@inline function mixture_heat_capacity(q::WarmCondensate, thermo::AT)
    cᵖᵈ = dry_air_heat_capacity(thermo)
    cᵖᵛ = vapor_heat_capacity(thermo)
    cᵖˡ = thermo.condensation.liquid_heat_capacity
    qᵗ = q.total
    qᵛ = q.vapor
    qˡ = q.liquid
    return cᵖᵈ * (1 - qᵗ) + cᵖᵛ * qᵛ + cᵖˡ * qˡ
end

@inline function mixture_gas_number(q::WarmCondensate, thermo::AT)
=#

#####
##### state thermodynamics for a Boussinesq model
#####

# Organizing information about the state is a WIP
struct ThermodynamicState{FT}
    θ :: FT
    q :: FT
    z :: FT
end

struct ReferenceState{FT}
    p₀ :: FT # base pressure: reference pressure at z=0
    θ :: FT  # constant reference potential temperature
end

function ReferenceState(FT = Oceananigans.defaults.FloatType;
                        base_pressure = 101325,
                        potential_temperature = 288)

    return ReferenceState{FT}(convert(FT, base_pressure),
                              convert(FT, potential_temperature))
end

"""
    reference_density(ref, thermo)

Compute the reference density associated with the reference pressure and potential temperature.
The reference density is defined as the density of dry air at the reference pressure and temperature.
"""
@inline function reference_density(z, ref, thermo)
    Rᵈ = dry_air_gas_constant(thermo)
    p = reference_pressure(z, ref, thermo)
    return p / (Rᵈ * ref.θ)
end

@inline function base_density(ref, thermo)
    Rᵈ = dry_air_gas_constant(thermo)
    return ref.p₀ / (Rᵈ * ref.θ)
end

@inline function reference_specific_volume(z, ref, thermo)
    Rᵈ = dry_air_gas_constant(thermo)
    p = reference_pressure(z, ref, thermo)
    return Rᵈ * ref.θ / p
end

@inline function reference_pressure(z, ref, thermo)
    cᵖᵈ = thermo.dry_air.heat_capacity
    Rᵈ = dry_air_gas_constant(thermo)
    ϰᵈ⁻¹ = Rᵈ / cᵖᵈ
    g = thermo.gravitational_acceleration
    return ref.p₀ * (1 - g * z / (cᵖᵈ * ref.θ))^ϰᵈ⁻¹
end

@inline function saturation_specific_humidity(T, z, ref::ReferenceState, thermo)
    ρ = reference_density(z, ref, thermo)
    return saturation_specific_humidity(T, ρ, thermo)
end

@inline function exner_function(state, ref, thermo)
    Rᵐ = mixture_gas_number(state.q, thermo)
    cᵖᵐ = mixture_heat_capacity(state.q, thermo)
    inv_ϰᵐ = Rᵐ / cᵖᵐ
    pᵣ = reference_pressure(state.z, ref, thermo)
    return (pᵣ / ref.p₀)^inv_ϰᵐ
end

#####
##### Saturation adjustment
#####

condensate_specific_humidity(T, state, ref, thermo) =
    condensate_specific_humidity(T, state.q, state.z, ref, thermo)

function condensate_specific_humidity(T, q, z, ref, thermo)
    qᵛ★ = saturation_specific_humidity(T, z, ref, thermo)
    return max(0, q - qᵛ★)
end

# Solve
# θ = T/Π ( 1 - ℒ qˡ / (cᵖᵐ T))
# for temperature T with qˡ = max(0, q - qᵛ★).
# root of: f(T) = T² - Π θ T - ℒ qˡ / cᵖᵐ
@inline function temperature(state::ThermodynamicState{FT}, ref, thermo) where FT
    state.θ == 0 && return zero(FT)

    # Generate guess for unsaturated conditions
    Π = exner_function(state, ref, thermo)
    T₁ = Π * state.θ
    qˡ₁ = condensate_specific_humidity(T₁, state, ref, thermo)
    qˡ₁ <= 0 && return T₁

    # If we made it this far, we have condensation
    r₁ = adjustment_residual(T₁, Π, qˡ₁, state, thermo)

    ℒ = thermo.condensation.reference_vaporization_enthalpy
    cᵖᵐ = mixture_heat_capacity(state.q, thermo)
    T₂ = (T₁ + sqrt(T₁^2 + 4 * ℒ * qˡ₁ / cᵖᵐ)) / 2
    qˡ₂ = condensate_specific_humidity(T₂, state, ref, thermo)
    r₂ = adjustment_residual(T₂, Π, qˡ₂, state, thermo)

    # Saturation adjustment
    R = sqrt(max(T₂, T₁))
    ϵ = convert(FT, 1e-4)
    δ = ϵ * R 
    iter = 0

    while abs(r₂ - r₁) > δ
        # Compute slope
        ΔTΔr = (T₂ - T₁) / (r₂ - r₁)

        # Store previous values
        r₁ = r₂
        T₁ = T₂

        # Update
        T₂ -= r₂ * ΔTΔr
        qˡ₂ = condensate_specific_humidity(T₂, state, ref, thermo)
        r₂ = adjustment_residual(T₂, Π, qˡ₂, state, thermo)
        iter += 1
    end

    return T₂
end

@inline function adjustment_residual(T, Π, qˡ, state, thermo)
    ℒ = thermo.condensation.reference_vaporization_enthalpy
    cᵖᵐ = mixture_heat_capacity(state.q, thermo)
    return T^2 - ℒ * qˡ / cᵖᵐ - Π * state.θ * T
end

@inline function specific_volume(state, ref, thermo)
    T = temperature(state, ref, thermo)
    Rᵐ = mixture_gas_number(state.q, thermo)
    pʳ = reference_pressure(state.z, ref, thermo)
    return Rᵐ * T / pʳ
end
