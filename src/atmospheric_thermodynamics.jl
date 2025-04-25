struct IdealGas{FT}
    molar_mass :: FT
    isentropic_exponent :: FT # κ = γ / (γ - 1)
end

function IdealGas(FT = Oceananigans.defaults.FloatType;
                  molar_mass = 0.02897,
                  isentropic_exponent = 2/7)

    return IdealGas{FT}(convert(FT, molar_mass),
                        convert(FT, isentropic_exponent))
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
                                  dry_air_isentropic_exponent = 2/7,
                                  vapor_molar_mass = 0.018015,
                                  vapor_isentropic_exponent = 4.03,
                                  condensation = Condensation(FT))

    dry_air = IdealGas(FT; molar_mass = dry_air_molar_mass,
                           isentropic_exponent = dry_air_isentropic_exponent)

    vapor = IdealGas(FT; molar_mass = vapor_molar_mass,
                         isentropic_exponent = vapor_isentropic_exponent)

    return AtmosphereThermodynamics(convert(FT, molar_gas_constant),
                                    convert(FT, gravitational_acceleration),
                                    dry_air, vapor, condensation)
end

#####
##### Computations
#####

const AT = AtmosphereThermodynamics
const IG = IdealGas

@inline gas_heat_capacity(R, ig::IG)      = ig.isentropic_exponent * R / ig.molar_mass
@inline dry_air_heat_capacity(thermo::AT) = gas_heat_capacity(thermo.molar_gas_constant, thermo.dry_air)
@inline vapor_heat_capacity(thermo::AT)   = gas_heat_capacity(thermo.molar_gas_constant, thermo.vapor)
@inline vapor_gas_constant(thermo::AT)    = thermo.molar_gas_constant / thermo.vapor.molar_mass
@inline dry_air_gas_constant(thermo::AT)  = thermo.molar_gas_constant / thermo.dry_air.molar_mass

const NonCondensingAtmosphereThermodynamics{FT} = AtmosphereThermodynamics{FT, Nothing} 

@inline function saturation_vapor_pressure(T, thermodynamics_constants)
    ℒ₀  = thermodynamics_constants.condensation.reference_vaporization_enthalpy
    T₀  = thermodynamics_constants.condensation.energy_reference_temperature
    Tᵗʳ = thermodynamics_constants.condensation.triple_point_temperature
    pᵗʳ = thermodynamics_constants.condensation.triple_point_pressure
    cᵖℓ = thermodynamics_constants.condensation.liquid_heat_capacity
    cᵖv = vapor_heat_capacity(thermodynamics_constants)
    Rᵛ  = vapor_gas_constant(thermodynamics_constants)
    Δc  = cᵖℓ - cᵖv
    # T = max(T, Tᵗʳ)
    return pᵗʳ * (T / Tᵗʳ)^(Δc / Rᵛ) * exp((ℒ₀ - Δc * T₀) * (1/Tᵗʳ - 1/T) / Rᵛ)
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
    cᵖᵈ = dry_air_heat_capacity(thermo)
    cᵖᵛ = vapor_heat_capacity(thermo)
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
    p₀ :: FT # reference pressure at z=0
    θ :: FT  # constant reference potential temperature
end

"""
    reference_density(ref, thermo)

Compute the reference density associated with the reference pressure and potential temperature.
The reference density is defined as the density of dry air at the reference pressure and temperature.
"""
@inline function reference_density(ref, thermo)
    Rᵈ = dry_air_gas_constant(thermo)
    return ref.p₀ / (Rᵈ * ref.θ)
end

@inline function reference_pressure(z, ref, thermo)
    cᵖᵈ = dry_air_heat_capacity(thermo)
    Rᵈ = dry_air_gas_constant(thermo)
    ϰᵈ = thermo.dry_air.isentropic_exponent
    g = thermo.gravitational_acceleration
    return ref.p₀ * (1 - g * z / (cᵖᵈ * ref.θ))^ϰᵈ
end

@inline function saturation_specific_humidity(T, ref::ReferenceState, thermo)
    ρ₀ = reference_density(ref, thermo)
    return saturation_specific_humidity(T, ρ₀, thermo)
end

@inline function exner_function(state, ref, thermo)
    Rᵐ = mixture_gas_number(state.q, thermo)
    cᵖᵐ = mixture_heat_capacity(state.q, thermo)
    inv_ϰᵐ = Rᵐ / cᵖᵐ
    pᵣ = reference_pressure(state.z, ref, thermo)
    return (pᵣ / ref.p₀)^inv_ϰᵐ
end

function condensate_specific_humidity(T, q, ref, thermo)
    q★ = saturation_specific_humidity(T, ref, thermo)
    return max(0, q - q★)
end

# Solve
# θ = T/Π ( 1 - ℒ qℓ / (cᵖᵐ T))
# for temperature T.
# root of: f(T) = T² - Π θ T - ℒ max(0, q - q★) / cᵖᵐ
@inline function compute_temperature(state::ThermodynamicState{FT}, ref, thermo) where FT
    # Generate guess for unsaturated conditions
    Π = exner_function(state, ref, thermo)
    T₁ = Π * state.θ
    qˡ₁ = condensate_specific_humidity(T₁, state.q, ref, thermo)
    qˡ₁ <= 0 && return T₁

    # If we're here, we have condensation
    r₁ = adjustment_residual(T₁, Π, qˡ₁, state, thermo)

    ℒ = thermo.condensation.reference_vaporization_enthalpy
    cᵖᵐ = mixture_heat_capacity(state.q, thermo)
    T₂ = (T₁ + sqrt(T₁^2 - 4 * ℒ * qˡ₁ / cᵖᵐ)) / 2
    qˡ₂ = condensate_specific_humidity(T₂, state.q, ref, thermo)

    r₂ = adjustment_residual(T₂, Π, qˡ₂, state, thermo)

    while r₂ > 1e-4
        # Compute slope
        ΔTΔr = (T₂ - T₁) / (r₂ - r₁)
        # Store
        r₁ = r₂
        T₁ = T₂
        # Update
        T₂ -= r₂ * ΔTΔr
        qˡ₂ = condensate_specific_humidity(T₂, state.q, ref, thermo)
        r₂ = adjustment_residual(T₂, Π, qˡ₂, state, thermo)
    end

    return T₂
end

@inline function adjustment_residual(T, Π, qˡ, state, thermo)
    ℒ = thermo.condensation.reference_vaporization_enthalpy
    cᵖᵐ = mixture_heat_capacity(state.q, thermo)
    return T^2 - ℒ * qˡ / cᵖᵐ - Π * state.θ * T
end
