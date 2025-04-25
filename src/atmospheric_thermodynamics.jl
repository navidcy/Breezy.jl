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

@inline gas_heat_capacity(R, ig::IG)          = ig.isentropic_exponent * R / ig.molar_mass
@inline dry_air_heat_capacity(td::AT)         = gas_heat_capacity(td.molar_gas_constant, td.dry_air)
@inline vapor_heat_capacity(td::AT)           = gas_heat_capacity(td.molar_gas_constant, td.vapor)
@inline vapor_gas_constant(td::AT)            = td.molar_gas_constant / td.vapor.molar_mass
@inline dry_air_gas_constant(td::AT) = td.molar_gas_constant / td.dry_air.molar_mass

const NonCondensingAtmosphereThermodynamics{FT} = AtmosphereThermodynamics{FT, Nothing} 

@inline function saturation_vapor_pressure(T, thermodynamics_constants)
    @show T
    ℒ₀  = thermodynamics_constants.condensation.reference_vaporization_enthalpy
    T₀  = thermodynamics_constants.condensation.energy_reference_temperature
    Tᵗʳ = thermodynamics_constants.condensation.triple_point_temperature
    pᵗʳ = thermodynamics_constants.condensation.triple_point_pressure
    cᵖℓ = thermodynamics_constants.condensation.liquid_heat_capacity
    cᵖv = vapor_heat_capacity(thermodynamics_constants)
    Rᵛ  = vapor_gas_constant(thermodynamics_constants)
    Δc  = cᵖℓ - cᵖv
    T⁺ = max(T, Tᵗʳ)

    return pᵗʳ * (T⁺ / Tᵗʳ)^(Δc / Rᵛ) * exp((ℒ₀ - Δc * T₀) * (1/Tᵗʳ - 1/T⁺) / Rᵛ)
end

@inline function saturation_specific_humidity(T, ρ, td)
    p★ = saturation_vapor_pressure(T, td)
    Rᵛ = vapor_gas_constant(td)
    return p★ / (ρ * Rᵛ * T)
end

# p = ρ R(q) T(θ, q)
# θ(q, T)
# variable sets
# prognostic (θ, q, ρ)
#
# auxiliary
# p → T
# T → p

struct MoistState{FT}
    ρ :: FT
    θ :: FT
    q :: FT
    z :: FT
end

struct MoistAuxiliaryState{FT}
    p :: FT
    T :: FT
    q :: FT
end

@inline function specific_gas_constant(q, thermo::AT)
    Rᵈ = dry_air_gas_constant(thermo)   
    Rᵛ = vapor_gas_constant(thermo)   
    return Rᵈ * (1 - q) + Rᵛ * q
end

@inline function heat_capacity(q, thermo::AT)
    cᵖᵈ = dry_air_heat_capacity(thermo)
    cᵖᵛ = vapor_heat_capacity(thermo)
    return cᵖᵈ * (1 - q) + cᵖᵛ * q
end

function density(T, p, q, thermo)
    Rᵐ = specific_gas_constant(q, thermo)
    return p / (Rᵐ * T)
end

# Over a liquid surface
@inline function saturation_specific_humidity(T, p, q, td)
    ρ = density(T, p, q, td)
    p★ = saturation_vapor_pressure(T, td)
    Rᵛ = vapor_gas_constant(td)
    return p★ / (ρ * Rᵛ * T)
end

struct ReferenceState{FT}
    p :: FT
    θ :: FT
end

@inline function reference_pressure(z, ref, thermo)
    cᵖᵈ = dry_air_heat_capacity(thermo)
    Rᵈ = dry_air_gas_constant(thermo)
    ϰᵈ = thermo.dry_air.isentropic_exponent
    g = thermo.gravitational_acceleration
    return ref.p * (1 - g * z / (cᵖᵈ * ref.θ))^ϰᵈ
end

@inline function exner_function(prog, ref, thermo)
    Rᵐ = specific_gas_constant(prog.q, thermo)
    cᵖᵐ = heat_capacity(prog.q, thermo)
    inv_ϰᵐ = Rᵐ / cᵖᵐ
    pᵣ = reference_pressure(prog.z, ref, thermo)
    return (pᵣ / ref.p)^inv_ϰᵐ
end

# Solve
# θ = T/Π ( 1 - ℒ max(0, q - q★) / (cᵖᵐ T))
# root of: f(T) = T² - Π θ T - ℒ max(0, q - q★) / cᵖᵐ
@inline function temperature(prog::MoistState{FT}, ref, thermo::AT) where FT
    # Generate guess for unsaturated conditions
    Π = exner_function(prog, ref, thermo)
    T₁ = Π * prog.θ
    q★ = saturation_specific_humidity(T₁, prog.ρ, thermo)
    prog.q < q★ && return T₁

    # Generate guess for saturated conditions
    ℒ = thermo.condensation.reference_vaporization_enthalpy
    cᵖᵐ = heat_capacity(prog.q, thermo)
    qℓ⁺ = 1e-4 #max(0, prog.q - q★)
    T₂ = T₁ + 1e-2 #sqrt(T₁^2 + 4 * ℒ * qℓ⁺ / cᵖᵐ) / 2

    r₁ = adjustment_residual(T₁, Π, prog, thermo)
    r₂ = adjustment_residual(T₂, Π, prog, thermo)

    for i = 1:2
        # Compute slope
        dTdr = (T₂ - T₁) / (r₂ - r₁)
        # Store
        r₁ = r₂
        T₁ = r₂
        # Update
        T₂ += r₂ * dTdr
        r₂ = adjustment_residual(T₂, Π, prog, thermo)
    end

    return T₂
end

@inline function adjustment_residual(T, Π, prog, thermo)
    q★ = saturation_specific_humidity(T, prog.ρ, thermo)
    ℒ = thermo.condensation.reference_vaporization_enthalpy
    cᵖᵐ = heat_capacity(prog.q, thermo)
    return T^2 - Π * prog.θ * T - ℒ * max(0, prog.q - q★) / cᵖᵐ
end
