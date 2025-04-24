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
    dry_air :: IdealGas{FT}
    vapor :: IdealGas{FT}
    condensation :: C
end

"""
    AtmosphereThermodynamics(FT = Oceananigans.defaults.FloatType;
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
                                    dry_air, vapor, condensation)
end

#####
##### Computations
#####

const AT = AtmosphereThermodynamics
const IG = IdealGas

@inline gas_heat_capacity(R, ig::IG) = ig.isentropic_exponent * R / ig.molar_mass
@inline dry_air_heat_capacity(td::AT)         = gas_heat_capacity(td.molar_gas_constant, td.dry_air)
@inline vapor_specific_heat_capacity(td::AT)  = gas_heat_capacity(td.molar_gas_constant, td.vapor)
@inline vapor_specific_gas_constant(td::AT)   = td.molar_gas_constant / td.vapor.molar_mass
@inline dry_air_specific_gas_constant(td::AT) = td.molar_gas_constant / td.dry_air.molar_mass

const NonCondensingAtmosphereThermodynamics{FT} = AtmosphereThermodynamics{FT, Nothing} 

@inline function saturation_vapor_pressure(T, thermodynamics_constants)
    ℒ₀  = thermodynamics_constants.condensation.reference_vaporization_enthalpy
    T₀  = thermodynamics_constants.condensation.energy_reference_temperature
    Tᵗʳ = thermodynamics_constants.condensation.triple_point_temperature
    pᵗʳ = thermodynamics_constants.condensation.triple_point_pressure
    cᵖℓ = thermodynamics_constants.condensation.liquid_heat_capacity
    cᵖv = vapor_specific_heat_capacity(thermodynamics_constants)
    Rᵛ  = vapor_specific_gas_constant(thermodynamics_constants)
    Δc  = cᵖℓ - cᵖv

    return pᵗʳ * (T / Tᵗʳ)^(Δc / Rᵛ) * exp((ℒ₀ - Δc * T₀) * (1/Tᵗʳ - 1/T) / Rᵛ)
end

@inline function saturation_specific_humidity(T, ρ, td)
    p★ = saturation_vapor_pressure(T, td)
    Rᵛ = vapor_specific_gas_constant(td)
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

struct MoistPrognosticState{FT}
    ρ :: FT
    θ :: FT
    q :: FT
end

struct MoistAuxiliaryState{FT}
    p :: FT
    T :: FT
    q :: FT
end

@inline function mixture_specific_gas_constant(q, thermo::AT)
    Rᵈ = dry_air_specific_gas_constant(thermo)   
    Rᵛ = vapor_specific_gas_constant(thermo)   
    return Rᵈ * (1 - q) + Rᵛ * q
end

@inline mixture_specific_gas_constant(aux::MoistAuxiliaryState, thermo::AT) =
    mixture_specific_gas_constant(aux.q, thermo)

function density(aux::MoistAuxiliaryState, thermo)
    Rᵐ = mixture_specific_gas_constant(aux, thermo)
    return aux.p / (Rᵐ * aux.T)
end

@inline function saturation_specific_humidity(aux::MoistAuxiliaryState, td)
    ρ = density(aux, td)
    p★ = saturation_vapor_pressure(aux.T, td)
    Rᵛ = vapor_specific_gas_constant(td)
    return p★ / (ρ * Rᵛ * aux.T)
end
