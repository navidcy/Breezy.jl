struct IdealGasThermodynamicsConstants{FT}
    gas_constant :: FT
    molar_mass :: FT
    isentropic_exponent :: FT
end

function IdealGasThermodynamicsConstants(FT = Oceananigans.defaults.FloatType;
                                        gas_constant = 8.3144598,
                                        isentropic_exponent = 2/7,
                                        molar_mass = 0.02897)

    return IdealGasThermodynamicsConstants{FT}(convert(FT, isentropic_exponent),
                                              convert(FT, molar_mass))
end

struct VaporThermodynamicsConstants{FT}
    molar_mass :: FT
    heat_capacity  :: FT
end

"""
    VaporThermodynamicsConstants

Build thermodynamic constants for a condensable vapor component in a mixture of gases.
The defaults apply to water vapor.
"""
function VaporThermodynamicsConstants(FT = Oceananigans.defaults.FloatType;
                                     molar_mass   = 1859,
                                     heat_capacity = 0.018015)

    return GasConstituentThermodynamicsConstants{FT}(convert(FT, molar_mass),
                                                    convert(FT, heat_capacity))
end

struct CondensationThermodynamicssConstants{FT}
    reference_vaporization_enthalpy :: FT
    energy_reference_temperature :: FT
    triple_point_temperature :: FT
    triple_point_pressure :: FT
    liquid_heat_capacity :: FT
end

function CondensationThermodynamicssConstants(FT = Oceananigans.defaults.FloatType;
                                             reference_vaporization_enthalpy = 2500800,
                                             energy_reference_temperature = 273.16,
                                             water_vapor_heat_capacity = 1859,
                                             triple_point_temperature = 273.16,
                                             triple_point_pressure = 611.657)

    return CondensationThermodynamicssConstants{FT}(convert(FT, reference_vaporization_enthalpy),
                                                   convert(FT, energy_reference_temperature),
                                                   convert(FT, water_vapor_heat_capacity),
                                                   convert(FT, triple_point_temperature),
                                                   convert(FT, triple_point_pressure))
end


const DATC = IdealGasThermodynamicsConstants
const WVTC = WaterVaporThermodynamicsConstants
const CTC =  CondensationThermodynamicssConstants

Base.eltype(::DATC{FT}) where FT = FT
Base.eltype(::WVTC{FT}) where FT = FT
Base.eltype(::CTC{FT})  where FT = FT

struct ThermodynamicsConstants{FT, C}
    dry_air :: IdealGasThermodynamicsConstants{FT}
    water_vapor :: VaporThermodynamicsConstants{FT}
    condensation :: C
end

vapor_specific_gas_constant(tc::ThermodynamicsConstants) = tc.dry_air.gas_constant / tc.water_vapor.molar_mass

const DryAirThermodynamicsConstants{FT} = DryAirThermodynamicsConstants{FT, Nothing} 

@inline function saturation_vapor_pressure(T, thermodynamics_constants)
    ℒ₀  = thermodynamics_constants.condesation.reference_vaporization_enthalpy
    T₀  = thermodynamics_constants.condesation.energy_reference_temperature
    Tᵗʳ = thermodynamics_constants.condesation.triple_point_temperature
    pᵗʳ = thermodynamics_constants.condesation.triple_point_pressure
    cᵖℓ = thermodynamics_constants.condesation.liquid_heat_capacity
    cᵖv = thermodynamics_constants.water_vapor.heat_capacity
    Rᵛ  = vapor_specific_gas_constant(thermodynamics_constants)
    Δc  = cᵖℓ - cᵖv

    return pᵗʳ * (T / Tᵗʳ)^(Δc / Rᵛ) * exp((ℒ₀ - Δc * T₀) * (1/Tᵗʳ - 1/T) / Rᵛ)
end

@inline function saturation_specific_humidity(T, ρ, thermodynamics_constants)
    p★ = saturation_vapor_pressure(T, thermodynamics_constants)
    Rᵛ = vapor_specific_gas_constant(thermodynamics_constants)
    return p★ / (ρ * Rᵛ * T)
end