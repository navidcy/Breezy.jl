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

struct Saturation{FT}
    energy_reference_temperature :: FT
    triple_point_temperature :: FT
    triple_point_pressure :: FT
end

"""
    Saturation(FT = Oceananigans.defaults.FloatType;
               energy_reference_temperature,
               triple_point_temperature,
               triple_point_pressure)

Returns `Saturation` with specified parameters converted to `FT`.

The `Saturation` struct contains reference values used in the Clausius-Clapeyron relation to calculate
saturation vapor pressure. The Clausius-Clapeyron equation describes how the saturation vapor pressure
changes with temperature for a liquid-vapor phase transition:

The Clausius-Clapeyron relation describes the pressure-temperature relationship during phase transitions:

```math
    d/dT log(p_sat) = L / (Rv * T²)
```
where:

- p_sat is the saturation vapor pressure
- T is temperature 
- L is the latent heat of vaporization
- R_v is the gas constant for water vapor

For water vapor, this integrates to:

    ln(p_sat/p_triple) = (L/R_v) * (1/T_triple - 1/T)

or equivalently:

    p_sat = p_triple * exp((L/R_v) * (1/T_triple - 1/T))

where:
- p_triple is the triple point pressure
- T_triple is the triple point temperature

This relation is used to calculate saturation vapor pressure, which determines the maximum possible
water vapor content of air at a given temperature and pressure.
"""
function Saturation(FT = Oceananigans.defaults.FloatType;
                    energy_reference_temperature = 273.16,
                    triple_point_temperature = 273.16,
                    triple_point_pressure = 611.657)

    return Saturation{FT}(convert(FT, energy_reference_temperature),
                          convert(FT, triple_point_temperature), 
                          convert(FT, triple_point_pressure))
end

struct PhaseTransition{FT}
    latent_heat :: FT
    heat_capacity :: FT
end

"""
    PhaseTransition(FT = Oceananigans.defaults.FloatType; latent_heat, heat_capacity)

Returns `PhaseTransition` with specified parameters converted to `FT`.

Two examples of `PhaseTransition` are condensation and deposition.
During condensation, water molecules in the gas phase cluster together and slow down to form liquid with `heat_capacity`,
The lost of molecular kinetic energy is called the `latent_heat`.

Likewise, during deposition, water molecules in the gas phase cluster into ice crystals.

Arguments
=========
- `FT`: Float type to use (defaults to Oceananigans.defaults.FloatType)
- `latent_heat`: Latent heat lost during phase transition from the gaseous state.
- `heat_capacity`: Heat capacity of the transitional phase.
"""
function PhaseTransition(FT = Oceananigans.defaults.FloatType; latent_heat, heat_capacity)
    return PhaseTransition{FT}(convert(FT, latent_heat),
                               convert(FT, heat_capacity))
end

water_condensation(FT) = PhaseTransition(FT; latent_heat=2500800, heat_capacity=4181)
water_deposition(FT)   = PhaseTransition(FT; latent_heat=2834000, heat_capacity=2108)

struct AtmosphereThermodynamics{FT, S, C, F}
    molar_gas_constant :: FT
    gravitational_acceleration :: FT
    dry_air :: IdealGas{FT}
    vapor :: IdealGas{FT}
    saturation :: S
    condensation :: C
    deposition :: F
end

"""
    AtmosphereThermodynamics(FT = Oceananigans.defaults.FloatType;
                             gravitational_acceleration = 9.81,
                             molar_gas_constant = 8.314462618,
                             dry_air_molar_mass = 0.02897,
                             dry_air_heat_capacity = 1005,
                             vapor_molar_mass = 0.018015,
                             vapor_heat_capacity = 1850,
                             saturation = Saturation(FT),
                             condensation = water_condensation(FT),
                             deposition = water_deposition(FT))

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
                                  saturation = Saturation(FT),
                                  condensation = water_condensation(FT),
                                  deposition = water_deposition(FT))

    dry_air = IdealGas(FT; molar_mass = dry_air_molar_mass,
                           heat_capacity = dry_air_heat_capacity)

    vapor = IdealGas(FT; molar_mass = vapor_molar_mass,
                         heat_capacity = vapor_heat_capacity)

    return AtmosphereThermodynamics(convert(FT, molar_gas_constant),
                                    convert(FT, gravitational_acceleration),
                                    dry_air, vapor, saturation, condensation, deposition)
end

const AT = AtmosphereThermodynamics
const IG = IdealGas

@inline vapor_gas_constant(thermo::AT)    = thermo.molar_gas_constant / thermo.vapor.molar_mass
@inline dry_air_gas_constant(thermo::AT)  = thermo.molar_gas_constant / thermo.dry_air.molar_mass

@inline function mixture_gas_constant(q, thermo::AT)
    Rᵈ = dry_air_gas_constant(thermo)   
    Rᵛ = vapor_gas_constant(thermo)   
    return Rᵈ * (1 - q) + Rᵛ * q
end

"""
    mixture_heat_capacity(q, thermo)

Compute the heat capacity of state air given the total specific humidity q
and assuming that condensate mass ratio qᶜ ≪ q, where qℓ is the mass ratio of
liquid condensate.
"""
@inline function mixture_heat_capacity(q, thermo::AT)
    cᵖᵈ = thermo.dry_air.heat_capacity
    cᵖᵛ = thermo.vapor.heat_capacity
    return cᵖᵈ * (1 - q) + cᵖᵛ * q
end
