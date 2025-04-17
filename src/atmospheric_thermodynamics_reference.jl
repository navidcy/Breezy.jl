struct ConstitutiveConstants{FT}
    gas_constant       :: FT
    dry_air_molar_mass :: FT
    water_molar_mass   :: FT
end

function Base.summary(p::ConstitutiveConstants{FT}) where FT
    return string("ConstitutiveConstants{$FT}(",
                    "R=", prettysummary(p.gas_constant),
                  ", Mᵈ=", prettysummary(p.dry_air_molar_mass),
                  ", Mᵛ=", prettysummary(p.water_molar_mass), ")")
end

Base.show(io::IO, p::ConstitutiveConstants) = print(io, summary(p))

"""
    ConstitutiveConstants(FT = Float64;
                          gas_constant       = 8.3144598,
                          dry_air_molar_mass = 0.02897,
                          water_molar_mass   = 0.018015)

Construct a set of thermodynamic constants define the density of moist air,

```math
ρ = p / Rᵐ(q) T,
```

where ``p`` is pressure, ``T`` is temperature, ``q`` defines the partition
of total mass into vapor, liqiud, and ice mass fractions, and
``Rᵐ`` is the effective specific gas constant for the mixture.

For more information see [reference docs].
"""
function ConstitutiveConstants(FT = Oceananigans.defaults.FloatType;
                               gas_constant       = 8.3144598,
                               dry_air_molar_mass = 0.02897,
                               water_molar_mass   = 0.018015)

    return ConstitutiveConstants{FT}(convert(FT, gas_constant),
                                     convert(FT, dry_air_molar_mass),
                                     convert(FT, water_molar_mass))
end

const CP{FT} = ConstitutiveConstants{FT} where FT

struct HeatCapacityConstants{FT}
    dry_air_adiabatic_exponent :: FT
    water_vapor_heat_capacity  :: FT
    liquid_water_heat_capacity :: FT
    water_ice_heat_capacity    :: FT
end

function Base.summary(p::HeatCapacityConstants{FT}) where FT
    return string("HeatCapacityConstants{$FT}(",
                    "κᵈ=", prettysummary(p.dry_air_adiabatic_exponent),
                  ", cᵖᵛ=", prettysummary(p.water_vapor_heat_capacity),
                  ", cᵖˡ=", prettysummary(p.liquid_water_heat_capacity),
                  ", cᵖⁱ=", prettysummary(p.water_ice_heat_capacity))
end

Base.show(io::IO, p::HeatCapacityConstants) = print(io, summary(p))

"""
    HeatCapacityConstants(FT = Float64,
                          dry_air_adiabatic_exponent = 2/7,
                          water_vapor_heat_capacity = 1859,
                          liquid_water_heat_capacity = 4181,
                          water_ice_heat_capacity = 2100)

Isobaric heat capacities.
"""
function HeatCapacityConstants(FT = Oceananigans.defaults.FloatType;
                               dry_air_adiabatic_exponent = 2/7,
                               water_vapor_heat_capacity = 1859,
                               liquid_water_heat_capacity = 4181,
                               water_ice_heat_capacity = 2100)

    return HeatCapacityConstants{FT}(convert(FT, dry_air_adiabatic_exponent),
                                     convert(FT, water_vapor_heat_capacity),
                                     convert(FT, liquid_water_heat_capacity),
                                     convert(FT, water_ice_heat_capacity))
end

const HCP{FT} = HeatCapacityConstants{FT} where FT

struct PhaseTransitionConstants{FT}
    reference_vaporization_enthalpy  :: FT
    reference_sublimation_enthalpy   :: FT
    reference_temperature            :: FT
    triple_point_temperature         :: FT
    triple_point_pressure            :: FT
    water_freezing_temperature       :: FT
    total_ice_nucleation_temperature :: FT
end

function Base.summary(p::PhaseTransitionConstants{FT}) where FT
    return string("PhaseTransitionConstants{$FT}(",
                    "ℒᵛ⁰=", prettysummary(p.reference_vaporization_enthalpy),
                  ", ℒˢ⁰=", prettysummary(p.reference_sublimation_enthalpy),
                  ", T⁰=", prettysummary(p.reference_temperature),
                  ", Tᵗʳ=", prettysummary(p.triple_point_temperature),
                  ", pᵗʳ=", prettysummary(p.triple_point_pressure),
                  ", Tᶠ=", prettysummary(p.water_freezing_temperature),
                  ", Tⁱⁿ=", prettysummary(p.total_ice_nucleation_temperature), ')')
end

Base.show(io::IO, p::PhaseTransitionConstants) = print(io, summary(p))

function PhaseTransitionConstants(FT = Oceananigans.defaults.FloatType;
                                  reference_vaporization_enthalpy = 2500800,
                                  reference_sublimation_enthalpy = 2834400,
                                  reference_temperature = 273.16,
                                  triple_point_temperature = 273.16,
                                  triple_point_pressure = 611.657,
                                  water_freezing_temperature = 273.15,
                                  total_ice_nucleation_temperature = 233)

    return PhaseTransitionConstants{FT}(convert(FT, reference_vaporization_enthalpy),
                                        convert(FT, reference_sublimation_enthalpy),
                                        convert(FT, reference_temperature),
                                        convert(FT, triple_point_temperature),
                                        convert(FT, triple_point_pressure),
                                        convert(FT, water_freezing_temperature),
                                        convert(FT, total_ice_nucleation_temperature))
end

const PTP{FT} = PhaseTransitionConstants{FT} where FT

struct MoistAirThermodynamicsConstants{FT}
    constitutive      :: ConstitutiveConstants{FT}
    heat_capacity     :: HeatCapacityConstants{FT}
    phase_transitions :: PhaseTransitionConstants{FT}
end

const PATP{FT} = MoistAirThermodynamicsConstants{FT} where FT

Base.eltype(::PATP{FT}) where FT = FT
Base.eltype(::CP{FT})   where FT = FT
Base.eltype(::HCP{FT})  where FT = FT
Base.eltype(::PTP{FT})  where FT = FT

Base.summary(::PATP{FT}) where FT = "MoistAirThermodynamicsConstants{$FT}"

function Base.show(io::IO, p::MoistAirThermodynamicsConstants)
    FT = eltype(p)

    cp = p.constitutive 
    hc = p.heat_capacity
    pt = p.phase_transitions

    return print(io, summary(p), ':', '\n',
        "├── ConstitutiveConstants{$FT}:", '\n',
        "│   ├── gas_constant (R):                      ", prettysummary(cp.gas_constant), '\n',
        "│   ├── dry_air_molar_mass (Mᵈ):               ", prettysummary(cp.dry_air_molar_mass), '\n',
        "│   └── water_molar_mass (Mᵛ):                 ", prettysummary(cp.water_molar_mass), '\n',
        "├── HeatCapacityConstants{$FT}:", '\n',
        "│   ├── dry_air_adiabatic_exponent (κᵈ):       ", prettysummary(hc.dry_air_adiabatic_exponent), '\n',
        "│   ├── water_vapor_heat_capacity (cᵖᵛ):       ", prettysummary(hc.water_vapor_heat_capacity), '\n',
        "│   ├── liquid_water_heat_capacity (cᵖˡ):      ", prettysummary(hc.liquid_water_heat_capacity), '\n',
        "│   └── water_ice_heat_capacity (cᵖⁱ):         ", prettysummary(hc.water_ice_heat_capacity), '\n',
        "└── PhaseTransitionConstants{$FT}", '\n',
        "    ├── reference_vaporization_enthalpy (ℒᵛ⁰): ", prettysummary(pt.reference_vaporization_enthalpy), '\n',
        "    ├── reference_sublimation_enthalpy  (ℒˢ⁰): ", prettysummary(pt.reference_sublimation_enthalpy), '\n',
        "    ├── reference_temperature (T⁰):            ", prettysummary(pt.reference_temperature), '\n',    
        "    ├── triple_point_temperature (Tᵗʳ):        ", prettysummary(pt.triple_point_temperature), '\n',
        "    ├── triple_point_pressure (pᵗʳ):           ", prettysummary(pt.triple_point_pressure), '\n',   
        "    ├── water_freezing_temperature (Tᶠ):       ", prettysummary(pt.water_freezing_temperature), '\n',
        "    └── total_ice_nucleation_temperature (Tⁱ): ", prettysummary(pt.total_ice_nucleation_temperature))
end

function MoistAirThermodynamicsConstants(FT = Oceananigans.defaults.FloatType;
                                         constitutive = ConstitutiveConstants(FT),
                                         phase_transitions = PhaseTransitionConstants(FT),
                                         heat_capacity = HeatCapacityConstants(FT))

    return MoistAirThermodynamicsConstants(constitutive, heat_capacity, phase_transitions)
end

const PATP = MoistAirThermodynamicsConstants
