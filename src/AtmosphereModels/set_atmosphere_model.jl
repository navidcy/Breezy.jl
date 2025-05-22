import Oceananigans.Fields: set!
using Oceananigans.TimeSteppers: update_state!
using Oceananigans.BoundaryConditions: fill_halo_regions!

function set!(model::AtmosphereModel; kw...)
    for (name, value) in kw

        # Prognostic variables
        if name ∈ propertynames(model.momentum)
            ϕ = getproperty(model.momentum, name)
            set!(ϕ, value)
        elseif name ∈ propertynames(model.tracers)
            ϕ = getproperty(model.tracers, name)
        elseif name == :e
            ϕ = model.energy
        elseif name == :ρq
            ϕ = model.absolute_humidity
        end

        # Setting diagnostic variables
        if name == :θ
            θ = model.temperature # use scratch
            set!(θ, value)

            ρ₀ = model.formulation.reference_density
            cᵖ = model.thermodynamics.dry_air.heat_capacity
            ϕ = model.energy # use scratch
            value = ρ₀ * cᵖ * θ
        elseif name == :q
            q = model.specific_humidity # use scratch
            set!(q, value)

            ρ₀ = model.formulation.reference_density
            ϕ = model.absolute_humidity # use scratch
            value = ρ₀ * q
        end

        set!(ϕ, value)                
        fill_halo_regions!(ϕ, model.clock, fields(model))
    end

    update_state!(model)
end
