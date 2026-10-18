# Run from the repository root; see docs/src/examples/torus_kerr.md.
# Including this file defines helpers without launching a ray-tracing run.
module IonTorusExample

using Skylight
using HDF5
using DelimitedFiles

const CASES = (
    (id="C1", a=0.1, R0=0.0),
    (id="C2", a=0.1, R0=-0.001),
    (id="C3", a=0.99, R0=0.0),
    (id="C4", a=0.99, R0=-0.001),
)
const FREQUENCIES_HZ = [86e9, 230e9, 345e9]

function configuration(case; pixels=256)
    pixels >= 4 && pixels % 4 == 0 ||
        throw(ArgumentError("Horizontal pixels must be a positive multiple of four, at least four."))
    spacetime = FRKerrSpacetime(M=1.0, a=case.a, R0=case.R0)
    model = IonTorus(spacetime; λ=0.7, ϵc=1e-17, n=1.5, β=0.5, ξ=0.01,
        H_abundance=0.75, He_abundance=0.25, thermodynamics=:hydrostatic,
        radiative_process=SynchrotronAndBremsstrahlung(), self_absorption=true)
    camera = PinholeCamera(position=[0.0, 50.0, π/3, 0.0],
        horizontal_aperture_in_degrees=120.0, vertical_aperture_in_degrees=90.0,
        horizontal_number_of_pixels=pixels, vertical_number_of_pixels=3pixels÷4)
    return NonVacuumOTEConfigurations(spacetime=spacetime, camera=camera,
        radiative_model=model, unit_mass_in_solar_masses=4.14e6,
        observation_energies=PhysicalConstants.h .* FREQUENCIES_HZ)
end

function validate_domain(config; horizon_buffer=1e-3)
    model, spacetime, camera = config.radiative_model, config.spacetime, config.camera
    isfinite(horizon_buffer) && horizon_buffer>0 ||
        throw(ArgumentError("The horizon buffer must be finite and positive."))
    event_horizon_radius(spacetime)+horizon_buffer < model.rcusp ||
        throw(ArgumentError("The horizon cutoff would remove part of the torus."))
    model.router < camera.position[2] ||
        throw(ArgumentError("The camera must be outside the emitting torus."))
    return nothing
end

# Six-decade logarithmic preview, separate from unmodified scientific data.
function write_preview(path, intensity, case_id)
    nx,ny=size(intensity)
    peak=maximum(intensity)
    open(path,"w") do io
        println(io,"<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 $nx $(ny+12)\">")
        println(io,"<rect width=\"100%\" height=\"100%\" fill=\"black\"/>")
        for j in 1:ny, i in 1:nx
            value = peak>0 && intensity[i,j]>0 ? clamp(1+log10(intensity[i,j]/peak)/6,0,1) : 0.0
            gray=round(Int,255value)
            gray==0 && continue
            println(io,"<rect x=\"$(i-1)\" y=\"$(ny-j)\" width=\"1\" height=\"1\" fill=\"rgb($gray,$gray,$gray)\"/>")
        end
        println(io,"<text x=\"2\" y=\"$(ny+8)\" font-size=\"5\" fill=\"white\">$case_id: 230 GHz, log stretch (6 decades)</text></svg>")
    end
end

function run_case(case; pixels=256, output="iontorus_results", horizon_buffer=1e-3,
    dtmax=0.25, reltol=1e-8, abstol=1e-10, intensity_abstol=1e-18, dry_run=false)
    config=configuration(case; pixels)
    validate_domain(config; horizon_buffer)
    model=config.radiative_model
    println("$(case.id): a=$(case.a), R0=$(case.R0), ell=$(model.l0), ",
        "cusp=$(model.rcusp), center=$(model.rcenter), outer=$(model.router), Tec=$(model.Tec) K")
    dry_run && return config
    all(x -> isfinite(x) && x>0,(dtmax,reltol,abstol,intensity_abstol)) ||
        throw(ArgumentError("Step limit and solver tolerances must be finite and positive."))
    mkpath(output)
    base=joinpath(output,case.id)
    any(isfile, (base*".h5",base*"_spectrum.csv",base*"_230GHz.svg")) &&
        throw(ArgumentError("Outputs for $(case.id) already exist. Choose a new --output directory."))
    initial=initialize(config)
    cb,cbp=callback_setup(config; rhorizon_bound=horizon_buffer)
    # J=I_E/E^3 is enormous in CGS at radio energies. A scalar absolute
    # tolerance would demand meaningless precision on initially dark rays.
    # Specify the intensity tolerance in observable I_nu units instead.
    ode_abstol=vcat(fill(abstol,11),
        intensity_abstol ./ (PhysicalConstants.h .* config.observation_energies.^3))
    result=integrate(initial,config,cb,cbp; method=Tsit5(), τmax=Inf,
        dtmax, reltol, abstol=ode_abstol, maxiters=10^6, save_everystep=false)
    final=result.output_data
    intensity_E=observed_specific_intensities(initial,final,config)
    intensity_nu=PhysicalConstants.h .* intensity_E
    flux_nu=PhysicalConstants.h .* spectrum(initial,final,config)
    all(isfinite,intensity_nu) && all(>=(0),intensity_nu) &&
        all(isfinite,flux_nu) && all(>=(0),flux_nu) ||
        error("Non-finite or negative observable for $(case.id).")
    alpha,beta=axes_ranges(config.camera)
    nx,ny=length(alpha),length(beta)
    # Arrays use (frequency, horizontal angle, vertical angle).
    h5open(base*".h5","w") do file
        file["frequency_Hz"]=FREQUENCIES_HZ
        file["alpha_rad"]=collect(alpha)
        file["beta_rad"]=collect(beta)
        file["I_nu"]=reshape(intensity_nu,3,nx,ny)
        file["I_E"]=reshape(intensity_E,3,nx,ny)
        file["tau"]=reshape(final[9:11,:],3,nx,ny)
        file["F_nu"]=flux_nu
        file["initial_data"]=initial
        file["final_data"]=final
        file["camera_position"]=config.camera.position
        file["camera_tetrad"]=Skylight.default_tetrad(config.camera,config.spacetime)
        meta=attributes(file)
        for (key,value) in pairs((case_id=case.id, M=1.0, a=case.a, R0=case.R0,
            mass_unit_solar=4.14e6, lambda=model.λ, ell=model.l0,
            rho_c=model.ϵc, n=model.n, beta=model.β, xi=model.ξ,
            X=model.H_abundance, Y=model.He_abundance, Tec=model.Tec,
            cusp=model.rcusp, center=model.rcenter, outer=model.router,
            horizon_buffer=horizon_buffer, dtmax=dtmax, reltol=reltol, abstol=abstol,
            intensity_abstol=intensity_abstol,
            tau_max=Inf, thermodynamics="hydrostatic", self_absorption=1,
            emission="thermal synchrotron + bremsstrahlung; no Compton",
            synchrotron_cutoff_K=5e8, observer="static, r=50 M, inclination=60 degrees",
            I_nu_units="erg s^-1 cm^-2 sr^-1 Hz^-1",
            I_E_units="erg s^-1 cm^-2 sr^-1 erg^-1",
            F_nu_units="erg s^-1 cm^-2 Hz^-1, local camera flux"))
            meta[string(key)]=value
        end
    end
    open(base*"_spectrum.csv","w") do io
        println(io,"frequency_Hz,F_nu_erg_s-1_cm-2_Hz-1,nu_F_nu_erg_s-1_cm-2")
        writedlm(io,hcat(FREQUENCIES_HZ,flux_nu,FREQUENCIES_HZ.*flux_nu),',')
    end
    write_preview(base*"_230GHz.svg",reshape(intensity_nu[2,:],nx,ny),case.id)
    println("Saved $(case.id) to $(abspath(output)); max tau = $(maximum(final[9:11,:]))")
    return (; config, initial, final, intensity_nu, flux_nu)
end

function main(args=ARGS)
    options=Dict("case"=>"all","pixels"=>"256","output"=>"iontorus_results",
        "horizon-buffer"=>"0.001","dtmax"=>"0.25","reltol"=>"1e-8","abstol"=>"1e-10",
        "intensity-abstol"=>"1e-18")
    dry_run=false
    i=1
    while i<=length(args)
        arg=args[i]
        if arg=="--help"
            println("Usage: julia --project=. --threads=auto run/iontorus.jl [options]\n",
                "  --case all|C1|C2|C3|C4  --pixels 256  --output iontorus_results\n",
                "  --horizon-buffer 0.001  --dtmax 0.25  --reltol 1e-8  --abstol 1e-10\n",
                "  --intensity-abstol 1e-18  (CGS I_nu absolute tolerance)\n",
                "  --dry-run  (validate and print parameters without tracing)")
            return
        elseif arg=="--dry-run"
            dry_run=true
        elseif startswith(arg,"--") && haskey(options,arg[3:end]) && i<length(args)
            i+=1
            options[arg[3:end]]=args[i]
        else
            throw(ArgumentError("Unknown option or missing value: $arg"))
        end
        i+=1
    end
    selected=filter(c -> options["case"]=="all" || c.id==options["case"],CASES)
    isempty(selected) && throw(ArgumentError("Case must be all, C1, C2, C3 or C4."))
    for case in selected
        run_case(case; pixels=parse(Int,options["pixels"]), output=options["output"], dry_run,
            horizon_buffer=parse(Float64,options["horizon-buffer"]),
            dtmax=parse(Float64,options["dtmax"]), reltol=parse(Float64,options["reltol"]),
            abstol=parse(Float64,options["abstol"]),
            intensity_abstol=parse(Float64,options["intensity-abstol"]))
    end
end

end # module

if abspath(PROGRAM_FILE)==@__FILE__
    IonTorusExample.main()
end
