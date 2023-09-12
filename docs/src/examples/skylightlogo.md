# Skylight's logo

Skylight's logo is produced with Skylight itself, mapping Julia's logo onto a screen behind a Kerr black hole and ray tracing the image seen by an observer aligned with the black hole and the screen. You can play with the logo by changing the parameters of the configuration in the following snippet, where you need to replace `original_image = load("julia-logo.png")` with the path to your own image file:

```julia

using Skylight
using Images

function main(original_image,
    filename;
    spin = 0.0,
    d_obs_screen,
    rel_bh_pos,
    bh_size,
    xaperture,
    n_factor = 1)
    nx, ny = size(original_image)

    obs_x = d_obs_screen * (1 - rel_bh_pos)
    xscreen = obs_x - d_obs_screen
    # xaperture = 60
    ximage = 1.0 / bh_size

    Nx = floor(Int, nx * n_factor)
    Ny = floor(Int, ny * n_factor)

    spacetime = KerrSpacetimeKerrSchildCoordinates(M = 1.0, a = spin)
    camera = PinholeCamera(position = [0.0, obs_x, 0.0, 0.0],
        horizontal_aperture_in_degrees = xaperture, #rad2deg(70/distance),
        vertical_aperture_in_degrees = (ny / nx) * xaperture, #rad2deg(70/distance),
        horizontal_number_of_pixels = Nx,
        vertical_number_of_pixels = Ny)
    model = VerticalScreen(x = xscreen,
        horizontal_side = ximage,
        vertical_side = (ny / nx) * ximage)
    configurations = VacuumOTEConfigurations(spacetime = spacetime,
        camera = camera,
        radiative_model = model,
        unit_mass_in_solar_masses = 1.0)
    initial_data = initialize(configurations)
    cb, cbp = callback_setup(configurations; rmax = 500.0, rhorizon_bound = 2e-1) #... or, define your own cb and cbp
    run = integrate(initial_data,
        configurations,
        cb,
        cbp;
        method = VCABM(),
        reltol = 1e-13,
        abstol = 1e-21)
    output_data = run.output_data

    # Create a new blank image with the same size as the original image
    new_image = imresize(original_image, (Nx, Ny))

    for j in 1:Ny
        for i in 1:Nx
            ipx = (j - 1) * Nx + i
            x, y, z = output_data[2:4, ipx]
            if !(x ≈ model.x) || abs(y) > 0.5 * model.horizontal_side ||
               abs(z) > 0.5 * model.vertical_side
                new_image[i, j] = zero(eltype(original_image))
                continue
            end

            in = floor(Int,
                1 + nx * (y + 0.5 * model.horizontal_side) / model.horizontal_side)
            jn = floor(Int, 1 + ny * (z + 0.5 * model.vertical_side) / model.vertical_side)

            if 1 <= in <= size(new_image, 1) && 1 <= jn <= size(new_image, 2)
                # Copy the pixel value from the original image to the new image
                new_image[i, j] = original_image[in, jn]
            else
                println("Pixel ($i, $j) is outside the new image, warning")
            end
        end
    end

    # Save the new image
    # display(new_image)
    save(filename, new_image)
end

original_image = load("julia-logo.png")
spin = 0.0
obs_x = 40.0
xscreen = -20.0
ximage = 30.0
xaperture = 60
n_factor = 1

main(original_image,
    "logo.png";
    d_obs_screen = 230.0,
    rel_bh_pos = 3 / 23,
    bh_size = 1 / 15.0,
    xaperture = 10,
    n_factor = 3)
```