# Camera

In the observer-to-emitter method, we integrate the radiative transfer equation backwards in time starting from the observation point. For this, a bundle of past-directed rays, forming what is usually called a virtual camera, is traced from the observer to the source, . See [Pinhole camera](@ref) for a detailed explanation of the camera construction and flux calculation methods.

## Examples 

In Skylight.jl, the camera is represented by the [`PinholeCamera`](@ref) type containing the parameters that define the camera position, the field of view, and the resolution. The camera can be constructed as follows:

```julia
camera = PinholeCamera(position = [0.0, 500, π/2-π/20, 0.0],
                        horizontal_aperture_in_degrees = 4,
                        vertical_aperture_in_degrees = 4,
                        horizontal_number_of_pixels = 600,
                        vertical_number_of_pixels = 600)
```

The camera position is given by the spacetime coordinates of the observation point, and must be in the same coordinates as those of the spacetime in which it is set. The apertures define the field of view, and ideally should be set to capture a full view of the emitting source. Optionally, a four-velocity can be set for the observation frame, with components given in the coordinate basis: 

```julia
camera = PinholeCamera(position = [0.0, 500, π/2-π/20, 0.0],
                        horizontal_aperture_in_degrees = rad2deg(80/500),
                        vertical_aperture_in_degrees = rad2deg(80/500),
                        horizontal_number_of_pixels = 600,
                        vertical_number_of_pixels = 600,
                        four_velocity = [1.0, -0.5, 0.0, 0.0])
```

Otherwise, the observation frame will be assumed to be the static frame.

Most ray-tracing applications use the so-called image plane instead of a pinhole camera. This is an approximation in which the spacetime is assumed to be asymptotically flat, and the observer is at a large distance from the source, such that the light rays arrive mutually parallel. Such a camera construction is available in Skylight.jl as the [`ImagePlane`](@ref) type, and can be constructed as follows:

```julia
camera = ImagePlane(distance = 500.0,
    observer_inclination_in_degrees = 90,
    horizontal_side = 23.0,
    vertical_side = 23.0,
    horizontal_number_of_pixels = 100,
    vertical_number_of_pixels = 100,
    observation_times = [0.0, 100.0])
```

In this case, the observation frame is assumed to be static, thus allowing to consider multiple observation times self-consistently, with light rays starting from the same position. For more details on this setup, see the [Image plane](@ref) section. However, when using Skylight.jl, we advise to always use the pinhole camera as it gives more accurate results and the computational cost is not much higher than for the image plane.

## Types

```@docs
PinholeCamera
ImagePlane
```

## Details of the construction 

### Pinhole camera

The energy-momentum tensor of the radiation field reads

$T^{ab} = \int k^a k^b \left( \frac{I_\nu}{\nu^3} \right) dK,$

where $dK = \nu d\nu d\Omega$ is the invariant volume element in the light-like momentum space. This can also be written as

$T^{ab} = \int_{S^2} \int_0^\infty n^a n^b I_\nu d\nu d\Omega,$

where $n^a = k^a / \nu$. The radiative flux measured in the direction of $\bar{n}^a$ in a reference frame with four-velocity $u^a$ is defined as $T^{ab} u_a \bar{n}_b$.

To discretize this integral at position $x^\mu$, we take an orthonormal tetrad given by

$e_0 = \partial_t, \quad e_1 = -\partial_r, \quad e_2 = \partial_\phi, \quad e_3 = -\partial_\theta,$

where $(r, \theta, \phi)$ are some topologically spherical coordinates in spacetime. Here we are assuming the spacetime coordinates can be transformed to a topologically spherical system consisting of a temporal coordinate $t$ and three spatial coordinates. This might seem restrictive in theory, but it is generic in practice. After calculating the set of vectos above, we orthonormalize it with respect to $g_{\alpha \beta}(x^\mu)$, the metric at that position. In flat spacetime with usual spherical coordinates, the interpretation of the resulting tetrad is simple: $e_1$ points towards the origin, $e_2$ is azimuthal, and $e_3$ is polar. The construction is somewhat arbitrary, but in the end the choice of tetrad is not important; it is merely an artifact for parameterizing the flux integral and defining an integration domain as small as possible. This construction is general enough to allow arbitrary position and four-velocity in any spacetime, but it is particularly convenient for large distances in asymptotically flat spacetimes where the source can be covered by a small spherical sector centered in the radial direction.

To calculate the flux integral, we take angular coordinates $(\alpha, \beta)$ in $S^2$ such that the components in the tetrad of a vector can be written as

$k^0 = \nu\,, \quad k^1 = \nu \cos \alpha \cos \beta\,, \quad k^2 = \nu \sin \alpha \cos \beta\,, \quad k^3 = \nu \sin \beta\,.$

The coordinates vary over $-\pi \le \alpha < \pi$ and $-\pi/2 \le \beta \le \pi/2$. In particular, $(\alpha, \beta) = (0,0)$ maps to the direction of $e_1$. Additionally, in these coordinates, $d\Omega = \cos \beta d\alpha d\beta$. In practice, for the angular integration, we do not need to cover the entire celestial sphere from the observation point, but rather take ranges of $(\alpha, \beta)$ that cover the image of the emitting source.

We then take a uniform angular grid given by

$\alpha_i = -\frac{s_\alpha}{2} + \left(i - \frac{1}{2}\right) \Delta \alpha, \quad \beta_j = -\frac{s_\beta}{2} + \left(j - \frac{1}{2}\right) \Delta \beta, \\
\Delta \alpha = \frac{s_\alpha}{N_\alpha}, \quad \Delta \beta = \frac{s_\beta}{N_\beta}, \quad 1 \le i \le N_\alpha, \quad 1 \le j \le N_\beta,$

where $s_\alpha, s_\beta$ are the horizontal and vertical angular apertures, respectively, and $N_\alpha, N_\beta$ are the numbers of points per side of the grid. We take a partition into spherical sectors $D_{ij} = [\alpha_i - \Delta \alpha/2, \alpha_i + \Delta \alpha/2] \times [\beta_j - \Delta \beta / 2, \beta_j + \Delta \beta /2]$ centered at $(\alpha_i, \beta_j)$, whose solid angle is

$\Delta \Omega_{ij} = \int_{D_{ij}} \cos \beta d\alpha d\beta = 2 \cos(\beta_j) \sin \left(\frac{\Delta \beta}{2} \right) \Delta \alpha,$

Finally, we numerically approximate the flux integral as

$\sum_{ijk} I_{ijk} n^a_{ij} n^b_{ij} \Delta \Omega_{ij} \Delta \nu,$

where $I_{ijk} = I_{\nu_k}(n^a_{ij})$ and the components of the vectors $n^a_{ij}$ in the tetrad are

$(1, \cos \alpha_i \cos \beta_j, \sin \alpha_i \cos \beta_j, \sin \beta_j),$

which are transformed to the coordinate frame before contraction with $u_a$ and $\bar{n}_a$.

### Image plane

For distant observation points in asymptotically flat spacetimes, it is reasonable to approximate that the light rays from the source arrive essentially parallel. This is the approach that many ray-tracing applications adopt, including in the [original Skylight paper](https://academic.oup.com/mnras/article-abstract/515/1/1316/6631564?login=false). Although this approximation is generally sufficient for most practical applications, it is not necessarily less computationally demanding than the more accurate pinhole camera approach.

Consider a Cartesian coordinate system $(t,x,y,z)$ such that the observation point is located on the $x$-$z$ plane at a distance $D$ from the origin, and the inclination with respect to the $z$-axis is $\xi$. In the static frame, $u^a = (\partial_t)^a$, we have $u^a n_a = 1$, where $n^a$ is the unit vector in the propagation direction of a photon. Since the observation point is distant, the image of the emitting source occupies some small spherical sector. If we take $\bar{n}^a$ pointing to the center of such sector, then $\bar{n}^a n_a \approx 1$ within the source's image. Given that the rays are approximately parallel, we can approximate the flux as

$F_\nu = \frac{1}{D^2} \int_{\mathcal{S}} I_\nu da db \,,$

where $\mathcal{S}$ is the image plane and $(a,b)$ are rectangular coordinates on $\mathcal{S}$. These coordinates are related to the Cartesian coordinates through

$x = -b \cos \xi + D \sin \xi\,, \quad y = a\,, \quad z = b \sin \xi + D \cos \xi \,.$

To calculate this integral, we approximate it numerically as a Riemann sum from a rectangular grid on the image plane, similarly to the method described in the [Pinhole camera](@ref) section. Each grid point is taken as the initial position of a geodesic, with the initial spatial momentum normal to the image plane, and the temporal component of the momentum set such that the resulting four-vector is null. The specific intensity for each ray is also calculated in the same manner as described previously.
