The energy momentum tensor of the radiation field in terms of the specific intensity is given by

$$
T{ab} = \int k^a k^b \left( \frac{I_\nu}{\nu^3} \right) dV_k 
$$

where $dV_k = \nu d\nu d\Omega$. Equivalently,

$$
T^{a b} = \int_{S^2} \int_0^\infty n^a n^b I_\nu d\nu d\Omega 
$$

where $n^a = k^a / \nu$. The radiative flux through a surface element with normal $\bar{n}^a$ measured by an osberver with four-velocity $u^a$ is $T^{a b} u_a \bar{n}_b$.

For discretizing this integral at position $x^\mu$, we choose an orthonormal tetrad. We take
$$
e_0 = \partial_t \\
e_1 = -\partial_r \\
e_2 = \partial_\phi \\
e_3 = -\partial_\theta \\
$$
and orthonormalize it, where $(r,\theta, \phi)$ are some topolgically spherical coordinates on the spacetime.
Note that this assumes there is a timelike coordinate $t$ and three spacelike coordinates, and that there exist some meaningful way of transforming to spherical-like coordinates. In flat spacetime, the interpretation of this tetrad is simple: $e_1$ points towards the origin, and $e_2, e_3$ are parallel to the equatorial plane and the $z$-axis respectively. Note that the triad we use is not direct. Note also that the choice of tetrad does not actually matter. The procedure we choose is sufficiently general to allow (almost) arbitrary spacetime position, surface element, and four-velocity, but it is particularly well suited for large distances in asymptotically flat spacetimes where the source can be covered by a small spherical sector centered around $e_1$.
Finally, we take coordinates $(\alpha, \beta)$ on $S^2$ such that 
$$
\alpha = \varphi \\
\beta = \pi/2 -{\vartheta}
$$
where $(\vartheta, \varphi)$ are the usual angular coordinates on $S^2$. The coordinates range over $-π \le \alpha \lt \pi$ and $-\pi/2 \le \beta \le \pi/2$. The tetrad components of a vector in these coordinates can be written as
$$
k^0 = \nu \\ 
k^1 = \nu \cos \alpha \cos \beta \\
k^2 = \nu \sin \alpha \cos \beta \\
k^3 = \nu \sin \beta 
$$ 

Thus, in particular $(\alpha, \beta) = (0,0)$ map to $e_1$. In these coordinates, $d\Omega = \cos \beta d\alpha d\beta$. We only need to take the $(\alpha, \beta)$ coordinate ranges large enough to cover the image of the emitting source. 

 
Finally, we discretize discretize the integral as

$$
\sum_{i j k} I_{ijk} n^a_{ij} n^b_{ij} \Delta \Omega_{ij} \Delta \nu \, ,    
$$

where $I_{ijk} = I_{\nu_k}(n^a_{ij})$, the tetrad components of $n^a_{ij}$ are 
$$
[1, \cos \alpha_i \cos \beta_j, \sin \alpha_i \cos \beta_j, \sin \beta_j]\, ,
$$
(we have to transform them to the coordinate frame before contraction), the solid angle of each section is
$$
\Delta \Omega_{ij} = \int_{D_{ij}} \cos \beta d\alpha d\beta = 2 \cos(\beta_j) \sin \left(\frac{\Delta \beta}{2} \right) \Delta \alpha,
$$
and $D_{ij} = [\alpha_i-\Delta \alpha/2, \alpha_i+\Delta \alpha/2] \times [\beta_j - \Delta \beta / 2, \beta_j + \Delta \beta /2]$. We took a uniform grid such that 
$$
\alpha_i = -s_\alpha/2+(i-1/2)\Delta \alpha \\ 
\beta_j = -s_\beta/2+(j-1/2)\Delta \beta \\
\Delta \alpha = s_\alpha/N_\alpha \\
\Delta \beta = s_\beta/N_\beta \\
1 \le i \le N_\alpha \\ 
1 \le j \le N_\beta
$$. 

where $s_\alpha, s_\beta$ are the horizontal and vertical aperture angles respectively, and $N_\alpha$, $N_\beta$ are the numbers of pixels. The total solid angle is $\Delta \Omega = 2 \sin(s_\beta/2)s_\alpha$
  