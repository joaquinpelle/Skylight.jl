The energy momentum tensor of the radiation field in terms of the specific intensity is given by

$$
T^{a b} = \int_{S^2} \int_0^\infty n^a n^b I_\nu d\nu d\Omega 
$$

where $n^a = k^a / \nu$, and the integral runs over all the frequency spectrum and directions. The radiative flux observed by an osberver at four-velocity $u^a$ through a surface element with normal $\bar{n}^a$ is $T^{a b} u_a \bar{n}_b$.

For discretizing this integral at a position $x^\mu$, we first construct an orthonormal tetrad. We choose the timelike vector $u^\mu$ by normalizing $\partial_t$. For the first spacelike vector we take [0, -x, -y, -z] (thus "pointing towards the origin", which is precise when the spacetime is flat), project it orthogonally to $u^\mu$ and normalize it. For the rest of the tetrad, we choose two other vectors randomly and orthonormalize them with respect to the first two.  

Finally, we take coordinates $(\alpha, \beta)$ on $S^2$ such that $\alpha = -\varphi$ and $\beta = \pi/2 -\theta$ where $\theta$ and $\varphi$ are usual spherical coordinates. In this manner, the coordinates are nicely centered in the direction of the first spatial vector (which "points towards the center"), and we only need to take the $(\alpha, \beta)$ coordinate ranges large enough to cover the image of the emitting source. The volume element in these coordinates reads $d\Omega = \cos \beta d\alpha d\beta$. 

Notice that the choice of the tetrad and angular coordinates doesn't really matter. The procedure we choose is sufficiently general to allow an arbitrary spacetime position, surface element, and four-velocity, but particularly well suited for large distances when the spacetime is approximately flat and the source can be covered by a small spherical sector nicely parameterized with $(\alpha, \beta)$ coordinates.

Finally, we discretize discretize the integral as

$$
\sum_{i j k} I_{ijk} n^a_{ij} n^b_{ij} \Delta \Omega_{ij} \Delta \nu \, ,    
$$

where $I_{ijk} = I_{\nu_k}(n^a_{ij})$, the tetrad components of $n^a_{ij}$ are 
$$
[1, \cos \alpha_i \cos \beta_j, -\sin \alpha_i \cos \beta_j, \sin \beta_j]\, ,
$$
the solid angle of each section is
$$
\Delta \Omega_{ij} = \int_{D_{ij}} \cos \beta d\alpha d\beta = [\sin(\beta_{j+1}) - \sin(\beta_{j})] \Delta \alpha \Delta \beta,
$$
and $D_{ij} = [\alpha_i, \alpha_{i+1}] \times [\beta_j, \beta_{j+1}]$. We took a uniform grid such that 
$$
\alpha_i = -\alpha_{max}+(i-1)\Delta \alpha \\ 
\beta_j = -\beta_{max}+(j-1)\Delta \beta \\
\Delta \alpha = 2\alpha_{max}/(N_\alpha-1) \\
\Delta \beta = 2\beta_{max}/(N_\beta-1) \\
1 \le i \le N_\alpha \\ 
1 \le j \le N_\beta
$$
  