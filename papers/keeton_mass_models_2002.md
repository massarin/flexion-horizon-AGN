Version 4 June 2002

A Catalog of Mass Models for Gravitational Lensing

Charles R. Keeton1

Astronomy and Astrophysics Department, University of Chicago,
5640 S. Ellis Ave., Chicago, IL 60637

ABSTRACT

Many diﬀerent families of mass models are used in modern applications of strong
gravitational lensing.
I review a wide range of popular models, with two points of
emphasis: (1) a discussion of strategies for building models suited to a particular lensing
problem; and (2) a summary of technical results for a canonical set of models. All of the
models reviewed here are included in publicly-available lensing software called gravlens.

1.

Introduction

As applications of strong gravitational lensing have become more sophisticated, the variety
and complexity of mass distributions used for lensing studies have increased. Gone are the days
when the singular isothermal sphere was all you needed to know. Now we have softened power law
ellipsoids, pseudo-Jaﬀe models, generalized NFW models, exponential disks embedded in isothermal
halos, and so on. Plus, all of these models are likely to be perturbed by other galaxies, groups, or
clusters near the lens galaxy or along the line of sight.

I have developed software called gravlens for a variety of lensing and lens modeling applica-
tions.2 In the course of writing the code I have collected most of the mass models used for lensing
studies, and this paper reviews those models. The outline is as follows. Section 2 discusses circular
and elliptical symmetry in lensing mass distributions, and argues that a canonical set of circular
and elliptical models provides a useful basis set for building much more complex composite mod-
els. Section 3 gives some suggestions (propaganda, really) about selecting models appropriate for
your application. Section 4 presents the general equations that describe the lensing properties of
a given mass distribution. Finally, Section 5 presents a detailed catalog of results for a number of
models. The mass models discussed in this paper (and included in the gravlens software) are listed
in Table 1.

1Hubble Fellow

2The software is discussed in a separate paper (Keeton 2001), and is available to the community via the web site

of the CfA/Arizona Space Telescope Lens Survey, at http://cfa-www.harvard.edu/castles.

– 2 –

2. Circular, Elliptical, and Composite Lens Models

The lensing properties of any mass distribution can be written in terms of two-dimensional
integrals over the surface mass density (see §4).
In general the integrals cannot be evaluated
analytically, but many lensing applications oﬀer simpliﬁcations due to symmetry. Applications
such as microlensing may permit the use of a lens with circular symmetry, in which case the lensing
properties can usually be found analytically. In other applications, such as lensing by galaxies, it
may be reasonable to assume elliptical symmetry, which allows the lensing properties to be written
as a set of one-dimensional integrals (see §4); the integrals can sometimes be evaluated analytically
and are always amenable to fast numerical techniques.

The geometric simpliﬁcations may not hold in increasingly sophisticated lensing applications.
For example, companion stars or planets can break circular symmetry in microlensing (e.g., Mao
& Paczy´nski 1991; Gould & Loeb 1992); neighboring galaxies (e.g., Young et al. 1981; Schechter et
al. 1997; Koopmans & Fassnacht 1999) or the lens galaxy’s internal structure (e.g., Maller, Flores
& Primack 1997; Bernstein & Fischer 1999; Keeton et al. 2000) may break elliptical symmetry in
lensing by galaxies; and individual cluster members can ruin the symmetry of the smooth back-
ground potential in lensing by clusters (e.g., Bartelmann & Steinmetz 1995; Tyson, Kochanski &
Dell’Antonio 1998). Nevertheless, in all of these examples the total mass distribution can be written
as a combination of circular or elliptical components, either placed at diﬀerent positions to repre-
sent diﬀerent objects, or combined at the same position to mimic an object with more complicated
internal structure (also see Schneider & Weiss 1991). Because the Poisson equation is linear, the
lensing potential of such composite models is simply the sum of the component potentials. In other
words, arbitrary combinations of circular and elliptical models yield a wide range of complex mass
distributions whose lensing properties are nevertheless easy to compute.

While composite models provide a great deal of freedom and complexity, they are not com-
pletely general. This limitation is eliminated in the elegant algorithm by Saha & Williams (1997;
also Williams & Saha 2000) for ﬁnding non-parametric lens models. The approach is to introduce a
set of mass pixels and construct a large linear programming problem for determining their masses.
The problem is severely underconstrained, but by requiring positive-deﬁnite masses and imposing
some smoothness criteria it is possible to ﬁnd a wide but ﬁnite range of non-parametric models
consistent with the data.

A limitation of the Saha & Williams algorithm is that the constraint of a positive-deﬁnite
surface density is weaker than the constraint of a positive 3-d mass density — or, better yet, a
positive-deﬁnite quasi-equilibrium distribution function. In other words, while parametric models
may provide too little freedom, the Saha & Williams non-parametric models provide too much. For
example, their method does not determine whether models are consistent with stellar dynamics,
and many of the models found by the method probably are not consistent. The gravlens software
includes only parametric models that have physical motivations apart from lensing, but it oﬀers
the ability to combine them in arbitrary ways to achieve extensive and reasonable complexity.

– 3 –

3. Selecting Models

Selecting the class of models to use for a particular application is a key part of lens modeling.
When modeling lenses produced by galaxies, a simple and useful place to start is an isothermal
model, i.e. a model with density ρ ∝ r−2 and a ﬂat rotation curve. Spiral galaxy rotation curves
(e.g., Rubin, Ford & Thonnard 1978, 1980), stellar dynamics of elliptical galaxies (e.g., Rix et al.
1997), X-ray halos of elliptical galaxies (Fabbiano 1989), models of some individual lenses (e.g.,
Kochanek 1995; Cohn et al. 2001), joint lensing and dynamical analyses (e.g., Treu & Koopmans
2002; Koopmans & Treu 2003), and lens statistics (e.g., Maoz & Rix 1993; Kochanek 1993, 1996)
are all consistent with roughly isothermal proﬁles.

However, an isolated isothermal ellipsoid rarely yields a good quantitative ﬁt to observed
lenses (e.g., Keeton et al. 1997, 1998; Witt & Mao 1997). In general, adding parameters to the
radial proﬁle of the galaxy fails to produce a good ﬁt, but adding parameters to the angular
structure of the potential dramatically improves the ﬁt (e.g., Keeton & Kochanek 1997; Keeton et
al. 1997). The additional angular structure comes from the tidal perturbations of objects near the
main lens galaxy or along the line of sight. In other words, the fact that few galaxies are truly
isolated means that lens models generically require two independent sources of angular structure:
an ellipsoidal galaxy plus external perturbations. The combination of angular terms can make
it diﬃcult to disentangle the shape of the galaxy and the nature of the external perturbations,
and it is extremely important to understand any degeneracies between the two sources of angular
structure before drawing conclusions from the models (see Keeton et al. 1997).

To move beyond isothermal models and explore other radial proﬁles, softened power law lens
models have traditionally been very popular. However, these models have ﬂat cores, while many
early-type galaxies have cuspy luminosity distributions (e.g., Faber et al. 1997), and dark matter
halos in cosmological simulations have cuspy mass distributions (e.g., Navarro, Frenk & White
1996, 1997; Moore et al. 1998, 1999). The lack of central or “odd” images in most observed galaxy
lenses also limits the extent to which galaxies can have ﬂat cores (e.g., Wallington & Narayan 1993;
Rusin & Ma 2001). Cohn et al. (2001) thus argue that softened power law models are outdated and
should be replaced with families of cuspy lens models. There are three traditional families of models
in which the proﬁle of the cusp is ﬁxed: NFW (Navarro, Frenk & White 1996, 1997), Hernquist
(1990), and de Vaucouleurs (1948) models. In addition, there are three families of models in which
the cusp is taken to be an arbitrary power law: Nuker law models (Lauer et al. 1995; Byun 1996;
Keeton 2003), and two variants of NFW models (see §5 for deﬁnitions; Jing & Suto 2000; Keeton
& Madau 2001; Mu˜noz, Kochanek & Keeton 2001; Wyithe, Turner & Spergel 2001).

An important reason to study families of cuspy lens models is to test the prediction from
the Cold Dark Matter paradigm that halos at a wide range of masses are consistent with a uniﬁed
family of halo models (e.g., Navarro et al. 1996, 1997). The distribution of image separations among
the known lenses rules out the hypothesis that galaxies and clusters can be described by identical
lens models (Keeton 1998; Porciani & Madau 2000; Kochanek & White 2001). However, it may be

– 4 –

possible to resolve the conﬂict with a model in which galaxies and clusters start out with similar
halo proﬁles, but baryonic processes such as cooling modify the inner proﬁles of galaxies. Kochanek
& White (2001) show that such a model can match the observed image separation distribution,
but it remains to be seen whether the model agrees with the lack of central or “odd” images in
observed lenses, or with detailed models of individual lenses.

4. General Equations

For a mass distribution with surface mass density κ(x) = Σ(x)/Σcr in units of the critical
surface density for lensing, the two-dimensional lensing potential is (e.g., Schneider, Ehlers & Falco
1992)

φ(x) =

1
π (cid:90) ln |x − y| κ(y) dy .

(1)

The other lensing properties can be derived from the potential. The deﬂection angle ∇φ determines
the positions of images via the lens equation,

where u is the source position. The magniﬁcation tensor,

u = x − ∇φ(x) ,

µ ≡ (cid:18)

−1

∂u
∂x

(cid:19)

= (cid:34)

1 − φxx −φxy
1 − φyy
−φxy

−1

(cid:35)

,

(2)

(3)

determines the distortions and brightnesses of images. (Subscripts denote partial diﬀerentiation,
φij ≡ ∂2φ/∂xi∂xj.) Many lensing applications involve only the locations and brightnesses of the
images, and thus require only the deﬂection and magniﬁcation components. Applications that
involve the time delays, such as lensing measurements of the Hubble constant H0, also require the
potential.

If the mass distribution has circular symmetry, the deﬂection vector is purely radial and has

an amplitude given by the 1-d integral

φr(r) =

r

2
r (cid:90)

0

u κ(u) du =

1
πΣcr

Mcyl(r)
r

,

(4)

where Mcyl(r) is the mass enclosed by a cylinder of radius r (the projected mass), which is often
easily evaluated. The potential and magniﬁcation components can be obtained by integrating or
diﬀerentiating φr.

More general is the case of elliptical symmetry, in which the surface mass density has the form

κ = κ(ξ), where

ξ2 = x2 + y2/q2 ,

(5)

where q is the projected axis ratio and ξ is an ellipse coordinate. This is the functional form in
a coordinate system with the ellipse centered on the origin and aligned along the x-axis; other

– 5 –

coordinate systems can be reached by suitable translation and rotation. With elliptical symmetry
the lensing properties can be written as a set of 1-d integrals (see Schramm 1990, although I have
changed variables in the integrals),

φ(x, y) =

q
2

I(x, y)

φx(x, y) = q x J0(x, y)
φy(x, y) = q y J1(x, y)
φxx(x, y) = 2 q x2 K0(x, y) + q J0(x, y)
φyy(x, y) = 2 q y2 K2(x, y) + q J1(x, y)
φxy(x, y) = 2 q x y K1(x, y)

where the integrals are

1

1

1

I(x, y) = (cid:90)
0

Jn(x, y) = (cid:90)
0

Kn(x, y) = (cid:90)
0

ξ(u)
u

φr (ξ(u))
[1 − (1 − q2)u]1/2
κ

ξ(u)2
[1 − (1 − q2)u]n+1/2

(cid:0)

(cid:1)

du

u κ(cid:48)

ξ(u)2

(cid:0)
[1 − (1 − q2)u]n+1/2

(cid:1)

du

where

ξ(u)2 = u (cid:18)x2 +

y2
1 − (1 − q2)u

(cid:19)

du

(6)

(7)

(8)

(9)

(10)

(11)

(12)

(13)

(14)

(15)

and κ(cid:48)(ξ2) = dκ(ξ2)/d(ξ2). Note from eq. (12) that the potential can be written as an integral over
the circular deﬂection function φr from eq. (4), but φr must be evaluated at the appropriate ellipse
coordinate ξ(u).

All of the previous expressions assume that the surface density κ is known. Some models
have 3-d density distributions for which the projection integral cannot be evaluated analytically.
In this case even a spherical lens model requires computationally expensive double integrals (the
projection integral followed by the lensing integral). However, the double integrals can be rewritten
as follows so the projection integral is replaced by the enclosed mass M (r), which can often be
computed analytically. (The mass M (r) is the mass in spheres, which is diﬀerent from Mcyl(r)
in eq. 4.) Writing κ(r) as an integral over ρ(r) — the projection integral — and substituting
ρ(r) = M (cid:48)(r)/(4πr2) where M (cid:48)(r) = dM/dr, we ﬁnd
1
2πΣcr
1
2πΣcrr (cid:90)

M (cid:48)(u)
√
u2 − r2
1
1 + y2 (cid:104)M (cid:48)

(cid:17) + M (cid:48)

1 + y−2

κ(r) =

1 + y2

(cid:17)(cid:105) ,

(16)

(17)

(cid:16)r

(cid:16)r

(cid:90)
r

du

dy

(cid:112)

(cid:112)

=

u

∞

1

,

0

where the second line represents a change of variables so the integral has a ﬁnite range, which is
convenient for numerical integration. Combining eqs. (4) and (16) gives the circular deﬂection as
∞

φr =

r

1
πΣcrr (cid:90)

0

du u (cid:90)
u

dv

v

M (cid:48)(v)
√
v2 − u2

,

(18)

– 6 –

r

du u

∞

∂
∂u (cid:90)

u

= − 1

0

πΣcrr (cid:90)
1
πΣcrr

(cid:34)M∞ − (cid:90)
r

=

dv M (cid:48)(v)
√

,

√

v2 − u2
v
v2 − u2
v

(cid:35) ,

∞

dv M (cid:48)(v)

where M∞ is the total mass. Integrating by parts then yields

φr =

=

∞

(cid:90)
r

r
πΣcr
1
πΣcrr (cid:90)

0

dv

v2

1

dy

,

M (v)
√
v2 − r2
1

(1 + y2)3/2 (cid:104)M (cid:16)r

1 + y2

(cid:112)

(cid:17) + y M (cid:16)r

(cid:112)

1 + y−2

(cid:17)(cid:105) ,

(19)

(20)

(21)

(22)

where again the second line represents a change of variables for numerical integration. Eqs. (16)
and (21), or eqs. (17) and (22) for numerical integration, represent the desired 1-d integrals for the
surface density and deﬂection. The magniﬁcation also requires φrr, which could be computed by
diﬀerentiating φr; however, it is easier to compute κ and φr and then use the identity r−1 φr +φrr =
2κ to determine φrr.

5. The Catalog

Table 1 lists a wide variety of popular lens models, all of which are available in the gravlens
software. This section summarizes what is known about the mass distributions and lensing prop-
erties of the models. Analytic results are given where available, which includes all but one of the
circular models and some of the elliptical models. Note that if the potential φ is regular at the
origin, it is normalized to have φ(0) = 0. Lensing is insensitive to an arbitrary constant added to
the potential.

I have tried to include relevant references. If you use results given here, please cite the original
references rather than this catalog. As for the unreferenced results, some of them are new, while
others are derived easily enough that references seem unnecessary. Use your own judgement about
citing such material.

Point mass: This model is inherently circular. A point mass M produces a lensing potential

where the Einstein radius is (in angular units)

φ = R2

E ln r ,

RE = (cid:114)

4GM
c2

Dls
DolDos

.

(23)

(24)

Softened power law potential: Sometimes for analytic simplicity it is convenient to put
the elliptical symmetry in the potential rather than in the density (e.g., Blandford & Kochanek

– 7 –

1987; Witt & Mao 1997, 2000; Evans & Hunter 2002),

φ = b

(cid:0)

s2 + x2 + y2/q2

α/2 − b sα .

(cid:1)

The deﬂection and magniﬁcation are given by simple derivatives.

Softened power law ellipsoid: This model has a projected surface density

κ(ξ) =

1
2

b2−α
(s2 + ξ2)1−α/2

,

(25)

(26)

which represents a ﬂat core with scale radius s, and then a power law decline with exponent α
deﬁned such that the mass grows as Mcyl(r) ∝ rα asymptotically. The core radius can be zero if
α > 0. The model gives a softened isothermal model for α = 1, a modiﬁed Hubble model for α = 0,
and a Plummer model for α = −2 (see Binney & Tremaine 1987). The circular model has lensing
properties

φ =

=

φr =

=

1

α2 b2−α rα
− 1
2α

1

α2 b2−α rα
b2−α
αr (cid:104)(cid:0)
b2
r

ln (cid:18)1 +

; 1 − α
2

; − s2

r2 (cid:21) − 1

α

2F1 (cid:20)− α
, − α
2
2
(cid:104)γE − Ψ (cid:16)

b2−α sα

− α

2 (cid:17)(cid:105)

(α > 0, s = 0)

s2 + r2

r2
s2 (cid:19)

(cid:1)

α/2 − sα

(cid:105)

(α (cid:54)= 0)

(α = 0)

b2−α sα ln (cid:16)

r
s (cid:17)

(27)

(28)

(29)

(30)

In the potential, 2F1[a, b; c; x] is a hypergeometric function, which can be written as or transformed
into a quickly converging series (see Press et al. 1992; Gradshteyn & Ryzhik 1994, §9.1). Also,
γE = 0.577216 . . . is Euler’s constant, and Ψ(x) = d[ln Γ(x)]/dx is the digamma function, or the
logarithmic derivative of the factorial function Γ(x). Analytic solutions for the elliptical model are
possible for α = 0, ±1, and two of these are given below. Barkana (1998) gives a fast numerical
algorithm for general softened power law ellipsoid models.

Isothermal ellipsoid, α = 1: This model describes mass distributions with ﬂat rotation

curves (outside the core). Its lensing properties are:

circular: φ = r φr − b s ln (cid:32)

s +

√

s2 + r2
2s

(cid:33)

φr =

b
r (cid:16)(cid:112)

s2 + r2 − s(cid:17)

elliptical: φ = x φx + y φy − b q s ln
+ b q s ln [(1 + q)s]

(ψ + s)2 + (1 − q2)x2
(cid:2)

1/2
(cid:3)

(31)

(32)

(33)

– 8 –

φx =

φy =

b q
1 − q2

b q
1 − q2

(cid:112)

(cid:112)

tan−1 (cid:34) (cid:112)

1 − q2 x
ψ + s

(cid:35)

tanh−1 (cid:34) (cid:112)

1 − q2 y
ψ + q2s

(cid:35)

(34)

(35)

where ψ2 = q2(s2 + x2) + y2. The elliptical solutions have been given by Kassiola & Kovner
(1993), Kormann, Schneider & Bartelmann (1994), and in the simple form quoted here by Keeton
& Kochanek (1998). In the limit of a singular (s = 0) and spherical (q = 1) model, b is the Einstein
radius of the model and is related to the 1-d velocity dispersion σ by

b = 4π (cid:16)

σ
c (cid:17)

2 Dls
Dos

(36)

(in angular units).

α = −1 ellipsoid: This model corresponds to an unnamed density proﬁle with Σ ∝ r−3

(ρ ∝ r−4) asymptotically. Its lensing properties are:

circular: φ =

φr =

elliptical: φ =

φx =

φy =

b3
sr
b3q
s
b3qx
sψ
b3qy
sψ

s +

ln (cid:32)

b3
s

√

s2 + r2
2s

(cid:33)

(cid:20)1 −

s√

s2 + r2 (cid:21)

ln

(ψ + s)2 + (1 − q2)x2
(cid:2)

ψ + q2s
(ψ + s)2 + (1 − q2)x2
ψ + s
(ψ + s)2 + (1 − q2)x2

1/2 − b3q
s
(cid:3)

ln [(1 + q)s]

(37)

(38)

(39)

(40)

(41)

where ψ2 = q2(s2 + x2) + y2. The elliptical solutions are given by Keeton & Kochanek (1998).

Pseudo-Jaﬀe ellipsoid: A standard Jaﬀe (1983) model has a 3-d density distribution ρ ∝
r−2(r + a)−2 where a is the break radius. For lensing it is useful to modify this model and write
ρ ∝ (r2 + s2)−1(r2 + a2)−1, where a is again the break radius and we have added a core radius
s < a. The projected surface density of the elliptical model has the form

κ(ξ) =

b
2

(cid:34)

(cid:112)

1
s2 + ξ2

−

(cid:112)

1
a2 + ξ2

(cid:35) ,

(42)

which is constant inside s, falls as R−1 between s and a, and falls as R−3 outside a; the total mass
is M = πΣcrqb(a − s). Eq. (42) deﬁnes the pseudo-Jaﬀe ellipsoid. In the limit a → ∞ it reduces
to the isothermal ellipsoid (α = 1). In the limit a → s it reduces to the α = −1 ellipsoid, although
the limit must be taken in ρ rather than in κ (i.e., the limit must be taken before the projection
integral is evaluated). The pseudo-Jaﬀe model is equivalent to a combination of two softened

– 9 –

isothermal ellipsoids, so its lensing properties can be computed with appropriate combinations of
eqs. (33)–(35).

King model: The King model can be approximated as a combination of two softened isother-

mal models (see Young et al. 1980; Barkana et al. 1999),

κ(ξ) =

2.12 b
0.75r2

s + ξ2

−

1.75 b
2.99r2

s + ξ2

.

(cid:112)

(cid:112)

(43)

It has a single scale radius rs. This approximation is convenient because it is written as the
diﬀerence of two softened isothermal ellipsoids, so its lensing properties can be computed with
appropriate combinations of eqs. (33)–(35).

De Vaucouleurs model: This is the prototypical constant mass-to-light ratio lens model (de

Vaucouleurs 1948), with surface mass density

κ(ξ) = κ0 exp (cid:104)

−k(ξ/Re)1/4

(cid:105) ,

(44)

where k = 7.66925001 and Re is the major-axis eﬀective (or half-mass) radius. The circular
deﬂection is (Maoz & Rix 1993)

φr = κ0

40320
k8

R2
e
r

(cid:34)1 −

(45)

e−ζ (cid:18)1 + ζ (cid:18)1 +

ζ
2

(cid:18)1 +

ζ
3

(cid:18)1 +

ζ
4

(cid:18)1 +

ζ
5

(cid:18)1 +

ζ
6

(cid:18)1 +

ζ
7

(cid:19)(cid:19)(cid:19)(cid:19)(cid:19)(cid:19)(cid:19)(cid:35),

where ζ = k (r/Re)1/4. The elliptical model can be computed numerically with eqs. (6)–(11).

Hernquist model: The Hernquist (1990) model is a 3-d density distribution with a projected

distribution that mimics the luminosity distribution of early-type galaxies. It has the form

ρ =

ρs
(r/rs)(1 + r/rs)3 ,

(46)

where rs is a scale length and ρs is a characteristic density. The projected surface mass density has
the form

where x = r/rs, κs = ρs rs/Σcr, and F(x) is the function

κ(r) =

κs
(x2 − 1)2 (cid:2)

−3 + (2 + x2)F(x)
(cid:3)

,

√

F(x) =




1√

x2−1

tan−1
1−x2 tanh−1

1√

x2 − 1
√
1 − x2

(x > 1)

(x < 1)
(x = 1)

1



F (cid:48)(x) =

1 − x2F(x)
x(x2 − 1)

.

A useful technical result is the derivative of this function,

(47)

(48)

(49)

The circular potential and deﬂection are

– 10 –

φ = κs r2

s (cid:20)ln

x2
4

+ 2F(x)(cid:21) ,

φr = 2 κs rs

x[1 − F(x)]
x2 − 1

,

(50)

(51)

where again x = r/rs. The elliptical model κ(ξ) can be computed numerically with eqs. (6)–(11).

NFW model: Cosmological N -body simulations (e.g., Navarro et al. 1996, 1997) suggest that

dark matter halos can be described by a “universal” density proﬁle with the form

ρ =

ρs
(r/rs)(1 + r/rs)2 .

(52)

For the spherical NFW model, the projected surface mass density, potential, and deﬂection are
(Bartelmann 1996; Meneghetti, Bartelmann & Moscardini 2001, 2003; Golse & Kneib 2002)

κ(r) = 2 κs

φ = 2 κs r2

φr = 4 κs rs

,

1 − F(x)
x2 − 1
s (cid:104)ln2 x
ln(x/2) + F(x)
x

2

− arctanh2

,

1 − x2

(cid:105) ,

(cid:112)

(53)

(54)

(55)

where x = r/rs, κs = ρs rs/Σcr, and the function F(x) is the same as in the Hernquist model. A
model with an elliptical density κ(ξ) can be computed numerically with eqs. (6)–(11). Alternatively,
a model with an elliptical potential can be deﬁned by using φ(ξ) and obtaining the deﬂection and
magniﬁcation components by simple diﬀerentiation (see Meneghetti et al. 2001, 2003; Golse &
Kneib 2002).

Cuspy NFW model: Moore et al. (1998, 1999) have suggested that the inner cusp of the
NFW proﬁle is too shallow, so Jing & Suto (2000), Keeton & Madau (2001), and Wyithe et al.
(2001) have studied a generalized NFW-type proﬁle of the form

ρ =

ρs
(r/rs)γ(1 + r/rs)3−γ ,

(56)

so the central cusp has ρ ∝ r−γ. The projected surface density cannot be computed analytically
even for a spherical halo. For a spherical model, eqs. (16) and (21) allow the surface density and
deﬂection to be written as

κ(r) = 2 κs x1−γ (cid:20)(1 + x)γ−3 + (3 − γ) (cid:90)

0

1

dy (y + x)γ−4

(cid:16)1 −

1 − y2

(cid:17)(cid:21) ,

(cid:112)

φr = 4 κs rs x2−γ ×

(57)

(58)

(cid:40)

1
3 − γ 2F1[3 − γ, 3 − γ; 4 − γ; −x] + (cid:90)

0

1

dy (y + x)γ−3 1 −

1 − y2
(cid:112)
y

(cid:41),

– 11 –

where x = r/rs, κs = ρs rs/Σcr, and 2F1 is the hypergeometric function. The integrals can be
evaluated for several values of γ to obtain simple analytic expressions. For γ = 1 the model reduces
to the standard NFW model (see eqs. 53 and 55). Two other analytic cases are as follows:

γ = 2 :

κ(r) = κs (cid:104)

π
x

φr = 4 κs rs (cid:20)

γ = 0 :

κ(r) = κs

− 2F(x)(cid:105) ,
π
2

1
x

+

ln (cid:16)
1 + x2[2 − 3F(x)]
(x2 − 1)2
x
2 (cid:17) +

(cid:20)2 ln (cid:16)

2 κs rs
x

,

φr =

x
2 (cid:17) +

1 − x2
x

F(x)(cid:21) ,

x2 + (2 − 3x2)F(x)
1 − x2

(59)

(60)

(61)

(62)

(cid:21) .

Elliptical NFW-cusp models κ(ξ) can be computed numerically with eqs. (6)–(11). In general

they are slow to compute because they require double integrals.

Cuspy halo models: To obtain a general cuspy model that is more amenable to lensing,

Mu˜noz et al. (2001) introduce a model with a proﬁle of the form

ρ =

ρs
(r/rs)γ[1 + (r/rs)2](n−γ)/2

,

(63)

where again rs is a scale length, and γ and n are the logarithmic slopes at small and large radii,
respectively. This model is a subset of the models whose physical properties were studied by Zhao
(1996). The central cusp must have γ < 3 for the mass to be ﬁnite. For (γ, n) = (1, 4) this is a
pseudo-Hernquist model, for (1, 3) it is a pseudo-NFW model, and for (2, 4) it is a singular pseudo-
1 + (r/rs)2 does not greatly
Jaﬀe model. Compared with eq. (56), replacing (1 + r/rs) with
change the proﬁle shape but does make it possible to solve the spherical model analytically (Mu˜noz
et al. 2001),

(cid:112)

κ(r) = κs B (cid:18)

n − 1
2

n (cid:54)= 3 :

φr = 2

κsrs
x

(cid:40)B (cid:18)

,

(cid:19)

1
2
(cid:0)
n − 3
2

,

(1−n)/2

1 + x2

(cid:1)

2F1 (cid:20)

n − 1
2

,

γ
2

;

n
2

;

1
1 + x2 (cid:21) ,

3 − γ
2

(cid:19)

(64)

(65)

−B (cid:18)

n − 3
2

,

3
2

(cid:19)

(cid:0)

1 + x2

(cid:1)

(3−n)/2

2F1 (cid:20)

n − 3
2

,

γ
2

;

n
2

;

1
1 + x2 (cid:21)(cid:41),

where x = r/rs, κs = ρs rs/Σcr, 2F1 is the hypergeometric function, B(a, b) = Γ(a)Γ(b)/Γ(a + b) is
the Euler beta function. For n = 3 eq. (65) is not valid (it has singularities of the form 0/0), so it
must be replaced with one of the two following expressions:

n = 3 :

φr = 2κsrs(cid:40)

1
x

ln(1 + x2) − G (cid:20)

γ
2

,

γ − 1
2

;

x2
1 + x2 (cid:21)

(66)

– 12 –

−x2−γ(1 + x2)(γ−3)/2 B (cid:18)

γ − 3
2

,

3
2

(cid:19) 2F1 (cid:20)

3
2

,

3 − γ
2

;

= 2

κsrs
x

(cid:26)ln(1 + x2) − G (cid:20)

γ
2

,

3
2

;

1

1 + x2 (cid:21) + Ψ (cid:18)

3
2

(cid:19) − Ψ (cid:18)

;

x2
1 + x2 (cid:21)(cid:41)

5 − γ
2
3 − γ
2

(cid:19)(cid:27)

where Ψ(x) is the digamma function (see eq. 27), and G(b, c; z) is the function

G(b, c; z) ≡ lim
a→0
b
c
∞

=

z +

=

(cid:88)
j=1

2F1[a, b, c; z] − 1
a
b(b + 1)
c(c + 1)

+

z2
2
b(b + 1) · · · (b + j − 1)
c(c + 1) · · · (c + j − 1)

b(b + 1)(b + 2)
c(c + 1)(c + 2)
zj
j

z3
3

+ · · ·

(67)

(68)

(69)

(70)

Note that eqs. (66) and (67) are equivalent and can be transformed into each other using identities
for the hypergeometric function (see Gradshteyn & Rhyzik 1994, §9.131). It is convenient to use
eq. (66) for x < 1 and eq. (67) for x > 1 in order to make the series for the G function converge
rapidly.

Elliptical cusp models κ(ξ) can be computed numerically with eqs. (6)–(11), or with a Fourier

series solution (Chae 2002).

Nuker law: Many early-type galaxies have surface brightness proﬁles that can be modeled as

a Nuker law (e.g., Lauer et al. 1995; Byun et al. 1996),

I(r) = 2(β−γ)/α Ib (cid:18)

−γ

r
rb

(cid:19)

(cid:20)1 + (cid:18)

α

(γ−β)/α

(cid:19)

(cid:21)

,

r
rb

(71)

where γ and β are inner and outer power law indices, respectively, rb is the radius where the break
in the power law occurs, α gives the sharpness of the break, and Ib is the surface brightness at the
break radius. If the luminosity distribution has circular symmetry and the mass-to-light ratio is Υ,
the lensing deﬂection is (Keeton 2003)

φr =

21+(β−γ)/α
2 − γ

κb rb (cid:18)

1−γ

r
rb

(cid:19)

2F1 (cid:20)

2 − γ
α

,

β − γ
α

; 1 +

2 − γ
α

; − (cid:18)

α

(cid:19)

(cid:21) ,

r
rb

(72)

where κb = ΥIb/Σcr is the surface mass density at the break radius in units of the critical density for
lensing, and 2F1 is the hypergeometric function. Ellipsoidal models can be computed numerically
with eqs. (6)–(11).

Exponential disk: The projected surface density is

κ(ξ) = q−1 κ0 exp [−ξ/Rd] ,

(73)

which represents a thin, circular disk with intrinsic central density κ0 and scale length Rd, seen in
projection with axis ratio q = | cos i| where i is the inclination angle (such that i = 0◦ is face-on

– 13 –

and i = 90◦ is edge-on). The circular deﬂection is

φr = 2 κ0

R2
d
r

(cid:20)1 − (cid:18)1 +

r
Rd

(cid:19) e−r/Rd(cid:21) .

(74)

The elliptical model κ(ξ) can be computed with eqs. (6)–(11), or it can be approximated with one
or more Kuzmin disks (see Keeton & Kochanek 1998).

Kuzmin disk: The α = −1 ellipsoid can be re-interpreted as the projection of a thin disk,
in which case it corresponds to a Kuzmin (1956) or Toomre (1962) Model I disk; see Keeton &
Kochanek (1998). Its projected surface density is

κ(ξ) = q−1 κ0 r3

s + ξ2
r2

−3/2

,

s (cid:0)
where κ0 is the intrinsic central surface density of the disk, and q = | cos i| is again the projected
axis ratio of the inclined disk. The only diﬀerence between the Kuzmin disk and the α = −1
ellipsoid is the normalization.

(cid:1)

(75)

Uniform density elliptical slab: Blah. The solution is given by Schramm (1990). Interior

to the disk the lens potential is

φ = κ0

ab
a + b

x2
a

(cid:18)

+

y2
b

(cid:19) ,

(76)

and the deﬂection and magniﬁcation components are given by simple derivatives. Exterior to the
disk the solution is more complicated but still analytic:

φx = 2κ0

φy = 2κ0

φxx = κ0

φyy = κ0

φxy = −κ0

x
a(cid:48) ,
y
b(cid:48) ,

ab
a(cid:48) + b(cid:48)
ab
a(cid:48) + b(cid:48)
2
ab
a(cid:48) (cid:20)
a(cid:48) + b(cid:48)
2
ab
b(cid:48) (cid:20)
a(cid:48) + b(cid:48)
∂λ
∂y

a b x
a(cid:48)3b(cid:48)

− x
a(cid:48)2b(cid:48)
− y
a(cid:48)b(cid:48)2

= −κ0

(cid:21) ,

(cid:21) ,

∂λ
∂x
∂λ
∂y
a b y
a(cid:48)b(cid:48)3

∂λ
∂x

,

b2 + λ, and λ is deﬁned by

x2
a2 + λ

+

y2
b2 + λ

= 1 ,

where a(cid:48) =

√

a2 + λ, b(cid:48) =

√

(77)

(78)

(79)

(80)

(81)

(82)

which leads to

λ =

1
2

(cid:26)x2 + y2 − a2 − b2 + (cid:104)(x2 + y2 − a2 − b2)2 + 4(x2b2 + y2a2 − a2b2)(cid:105)(cid:27) .

(83)

Boxy power law potential: Many elliptical galaxies have isophotes that are not strictly
ellipsoidal but rather are “boxy” or “disky” (REFS). Models with boxy or disky isodensity contours

– 14 –

cannot be treated with the formalism of §4. However, putting the boxiness or diskiness in the
isopotential contours produces a model that can be solved analytically. Zhao & Pronk (2001) and
Cardone et al. (2001) consider such models where, for simplicity, they use a self-similar power law
radial proﬁle. The potential then has the form

φ = b rα [1 − ε cos 2(θ − θε)]αβ,

(84)

where α is the radial power law index, ε is related to the ﬂattening of the isopotential contours,
and θε is the orientation angle of the potential. The parameter β controls the boxiness or diskiness:
β = 1/2 corresponds to strict elliptical symmetry, while β < 1/2 produces a boxy potential and
β > 1/2 produces a disky potential. The deﬂection and magniﬁcation components are given by
simple derivatives of the potential.

External perturbations: Objects near the main lens galaxy or along the line of sight often
perturb the lensing potential.
If the perturbation is weak it may be suﬃcient to expand the
perturbing potential as a Taylor series and keep only a few terms. In a coordinate system centered
on the lens galaxy, the expansion to 3rd order can be written as (see Kochanek 1991; Bernstein &
Fischer 1999)

φ ≈ φ0 + b · x +

r2
2 (cid:104)κ − γ cos 2(θ − θγ)(cid:105) +

r3
3 (cid:104)δ cos(θ − θδ) − ε cos 3(θ − θε)(cid:105) + . . .

(85)

The 0th order term φ0 represents an unobservable zero point of the potential and can be dropped.
The 1st order term b · x represents an unobservable uniform deﬂection and can also be dropped.
The 2nd order term κ represents the convergence from the perturbing mass and is equivalent to a
uniform mass sheet with density Σ/Σcr = κ. The only observable eﬀect of this term is to rescale
the time delay(s) by 1 − κ, which leads to the “mass sheet degeneracy” (e.g., Falco, Gorenstein
& Shapiro 1985); hence this term is often omitted from lens models and introduced a posteriori
using independent mass constraints (see, e.g., Bernstein & Fischer 1999). The 2nd order term γ
represents an external tidal shear with strength γ and direction θγ. The 3rd order term δ arises
from the gradient of the surface density κ(x) of the perturber; it has an amplitude δ = (3/4)|∇κ|
and a direction equal to the direction of ∇κ. The 3rd order term ε arises from the m = 3 multipole
moment of the perturbing mass. The constant coeﬃcients (κ, γ, δ, ε) are all evaluated at the position
of the lens galaxy, and the corresponding direction angles are written here as theory angles measured
counter-clockwise from the x-axis.

General multipole terms: Blah.

φ = − εm
m

rα cos(mθ) ,

(86)

where α is the radial power law index, m is the order of the multipole, and εm speciﬁes the amplitude
of the multipole term.

Acknowledgements: I would like to thank Chris Kochanek, Joanne Cohn, Jose Mu˜noz, David
Rusin, Brian McLeod, and Joseph Leh´ar for many helpful discussions about the models and about

– 15 –

the gravlens software. Support for this work has been provided by ONR-NDSEG grant N00014-
93-I-0774, NSF grant AST-9407122, NASA ATP grant NAG5-4062, Steward Observatory, and
Hubble Fellowship grant HST-HF-01141.01-A from the Space Telescope Science Institute, which is
operated by the Association of Universities for Research in Astronomy, Inc., under NASA contract
NAS5-26555.

Barkana, R. 1998, ApJ, 502, 531

REFERENCES

Barkana, R., Leh´ar, J., Falco, E. E., Grogin, N. A., Keeton, C. R., & Shapiro, I. I. 1999, ApJ, 523,

54

Bartelmann, M., & Steinmetz, M. 1995, A&A, 297, 1

Bartelmann, M. 1996, A&A, 313, 697

Bernstein, G., & Fischer, P. 1999, AJ, 118, 14

Binney, J., & Tremaine, S. 1987, Galactic Dynamics (Princeton: Princeton Univ. Press)

Blandford, R. D., & Kochanek, C. S. 1987, ApJ, 321, 658

Byun, Y.-I., et al. 1996, AJ, 111, 1889

Cardone, V. F., Capozziello, S., Re, V., & Piedipalumbo, E. 2001, A&A, 379, 72

Chae, K.-H. 2002, ApJ, 568, 500

Cohn, J. D., Kochanek, C. S., McLeod, B. A., & Keeton, C. R. 2001, ApJ, 554, 1216

de Vaucouleurs, G. 1948, Ann d’Ap, 11 247

Evans, N. W., & Hunter, C. 2002, ApJ, 575, 68

Fabbiano, G. 1989, ARA&A, 27, 87

Faber, S. M., Tremaine, S., Ajhar, E. A., Byun, Y.-I., Dressler, A., Gebhardt, K., Grillmair, C.,

Kormendy, J., Lauer, T. R., & Richstone, D. 1997, AJ, 114, 1771

Falco, E. E., Gorenstein, M. V., & Shapiro, I. I. 1985, ApJ, 289, L1

Golse, G., & Kneib, J.-P. 2002, A&A, 390, 821

Gould, A., & Loeb, A. 1992, ApJ, 396, 104

Gradshteyn, I. S., & Ryzhik, I. M. 1994, Table of Integrals, Series, and Products, Fifth Edition,

ed. A. Jeﬀrey (San Diego: Academic Press)

Hernquist, L. 1990, ApJ, 356, 359

Jaﬀe, W. 1983, MNRAS, 202, 995

Jing, Y. P., & Suto, Y. 2000, ApJ, 529, L69

Kassiola, A., & Kovner, I. 1993, ApJ, 417, 459

– 16 –

Keeton, C. R., Kochanek, C. S., & Seljak, U. 1997, ApJ, 482, 604

Keeton, C. R., & Kochanek, C. S. 1997, ApJ, 487, 42

Keeton, C. R., & Kochanek, C. S. 1998, ApJ, 495, 157

Keeton, C. R. 1998, PhD. thesis, Harvard University

Keeton, C. R., Kochanek, C. S., & Falco, E. E. 1998, ApJ, 509, 561

Keeton, C. R., Falco, E. E., Impey, C. D., Kochanek, C. S., Leh´ar, J., McLeod, B. A., Rix, H.-W.,

Mu˜noz, J. A., & Peng, C. Y. 2000, ApJ, 542, 74

Keeton, C. R. 2001, ApJ, preprint (astro-ph/0102340)

Keeton, C. R., & Madau, P. 2001, ApJ, 549, L25

Keeton, C. R. 2003, ApJ, 582, 17

Kochanek, C. S. 1991, ApJ, 382, 58

Kochanek, C. S. 1993, ApJ, 419, 12

Kochanek, C. S. 1995, ApJ, 445, 559

Kochanek, C. S. 1996, ApJ, 466, 638

Kochanek, C. S., & White, M. 2001, ApJ, 559, 531

Koopmans, L. V. E., & Fassnacht, C. D. 1999, ApJ, 527, 513

Koopmans, L. V. E., & Treu, T. 2003, ApJ, 583, 606

Kormann, R., Schneider, P., & Bartelmann, M. 1994, A&A, 284, 285

Kuzmin, G. 1956, AZh, 33, 27

Lauer, T. R., et al. 1995, AJ, 110, 2622

Maller, A. H., Flores, R. A., & Primack, J. R. 1997, ApJ, 486, 681

Mao, S., & Paczy´nski, B. 1991, ApJ, 374, L37

Maoz, D., & Rix, H.-W. 1993, ApJ, 416, 425

Meneghetti, M., Bartelmann, M., & Moscardini, L. 2001, preprint (astro-ph/0109250)

Meneghetti, M., Bartelmann, M., & Moscardini, L. 2003, mnras, 340, 105

Moore, B., Governato, F., Quinn, T., Stadel, J., & Lake, G. 1998, ApJ, 499, L5

Moore, B., Quinn, T., Governato, F., Stadel, J., & Lake, G. 1999, MNRAS, 310, 1147

Mu˜noz, J. A., Kochanek, C. S., & Keeton, C. R. 2001, ApJ, 558, 657

Navarro, J. F., Frenk, C. S., & White, S. D. M. 1996, ApJ, 462, 563

Navarro, J. F., Frenk, C. S., & White, S. D. M. 1997, ApJ, 490, 493

Porciani, C., & Madau, P. 2000, ApJ, 532, 679

– 17 –

Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. 1992, Numerical Recipes in

C: The Art of Scientiﬁc Computing, Second Edition (New York: Cambridge Univ. Press)

Rix, H.-W., de Zeeuw, P. T., Carollo, C. M., Cretton, N., & van der Marel, R. P. 1997, ApJ, 488,

702

Rubin, V. C., Ford, W. K., & Thonnard, N. 1978, ApJ, 225, L107

Rubin, V. C., Ford, W. K., & Thonnard, N. 1980, ApJ, 238, 471

Rusin, D., & Ma, C.-P. 2001, ApJ, 549, L33

Saha, P., & Williams, L. L. R. 1997, MNRAS, 292, 148

Schechter, P. L., et al. 1997, ApJ, 475, L85

Schneider, P., & Weiss, A. 1991, A&A, 247, 269

Schneider, P., Ehlers, J., & Falco, E. E. 1992, Gravitational Lenses (New York: Springer)

Schramm, T. 1990, A&A, 231, 19

Toomre, A. 1962, ApJ, 138, 385

Treu, T., & Koopmans, L. V. E. 2002, ApJ, 575, 87

Tyson, J. A., Kochanski, G. P., & Dell’Antonio, I. P. 1998, ApJ, 498, L107

Wallington, S., & Narayan, R. 1993, ApJ, 403, 517

Williams, L. L. R., & Saha, P. 2000, AJ, 119, 439

Witt, H. J., & Mao, S. 1997, MNRAS, 291, 211

Witt, H. J., & Mao, S. 2000, MNRAS, 311, 689

Wyithe, J. S. B., Turner, E. L., & Spergel, D. N. 2001, ApJ, 555, 504

Young, P., Gunn, J. E., Kristian, J., Oke, J. B., & Westphal, J. A. 1980, ApJ, 241, 507

Young, P., Gunn, J. E., Kristian, J., Oke, J. B., & Westphal, J. A. 1981, ApJ, 244, 736

Zhao, H. 1996, MNRAS, 278, 488

Zhao, H., & Pronk, D. 2001, MNRAS, 230, 401

This preprint was prepared with the AAS LATEX macros v5.0.

– 18 –

Table 1. Mass Models for Lensing

Model

Nr

Density ρ(r)

Surface Density κ(r)

Point mass

Power law or α-models

Isothermal (α = 1)
α = −1
Pseudo-Jaﬀe

King (approximate)

de Vaucouleurs
Hernquist
NFW
Cuspy NFW

Cusp
Nuker
Exponential disk

Kuzmin disk

0

2

1

1

2

1

1
1
1
2

3
4
1

1

(cid:0)

δ(x)
s2 + r2
s2 + r2
s2 + r2
s2 + r2
(cid:0)
· · ·

(cid:0)

(cid:0)

(α−3)/2

−1

−2

−1

(cid:1)

(cid:1)

(cid:1)

(cid:1)

a2 + r2

(cid:0)

−1

(cid:1)

−3

· · ·
r−1 (rs + r)
−2
r−1 (rs + r)
r−γ (rs + r)γ−3
r−γ
s + r2
r2
· · ·
· · ·

(cid:1)

(cid:0)

(γ−n)/2

· · ·

(α−2)/2

(cid:1)

(cid:1)

−1/2

−3/2

(cid:0)

(cid:0)

δ(x)
s2 + r2
s2 + r2
s2 + r2
s2 + r2
(cid:0)
2.12

(cid:0)

a2 + r2
−1/2

−1/2

(cid:1)

−1/2

(cid:1)

(cid:1)
0.75r2

−1/2 −
(cid:0)
s + r2
(cid:1)
s + r2
2.99r2
−7.67(r/Re)1/4
(cid:2)

(cid:0)

(cid:0)
−1.75

(cid:1)
(cid:3)

exp
see eq. (47)
see eq. (53)
see eq. (57)

see eq. (64)
see eq. (71)
exp[−r/Rd]
s + r2
r2

−3/2

(cid:0)

(cid:1)

Note. — Density proﬁles for lensing mass models; see §5 for detailed deﬁnitions,
including normalizations. Three-dimensional density proﬁles are not given for the King,
de Vaucouleurs, Nuker, exponential disk, and Kuzmin disk models because these models
are deﬁned by their surface densities. Column 2 indicates the number of parameters
associated with the radial proﬁle alone; each model would also have parameters for the
position and the mass scale, and elliptical models would have parameters for the ellipticity
and orientation. The proﬁles are given for spherical models; elliptical models are deﬁned
by κ(ξ) where ξ is an ellipse coordinate (see eq. 5). add boxypot


