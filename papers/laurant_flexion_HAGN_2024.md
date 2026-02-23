Laurent MAGRI-STELLA

Cosmic Bananas:
Measuring weak gravitational flexion in
Horizon-AGN hydrodynamical
cosmological simulation

MASTER’S THESIS
to achieve the university degree of
Fundamental Physics National Master Diploma

submitted to
Aix-Marseille Université

Supervised by
Raphaël GAVAZZI
Chargé de Recherche, CNRS
Laboratoire d’Astrophysique de Marseille
LAM (UMR7326)

Marseille, March - July 2024

Every model is wrong, but some are useful – George E.P. Box

Contents

1 Introduction

.
.
1.1 Motivation .
1.2 Research Questions and Contributions
.
1.3 Thesis Structure and Organisation .

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

2 Theory and formalism of weak gravitational lensing

2.1 A quick historical review of the first derivations and proofs of gravitational
.
.
.
.
.
2.2 The bases of gravitational lensing .

lensing .

.
.

.
.

.
.

.
.

.
.

.
.

.
.

.
.

.
.

.
.

.
.

.
.

.
.

.
.

.
.

.
.

.
.

.
.

.
.

.
.

.
.

.
.

.

.

.

.

.

.

.

.

.

.

.

.

3 Simulation, data and methods

3.1 The Horizon-AGN simulation .
.
3.2 Data definition and analysis
.
.
3.3 Mass density profiles

.
.
.
.
.
.
The Singular Isothermal Sphere (SIS) model
.
The Navarro-Frenk-White (NFW) model
.
.

3.4 Numerical lensing measurement .

3.3.1
3.3.2

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.

.

.

.

.

.

.

.

.

.

.
.
.
.
.
.

.

.
.

.
.
.

.
.
.

.
.
.

.
.
.

4 Results and interpretation
.
4.1 Tools and discussion .
.
.

4.1.1
4.1.2

.
Julia .
.
.
.
Error estimation .

.
.
.
.
.
.
.
.
.
.
.
4.2 Thin lens planes and lensing maps .
4.3 Numerical correlation computation and results .
.
4.4 The difference between satellite and central galaxies .
.
4.5
.
.

Fitting the SIS and NFW models .
.
Central galaxies .
4.5.1
.
Satellite galaxies
4.5.2

.
.
.
.
.

.
.
.
.
.

.
.
.
.

.
.
.
.

.
.
.
.

.
.
.
.

.
.
.
.

.
.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.

.
.

.
.

.
.

5 Euclid, SourceExtractor++ and Machine Learning
.
.
.

5.1 The Euclid Survey .
SourceExtractor++ .
5.2
.
.
5.3 ONNX .

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.
.
.

.

.

.

.

.
.
.

.
.
.

.
.
.

6 Conclusion and Future Work

Bibliography

iii

.
.
.
.
.
.

.
.
.
.
.
.
.
.
.

.
.
.

.
.
.
.
.
.

.
.
.
.
.
.
.
.
.

.
.
.

.
.
.
.
.
.

.
.
.
.
.
.
.
.
.

.
.
.

.
.
.
.
.
.

.
.
.
.
.
.
.
.
.

.
.
.

.
.
.
.
.
.

.
.
.
.
.
.
.
.
.

.
.
.

.
.
.
.
.
.

.
.
.
.
.
.
.
.
.

.
.
.

.
.
.
.
.
.

.
.
.
.
.
.
.
.
.

.
.
.

.
.
.
.
.
.

.
.
.
.
.
.
.
.
.

.
.
.

.
.
.
.
.
.

.
.
.
.
.
.
.
.
.

.
.
.

.
.
.
.
.
.

.
.
.
.
.
.
.
.
.

.
.
.

.
.
.
.
.
.

.
.
.
.
.
.
.
.
.

.
.
.

.
.
.
.
.
.

.
.
.
.
.
.
.
.
.

.
.
.

.
.
.
.
.
.

.
.
.
.
.
.
.
.
.

.
.
.

.
.
.

.
.

.
.
.
.
.
.

.
.
.
.
.
.
.
.
.

.
.
.

1
1
2
2

4

4
5

8
8
9
13
13
14
16

19
19
19
19
20
22
24
26
26
27

28
28
28
29

30

31

.
.

.
.

.
.

32
32
37

41

42

43

Contents

7 Appendix

A
B

.
Julia code for computing the lensing maps
Julia code for computing the galaxy-galaxy lensing correlation function .

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

List of Figures

List of Tables

List of Listings

iv

CHAPTER 1

Introduction

In our vast Universe, where galaxies twinkle and space-time bends, gravitational lensing emerges
as a compelling artifact predicted by the theory of General Relativity, even if it was also inscribed
in Newtonian theory. It is now regarded as one of the best cosmological probes, providing insight
into the intricate interaction between the nature of matter and gravity on cosmic scales. Among
its various manifestations, weak gravitational flexion stands out as a promising avenue for
exploration, holding the potential to elucidate the nature of dark matter. This Master thesis
stands as a short investigation into weak gravitational flexion, with the goal of refining our
understanding and leveraging it to constrain dark matter models.

While gravitational lensing has a rich historical context, with significant attention traditionally
directed towards shear, the exploration of flexion is an interesting deepening of this domain.
Our research builds upon existing knowledge, seeking not to introduce novel concepts but
rather to delve deeper into established frameworks. By extending our analysis from shear to
flexion, our aim is to extract additional information regarding the distribution of matter in the
universe, particularly concerning its implications for dark matter.

1.1 Motivation

The Euclid satellite, led by the European Space Agency (ESA), aims to map the geometry
of the dark universe by observing an unprecedented amount of galaxies over a third of the
sky. Its primary goal is to investigate dark energy and dark matter, which make up about
95% of the universe’s mass-energy content. Equipped with a visible imaging channel (VIS)
and a near-infrared spectrometer and photometer (NISP), Euclid will measure the shapes and
redshifts of distant galaxies with a precision never achieved before. This will provide crucial
data for understanding the universe’s expansion history and the growth of cosmic structures.
Gravitational lensing, including shear, convergence, and especially flexion which is the main

focus of this work, is essential for Euclid’s mission.

1

Chapter 1: Introduction

• Shear measures the distortion of galaxy shapes due to gravitational fields, providing

insights into the distribution of dark matter.

• Convergence refers to the magnification of light from distant objects, helping to map

mass distribution, including dark matter halos.

• Flexion captures higher-order distortions, offering more detailed information about the

gravitational field and substructures within dark matter halos.

By analyzing these lensing effects, Euclid will test theories of gravity on cosmic scales and
enhance our understanding of the universe’s fundamental properties.

Our exploration is guided by the overarching objective of constraining dark matter models,
using insights obtained from weak gravitational flexion to refine our comprehension of this
matter constituent. By observing this phenomenon, in simulated data as well as images from
space satellites, we hope to gain insights into the distribution and properties of dark matter,
thereby advancing our pursuit to unravel its mysteries.

To achieve these objectives, we adopt an interdisciplinary approach, selecting and using
tools and techniques from theoretical physics, computational modeling, and observational
cosmology. Through the integration of numerical simulations, statistical analyses, and ob-
servational data, we aim to build an efficient and complete study of the data provided.

1.2 Research Questions and Contributions

The primary questions we aim to understand in this work are the following :

• What insights can weak gravitational flexion offer into the distribution and properties

of dark matter on cosmic scales?

• What empirical results can be obtained to validate the compatibility of current models

with observational reality?

• To what extent can we expect the incorporation of flexion to enhance our results,

precision and interpretations from what we could do before using only shear?

1.3 Thesis Structure and Organisation

We will start by going through the theoretical grounds needed to work through gravitational
lensing and how does one express the different terms in the Jacobian to highlight what we
need to measure, and how we will do so. We will present the Horizon-AGN hydrodynamical
cosmological simulation, which cosmology it uses, as well as the numerical methods it uses to
generate a volume of Universe. We will describe the way we used Horizon-AGN to generate
some data to work on, explaining the nature of this data in the mean time. We will then go
on to present the methodology of our analysis, the tools and methods used and finally present
and interpret the results. We will also include a short section on how the internship work will

2

Chapter 1: Introduction

end, as it carries on for a bit after the report deadline, explaining how we intend to adapt this
whole analysis to Euclid by using Machine Learning tools to modify the well known and used
software SourceExtractor to include the flexion components to its parameter fitting list.

3

CHAPTER 2

Theory and formalism of weak
gravitational lensing

2.1 A quick historical review of the first derivations and

proofs of gravitational lensing

Quite notably for gravitational lensing, theory came first and then the observational proofs
followed, which is unusual for astronomical research, where observation is usually the first
step to understand a physical phenomenon.

Gravitational lensing, funnily enough, can be treated under a Newtonian point of view,
considering the movement of a light corpuscle (which is a theory supported by Newton
himself) along a hyperbolic trajectory of a test particle, travelling at the speed of light under
the action of Newtonian gravity. Note that under those assumptions, we are not talking
about photons, as those are only defined in relativistic field theories, but about light corpuscles
moving at a speed v = c. Considering a point mass of mass M responsible for the force felt by
a mass m test particle,

⃗F = −

GMm
r2

ˆr

one can write the velocity boost felt by the corpuscle as

∆⃗v = −

(cid:90) t

0

dτ

GM
r(τ)2 ˆr(τ)

(2.1)

(2.2)

which, as in thoroughly developed in [1], once worked out, considering a constant speed
of light on the parallel direction of travel and that the deflection of light is small, gives the
following result for the angle of deflection of the light corpuscle,

ΨN =

2GM
c2r

which is half the result predicted by Einstein’s theory of General Relativity

ΨGR =

4GM
c2r

(2.3)

(2.4)

4

Chapter 2: Theory and formalism of weak gravitational lensing

The Newtonian result was supposed by Newton in 1704 and predicted by Henry Cavendish and
Johann Georg von Soldner between 1784 and 1804, while the GR prediction of the value was
computed by Einstein in 1911, while completing his theory of relativity. The first observation
of the bending of the path of light was made by Eddington and Dyson in 1919 during a solar
eclipse, which allowed for stars near the Sun to be seen. Those observations showed that the
positions of those stars was slightly out of place, proving that the light was being lensed by
the Sun.

From this point, gravitational lensing was a bit forgotten and thought of as a "fun" artifact
from general relativity, studied and theorized by many physicists as Zwicky or Klimov. Then
came the 80s, when astronomers understood that using CCD imaging, they could measure
data from millions of stars and observe so-called microlensing events.

2.2 The bases of gravitational lensing

Whenever a physicist thinks of lenses, they might get reminded of geometrical optics classes,
as such, gravitational lensing also has its own lens equation.

Figure 2.1: Definition of lensing angles and distances. Figure by Michael SACHS

Using 2.1, we can write the following equation, which is just a difference of angles, up to a

factor of distances ratio, called the lens equation,

⃗β = ⃗θ −

Dds
Ds

ˆα( ⃗Ddθ )

(2.5)

we can use (2.5) to write the Jacobian of the coordinate transformation between image and

5

Chapter 2: Theory and formalism of weak gravitational lensing

source plans, or respectively, lensed and unlensed plans,

∂ βi
∂ θ j

= δi j −

∂ αi
∂ θ j

= Ai j

(2.6)

we can then decompose this matrix as a diagonal term and a trace-free term, which will contain
important quantities for the rest of the work, respectively the scalar convergence κ and the
2-spinor shear γ

A =

(cid:18)1 − κ
0

(cid:19)

0
1 − κ

− γ

(cid:18)cos(2φ )

sin(2φ )

(cid:19)

sin(2φ ) − cos(2φ )

(2.7)

where φ is the angle between ⃗α and the x-axis (???) The convergence is responsible for
increasing the size of the lensed object, the image, and conserving its surface brightness,
therefore magnifying the object and increasing its appearing luminosity. The shear on the
other hand causes more geometrical consequences on the image, it stretches it tangentially,
meaning it seems to align on the gradient of the lensing potential.

Before diving further into how we accurately measure those quantities, let us define two

new lensing fields that will be the core of the work of this report.

As you might have seen, eq.(2.7) is only an approximation to the first order of the matrix
Ai j defined in eq.(2.6). But since convergence and shear are not constant within a galaxy
image, we need to go to the next order to have a more precise description of the coordinate
transformation, i.e. to the second order in deflection angle derivatives, or equivalently, third
order in lensing potential derivatives,

Ai j =

∂ βi
∂ θ j

= δi j −

∂ αi
∂ θ j

−

∂ 2αi
∂ θ j∂ θk

dθk = δi j −

∂ 2ψi
∂ θ j∂ θk

dθk −

∂ 3ψi
∂ θ j∂ θk∂ θl

dθkdθl

(2.8)

As in [2], we can use a complex representation of our fields and a complex gradient operator,

∂ = ∂1 + i∂2

(2.9)

to express conveniently all of our lensing quantities, this also allows us to generate and thus
reexpress our fields in a new manner. Thus, the displacement field α becomes

α = ∂ ψ = α1 + iα2

(2.10)

applying this operator (which increases the field’s spin by 1) and its conjugate ∂ ∗, (which
lowers the spin by 1) on the lensing potential successively generates new expressions for κ
and γ and the new lensing quantities we will be interested in. Thus, the spin-0 convergence
field becomes

κ =

1
2

∂ ∗∂ ψ =

∂ ∗α

1
2

and the spin-2 shear field,

γ =

1
2

∂ ∂ ψ =

1
2

∂ α = γ1 + iγ2

6

(2.11)

(2.12)

Chapter 2: Theory and formalism of weak gravitational lensing

We can now express in a simple way the third derivatives of the lensing potential with this
formalism, giving us new fields F and G

F = ∂ κ =

G = ∂ γ =

∂ ∂ ∗∂ ψ = |F |eiφ

1
2
1
∂ ∂ ∂ ψ = |G |e3iφ
2

(2.13)

(2.14)

using that ∂ κ = ∂ ∗γ because ∂ commutes with its conjugate, eq.(??) and eq.(2.9), we can
express in an exploitable way the two new lensing fields F and G , called respectively 1-flexion
or first flexion, and 3-flexion or second flexion.

F = ∂1γ1 + ∂2γ2 + i(∂1γ2 − ∂2γ1) = F1 + iF2
G = ∂1γ1 − ∂2γ2 + i(∂1γ2 + ∂2γ1) = G1 + iG2

(2.15)
(2.16)

as a result, we are left with a total of 7 independant lensing fields, all expressed from various
derivative combinations of the deflection angle α, namely κ, γ1, γ2, F1, F2, G1 and G2. Those
4 new quantities called flexion terms account for the deformation and "arc-ed" aspect that can
happen when observing gravitational lenses. Studying the behaviour of those fields should
allow us to better understand the link between weak and strong gravitational lensing.

Figure 2.2: Effects of the different lensing fields on a Gaussian galaxy of radius 1 arcsec. 10%
convergence/shear and 0.28 arcsec−1 flexion (which is a very large value for this
quantity, chosen only for visualizing) are applied. Figure extracted from [2]

We will focus our analysis on those 7 quantities, but one could obviously go higher in orders
of derivatives, the next order would be composed of 3 fields of spin 4,2 and 0, and the next of
3 fields of spin 5,3 and 1 and so on.

7

CHAPTER 3

Simulation, data and methods

3.1 The Horizon-AGN simulation

The Horizon-AGN simulation is a hydrodynamical simulation of a volume of the universe.
This simulation, which offers many possibilities to study cosmological and astrophysical
phenomena, presented, for instance, in [3], uses the ramses code from [4], a N-body method
designed to study structure formation across the Universe, it uses AMR (adaptive mesh
refinement) to simulate, in [3] for example, up to 10243 dark matter particles in a square box
of side 100 h−1 Mpc. The mass resolution is thus 8 · 107 h−1M⊙. The chosen cosmology is a
standard flat Λ CDM, and its parameters are listed in the following table.

H0
70.4

Ωb
0.045

Ωdm
0.227

ΩΛ
0.728

σ8
0.81

ns
0.967

Table 3.1: List of parameters used for the flat ΛCDM cosmology defined in the Horizon-AGN
simulation and their values, all parameters are unitless, except for the Hubble
parameter H0, expressed in km s−1Mpc−1

In this simulation, the gas evolution is computed with a second-order unsplit Godunov
scheme for the Euler equation. The interpolated variables are reconstructed using the Harten-
Lax-van Leer-Contact Riemann solver (presented in [5]). The gas is heated via a uniform UV
background that originated from the reionization epoch at z = 10, while its cooling is allowed
via H and He collisions. Schmidt law is used to model star formation, and black hole formation
occurs when the gas density reaches 0.1 H.cm−3.

8

Chapter 3: Simulation, data and methods

Figure 3.1: Two dimensional graph of the past lightcone generated with Horizon-AGN. Each
blue cell is a replicate of the H-AGN simulation box. The orange vertical line marks
z = 1 and the tiling repeats up to z = 8. The black lines represents the aperture angle
of the cone which goes from 2.25 arcseconds to 1.00 arcsecond after z = 1. Figure
extracted from [6]

A raytracing work has been done by C. Gouin et al. in [6], based on similar previous works in
the field such as [7]. The idea is that once while the volume is being simulated, we simulate the
past lightcone of the observer, represented in figure 3.1. Which means we need to trace the path
of light rays that are being deflected by the forefront mass along their journey, from the source
to the observer. As explained in C. Gouin’s article, this can be done in two different ways, one
of them being to integrate the gravitational potential on a thin lens plane, meaning that we
dissect our line of sight into different thin lens planes and then integrate the gravitational
potential included in each of them.

α(θ ) =

(cid:90)

2
c2

∇⊥φ (θ , x)dx

3.2 Data definition and analysis

(3.1)

To work on measuring convergence, shear and flexion in the Horizon-AGN simulation, we
need to properly define what we want to look out and the procedure in order to do so. We are
interested in the galaxy-galaxy lensing statistic. It is a type of weak (although it sometimes
falls in the strong regime) gravitational lensing in which a forefront massive galaxy causes

9

Chapter 3: Simulation, data and methods

distortions on the shapes of background galaxies/stars. This type of lensing can be used for
various interests, for instance :

• Measuring mass density profiles :

Galaxy-galaxy lensing can provide information on the shape of mass profiles. Since
lensing doesn’t respond to the nature of the mass causing it (it is caused by gravity only),
we can probe a very interesting amount of mass environments, like baryons, which
proportion tends to dominate at small scales, at the heart of the galaxies, and dark matter,
which accumulates in the outer parts of the halos.

• Time evolution:

The mass is not constant throughout the life of a galaxy, thus, focusing the study on
different lens populations for different redshifts, we can observe the time evolution of
the mass properties of halos.

• Mass-to-light relation:

Not only can we restrict distances and masses, but we can have a look into the average
luminosity emitted by a stack of galaxies in specific filters. We recall that since lensing
is unaffected by the nature of the matter causing it, the measure would be a virial mass
to light ratio. Using then other probes to trace the ratio of baryonic to dark matter will
give precious information about the distribution of matter in galactic halos.

We are interested the most in the first bullet point: measuring mass density profiles. To do
so, we have at our disposition a selection of the thin lens planes we presented earlier, for each
plane we have 2 x 36 000 x 36 000 pixels, corresponding to the two different components of
the deflection field α = α1 + iα2, over a square grid of 36 000 x 36 000 pixels. We will use these
maps to compute the convergence, shear, and flexion fields. Recalling eqs. (2.11), (2.12) and
(2.16), we express each component in the following table, as derivatives of the deflection field,
we detail for instance also how we obtain the expression for κ.

κ =

=

=

1
2
1
2
1
2

∂ ∗α

(∂1 − i∂2)(α1 + iα2)

[(∂1α1 + ∂2α2) + i(∂1α2 − ∂2α1)]

(3.2)

But since κ is defined to be a scalar, the rotational term or equivalently the complex part of
κ needs to be zero, we obtain a relation whose validity will be tested and discussed, as it was
in [6],

∂1α2 = ∂2α1

(3.3)

Computing in a similar manner, we get all of the following expressions that we need for builing
our lensing quantities maps

10

Chapter 3: Simulation, data and methods

κ =

γ1 =

γ2 =

F1 =

F2 =

G1 =

G2 =

1
2
1
2
1
2
1
2
1
2
1
2
1
2

(∂1α1 + ∂2α2)

(∂1α1 − ∂2α2)

(∂1α2 + ∂2α1) = ∂1α2 = ∂2α1

(∂1∂1α1 − ∂1∂2α2 + ∂2∂2α1 + ∂2∂1α2)

(∂1∂1α2 + ∂1∂2α1 − ∂2∂1α1 + ∂2∂2α2)

(∂1∂1α1 − ∂1∂2α2 − ∂2∂2α1 − ∂2∂1α2)

(∂1∂1α2 + ∂1∂2α1 + ∂2∂1α1 − ∂2∂2α2)

(3.4)

Figure 3.2: Expressions as derivatives of the two components of the deflection fields of the

seven lensing quantities defined earlier and used

To compute the numerical derivatives, we will use a simple second order finite differences
scheme. Now that our quantities are thoroughly defined, we need to correlate the maps (which
vary in redshift, as they are different thin lens planes) to the lensed position (equivalently
position on the image plan) and appearance of the galaxies extracted from a catalog built from
the simulated lightcone defined earlier. To do so, we need to define what we call the tangential
component of all the lensing fields. Talking about tangential convergence has no sense, as it
is a scalar, so its projection on the tangential direction is simply itself. But for all of the other
lensing component, we need to take the real part of each lensing component λi , and we define

λi+ = ℜ(λi · eniπφ ), ∀(λi, ni) ∈ [(γ, 2), (F , 1), (G , 3)]

(3.5)

Proceeding to such a projection, we get the following expressions for each component we need
to measure. It will also be interesting to study the cross lensing components, as opposed to
the tangential ones, which are defined as their imaginary part.

11

Chapter 3: Simulation, data and methods

κ+ = κ
γ+ = −γ1 cos(2φ ) − γ2 sin(2φ ) = −γ1(cos2(φ ) − sin2(φ )) − γ2 · 2 sin(φ ) cos(φ )
F+ = −F1 cos(φ ) − F2 sin(φ )
G+ = −G1(4 cos3(φ ) − 3 cos(φ )) − G2(4 sin3(φ ) − 3 sin(φ ))

(3.6)

Figure 3.3: Tangential projection of the lensing quantities, generated by gravitational lensing,

used to effectively compute and measure the lensing quantities.

Those tangential expressions were obtained by using eq.(3.5), and then reexpressing the
components as functions of φ rather than nφ . The following cross lensing quantities are not
generated by gravitational lensing, so measuring them should give no signal; doing so is a
great way of checking for systematics.

γ× = γ1 sin(2φ ) − γ2 cos(2φ ) = γ12 sin(φ ) cos(φ ) − γ2(cos2(φ ) − sin2(φ ))
F× = F1 sin(φ ) − F2 cos(φ )
G× = −G2

(cid:0)4 cos3(φ ) − 3 cos(φ )(cid:1) + G1

(cid:0)−4 sin3(φ ) + 3 sin(φ )(cid:1)

(3.7)

Figure 3.4: Cross projection of the lensing quantities, which should yield no signal and thus

can be used as a robust way of checking for systematical errors.

Eq. (3.6) can be re-expressed in a final way, using the x,y coordinates on the sky plane,
considered flat because we know that the lightcone is opened with a 2.25° angle until redshift
1, and 1° afterwards, since the angles are small relative to the distances, we can easily compute
the angular diameter distance. Samely with eq. (3.7).

γ+ = −γ1

x2 − y2

r2 − γ2

2xy
r2

(3.8)

F+ =

−F1x − F2y
r
(cid:18)
4

x3
r3 − 3

G+ = −G1

(cid:19)

x
r

− G2

(cid:19)

(cid:18)
4

y3
r3 − 3

y
r

Figure 3.5: Tangential projection of the lensing quantities, expressed in flat sky plane coordi-

nate.

12

Chapter 3: Simulation, data and methods

x2 − y2
r2

γ× = γ1

F× =

2xy
r2 − γ2
F1y − F2x
r
(cid:18)
4

G× = −G2

x3
r3 − 3

(cid:19)

x
r

+ G1

(cid:19)

(cid:18)
4

y3
r3 − 3

y
r

(3.9)

Figure 3.6: Cross projection of the lensing quantities, expressed in flat sky plane coordinates.

All of those formulas allow us to measure correctly the flexion, shear and convergence in

our simulated lightcone data from the Horizon-AGN simulation.

3.3 Mass density profiles

We now need to define the various mass density profiles we will use to fit our gravitational
lensing measurements and thus derive from their original expressions, the formulas for
convergence κ, shear γ and the two flexion fields F and G .

3.3.1 The Singular Isothermal Sphere (SIS) model
This model is by far one of the simplest ones we can study, as it only implies that the mass
density profile scales as the inverse of the square of the distance from the lens’ center.

ρ ∝

1
r2

(3.10)

This type of profile, presented for instance in [8] is obtained by considering constant velocity
dispersion for the dark matter particles all around the matter halo. Defining the projected
surface mass density of the SIS,

Σ(ξ ) =

σ 2
v
2Gξ

(3.11)

where ξ is the distance from the center of the lens and σv the velocity dispersion aforemen-
tioned. Since the convergence is defined as κ = Σ/Σc , with Σc the critical density,

Σc =

c2
4πG

Ds
DlDls

(3.12)

13

Chapter 3: Simulation, data and methods

with Ds, Dl and Dls respectively the angular diameter distance from observer to source,
observer to lens, lens to source, we obtain finally for κ,

κ(θ ) =

θE
2θ

(3.13)

θE is called the Einstein angle and θ is the angular distance to the center of the lens, defined
respectively as,

θ =

ξ
Dl

θE = 4π

Dls
Ds

σ 2
v
c2

As derived in [9], the shear caused by such a mass profile is,

γ(θ ) = −

θE
2θ

e2iφ

(3.14)

(3.15)

where we recall that as is eqs. (3.5), (3.6) and (3.7) φ is the angle around the lens. We can
finally express the flexion components,

F = −

G =

θE
2θ 2 eiφ
3θE
2θ 2 e3iφ

(3.16)

where we recover properly the spin properties we talked about when introducing the two
different flexion fields. We see that fitting using the SIS model, we have only one free
parameter, θE, for the four different lensing fields, we will see the implications and differences
with the next model in section 4.

3.3.2 The Navarro-Frenk-White (NFW) model
This model quickly raises complexity in its expression, and especially when expressing the
lensing components. This model was empirically found and studied by the three scientists
that gave the model its name, Julio Navarro, Carlos Frenk and Simon White in [10]. The mass
density is given by the following formula,

ρ(x, z) =

∆cρc(z)
x(1 + x)2

(3.17)

where x = r/rs is the distance to the center, in units of a characteristic radius rs, one of the free
parameters of the model, ρc(z) is the critical density at redshift z and ∆c is a unit-less scaling

14

Chapter 3: Simulation, data and methods

density which is a function of another free parameter, c the concentration of the halo,

∆c =

200
3

c3
ln(1 + c) − c
1+c

(3.18)

for this profile, the convergence can be written as a complicated function of a new variable
y = ξ /rs

κ(y) = 2κs

f (y)
y2 − 1

(3.19)

where κs is related to the other constants as κs = ρc(z)∆crs/Σc and the function f (y) is given
by

f (y) =






1 − 2√

arctanh

1 − 2√

arctan

1−y2

y2−1

(cid:113) 1−y
1+y
(cid:113) y−1
y+1

y < 1

y > 1.

The shear is given by the following formula,

γ(y) = κsg(y)

where,

g(y) =






√
√

(1−y)/(1+y)
8 arctanh
y2
1−y2
3 + 4 ln( 1
10
2 )
√
(y−1)/(1+y)
8 arctan
√
y2−1

y2

+ 4

y2 ln (cid:0) y

2

(cid:1) − 2

(y2−1) + 4 arctanh

√

(1−y)/(1+y)

(y2−1)(1−y2)1/2

+ 4

y2 ln (cid:0) y

2

(cid:1) − 2

(y2−1) + 4 arctan

(y2−1)3/2

√

(y−1)/(1+y)

and the flexion fields are given by

F (y, φ ) = −

2Fs
(y2 − 1)2 [2y f (y) − h(y)]eiφ

G (y, φ ) = 2F∫

(cid:20) 8
y3 ln
where the functions h(y) and j(y) are given by

(3/y)(1 − 2y2) + j(y)
(y2 − 1)2

y
2

+

(cid:21)

e3iφ

h(y) =





(cid:112)

2y
1 − y2
2y
(cid:112)y2 − 1

arctanh

arctan

(cid:114) 1 − y
1 + y
(cid:114)y − 1
y + 1

−

1
y

,

−

1
y

,

y < 1

y > 1.

15

y < 1

y = 1.

y > 1

(3.20)

(3.21)

(3.22)

(3.23)

(3.24)

Chapter 3: Simulation, data and methods

j(y) =






(cid:16) 8
y3 − 20
(cid:16) 8
y3 − 20

(cid:17)
y + 15y
(cid:17)
y + 15y

2√

1−y2

2√

y2−1

arctanh

arctan

(cid:113) 1−y
1+y
(cid:113) y−1
y+1

if y < 1

if y > 1

(3.25)

and Fs = κsDl/rs It is worth to note that the variables y and φ are decoupled, and since we
are only interested for now in the radial mass profile, we can simply integrate over a full turn
around the lens to "free" the formulas of their angular components.

3.4 Numerical lensing measurement

Let us explain how we measure lensing in practice. Given a catalog of the galaxies present in
the Horizon-AGN lightcone, we have at our disposition a list of parameters for each galaxy.

Description
Main virial mass of the halo
Right ascension of the galaxy image
Total mass of the subhalo
Quality check flag
Declination
Right ascension
Total mass of the main halo
Main halo ID

Column Name
MVIRH_MAIN
RA_IMG
MTOTH_SUB
CHECK
DEC
RA
MTOTH_MAIN
IDH_MAIN
ID_PHOTOMETRY Photometry ID
ES2
E1
GAMMA2
ID
MASS_TRUE
ID_HzAGNFULL
E2
MAGNIF
ES1
z_true
KAPPA
MVIRH_SUB
GAMMA1
DEC_IMG
IDH_SUB

Ellipticity component s2
Ellipticity component e1
Shear component gamma2
Unique galaxy ID
True mass of the galaxy
ID in the HzAGN full simulation
Ellipticity component e2
Magnification due to lensing
Ellipticity component s1
True redshift
Convergence
Virial mass of the subhalo
Shear component gamma1
Declination of the galaxy image
Subhalo ID

Table 3.2: Available information in the header of the galaxy catalog within the Horizon-AGN

simulated light cone.

16

Chapter 3: Simulation, data and methods

Since we are interested in gravitational lensing and the influence of redshift and mass on it,

we will mainly be interested in the following columns :

• RA_IMG, the lensed right ascension coordinate of the galaxy

• DEC_IMG, the lensed declination coordinate of the galaxy

• z_true, the galaxy’s redshift in the simulated lightcone

• MTOTH_MAIN, the total mass of the main dark matter halo in which the galaxy is

contained

• MTOTH_SUB, the total mass of the possible dark matter sub-halo that the galaxy

represents

• IDH_SUB, the unique ID of the aforementioned possible sub-halo

• IDH_MAIN. the unique ID of the main halo containing the galaxy

When we say "possible dark matter sub-halo", we mean that we can divide the galaxies
into two categories, the main and satellite galaxies. The first category represents the galaxies
with the largest halos, that tend to "host" other substructures, so called sub-halos, or satellite
galaxies. To do such a categorization, we compare if for a given galaxy, its IDH_SUB is equal
to its IDH_MAIN, if its equal, then the galaxy is its own main halo, it is a main or host galaxy.
If not, then it is comprised within a larger dark matter halo, it is a satellite galaxy.

We also need to sort the galaxies by redshift. We are given a handful of thin lens planes
which are to be interpreted as an integral of all the forefront mass between the observer and
the corresponding redshift of the thin lens plane in question. A known result in lensing physics
is that the usual deflectors are located for a given source at around half the physical distance,
as specified in [11]. This means that for a given plane, we compute the comoving distance,
half it and then compute the redshift corresponding to this comoving distance, taking a slice
of arbitrary width around this redshift (±0.1 in our case) allows us to minimize noise from
galaxies that wouldn’t interfere with the given lens system, as they would be too far or too
close from the observer, while keeping a significant and satisfying amount of signal to work
with.

Once those two selections are done, we sort by mass, i.e. taking a specific range of galaxy
masses, be it the 500 most massive, or 500 least massive, to see the influence of mass on the
lensing signal. Now, for a given plane, our galaxies are correctly sorted and selected, we can
go on and measure the convergence, shear and flexion signals. To do so, we iterate over the
selected galaxies, and draw concentric circles around each galaxy, as to divide its surrounding
space in radii bins.

17

Chapter 3: Simulation, data and methods

Figure 3.7: The process of measuring the galaxy-galaxy lensing correlation statistic. Figure by

Masamune Oguri

Doing so, we can sort the pixels of the convergence/shear/flexion maps, computed with
Eq. (3.4) into the corresponding bins, and add successively the values of the tangential
components, defined in Eq. (3.8), then dividing each components’ bins by the numbers of
pixels in it. This will allow us to finally have a measurement of gravitational lensing caused
by galaxies on other galaxies by applying this process on all the lensing fields.

18

CHAPTER 4

Results and interpretation

This section will focus on a detailed analysis of measurements resulting from the framework
described earlier. We aim to present the data collected, followed by an interpretation of what
the results signify in the context of the study. Each subsection will address specific checkpoints
of the internship and situate it in regards to the set objective. The purpose is to provide a
comprehensive understanding of the outcomes, bridging the gap between raw data and its
implications for the rest of the work.

4.1 Tools and discussion

4.1.1 Julia
We will use the language Julia for the purpose of our work. Choosing between Julia, Python, C,
or C++ (among others) involves considering factors such as performance, ease of use, available
libraries, the ability to interface between them, and specific use cases. Julia offers high perfor-
mance close to C, thanks to JIT compilation, making it ideal for numerical computation and
scientific computing. It’s user-friendly like Python, with a syntax conducive to mathematical
operations, but its ecosystem, though growing, isn’t as extensive as Python’s.

Python is most likely the best in readability and simplicity, supported by a vast ecosystem
specialized in data analysis, machine learning, and web development. It is versatile, with the
ability to call C/C++ for performance tasks, but it generally offers slower execution speed.

C and C++ are known for their high performance and control over system resources, which
makes them suitable for system-level programming and applications that require maximum
efficiency. However, they have a steeper learning curve and more complex syntax compared
to Julia and Python, which can slow down development and hinder its readability.

For all these reasons, Julia is a strong candidate for projects that demand high computational

performance and mathematical modeling like ours.

4.1.2 Error estimation
To estimate the error bars in the measurements of convergence, shear and flexion, we will use
a Jackknife resampling. The decision between using the jackknife rather than other statistical

19

Chapter 4: Results and interpretation

methods to estimate error bars is based on various factors related to the data and the objectives
of the analysis. The Jackknife method is known for its computational simplicity and speed,
which makes it particularly suitable for large datasets because of its methodical approach of
excluding one (or more, in the case of block Jackknife) observation at a time. This method not
only facilitates bias correction, but is also a simple tool for variance estimation, ensuring our
ability to reproduce our results.

However, one might want to use the bootstrap, which is a pretty flexible tool. It allows
for intensive replacement resampling, creating numerous samples that effectively simulate
drawing new samples from the population. This method is particularly beneficial for complex
estimators, where it can provide more accurate error estimates by capturing the estimator’s
distribution more comprehensively. Despite being computationally intensive, bootstrap meth-
ods are relatively straightforward to implement and can be adapted to various estimators and
statistical measures.

Since our estimators are simple in nature and our samples are large but somewhat reason-
able when it comes to computing power needed to treat them, we will settle on the Jackknife
resampling.

4.2 Thin lens planes and lensing maps

As mentioned earlier, we have at our disposition a selection of thin lens planes, each containing
information such as their redshift, and two 36 000 pixels x 36 000 pixels maps, representing
the two independent deflection fields α1 and α2 , recalling Eq. (2.10). By Eq. (3.4), and using
a second order finite difference differentiation scheme, we can easily compute each of the
lensing maps κ, γ1, γ2, F1, F2, G1 and G2. We recall that since the angles are small the plane is
considered flat, and thus the derivatives are to be taken cartesian and one can understand ∂1,
the derivative with respect to the first coordinate, as ∂ /∂ x the derivative with respect to the x
coordinate, and samely with ∂2 as ∂ /∂ y.

20

Chapter 4: Results and interpretation

Figure 4.1: Convergence map for a redshift zplane ≈ 1.016, computed from a combination of

derivatives of the two deflection field maps.

(a) γ1 shear map, with its visible spin 2 prop-

erty

(b) γ2 shear map, one can notice that it is
seemingly the same as γ1 but rotated 90°

Figure 4.2: The two different shear fields γ1 and γ2

21

Chapter 4: Results and interpretation

(a) F1 flexion map, with its visible spin 1

property

(b) F2 flexion map, one can notice that it is
seemingly the same as F1 but rotated 90°

(c) G1 flexion map, and its noticeable spin 3

property

(d) G2 flexion map, also seems to be like G2 but

rotated 90°

Figure 4.3: The four different flexion fields F1, F2 and G1, G2. Computed by the code in

Appendix A.

These different maps were computed for every single deflection map, and then used to apply
the process of averaging the signal in radii bins around the galaxies of the lightcone catalog.
Doing so we get the results presented in the following section.

4.3 Numerical correlation computation and results

The whole code that was used for those computations can be found in Appendix B.

22

Chapter 4: Results and interpretation

(a) Log plot of convergence, κ is dimensionless

sionless

(b) Semi-log plot of convergence, κ is dimen-

(c) Log plot of shear, γ is dimensionless

(d) Semi-log plot of shear, γ is dimensionless

(e) Log plot of F flexion, F in 1/arcsec

(f) Semi-log plot of F flexion, F in 1/arcsec

(g) Log plot of G flexion, G in 1/arcsec

(h) Semi-log plot of G flexion, G in 1/arcsec

Figure 4.4: Galaxy-Galaxy lensing correlations for each component and arbitrary lens planes

23

Chapter 4: Results and interpretation

The previous figure stands to show some first approach to the lensing signal, for each
component, convergence, tangential shear, tangential 1-flexion and 3-flexion. With on the
left the loglog plots, and since for some radii bins, the signal gets too close to zero (but since
we don’t want to introduce some bias in our analysis, we get rid of bins where the signal goes
lower than zero, which explains for example the missing points in panel (c), (d), (e) and (f)
of the figure) the error bars get stretched out, so we included the right panels in semi-log for
better visibility of the signal.

4.4 The difference between satellite and central galaxies

Central galaxies are the primary galaxies located at the center of their dark matter halos. They
typically reside in the most massive and densest part of the halo.

• 1-Halo Term: For central galaxies, the 1-halo term dominates at small scales, showing a
strong lensing signal due to the dense dark matter distribution around the central galaxy.
This signal reflects the mass profile of the entire halo, including any substructures like
satellite galaxies within the same halo. We expect the lensing signal for central galaxies
be maximal at small separations from the galaxy and continuously decrease with a
relatively steep slope (compared to the satellites) towards a stable value at large scales,
where the 2-halo term now dominates.

• 2-Halo Term: At larger scales, the 2-halo term takes over, reflecting the contribution
from neighboring halos. Beyond a few hundred kiloparsecs, flattens out and extends,
indicating the influence of nearby halos and the large-scale structure.

On the other hand, satellite galaxies are those that orbit within the halo of a larger central

galaxy. They are not at the center of the halo and typically inhabit the outer regions.

• 1-Halo Term: For satellite galaxies, the 1-halo term still contributes at small scales,
but it reflects the dark matter distribution around the satellite within the larger halo of
the central galaxy. The lensing signal at small scales for satellite galaxies is expected to
decrease with a softer slope compared to central galaxies, because it represents the less
dense outer regions of the halo where satellites reside. The peak of the 1-halo term may
also shift slightly outward compared to central galaxies, because the main halo is thus
in the surroundings of the satellite.

• 2-Halo Term: Similar to central galaxies, the 2-halo term for satellites becomes signif-
icant at larger scales. The larger scale lensing signal will again show a flatter extension,
indicating the influence of surrounding halos. However, this signal might be somewhat
weaker compared to that of central galaxies due to the satellite’s position within the
larger halo.

24

Chapter 4: Results and interpretation

(a) Log plot of convergence, κ is dimensionless

sionless

(b) Semi-log plot of convergence, κ is dimen-

(c) Log plot of shear, γ is dimensionless

(d) Semi-log plot of shear, γ is dimensionless

(e) Log plot of F flexion, F in 1/arcsec

(f) Semi-log plot of F flexion, F in 1/arcsec

(g) Log plot of G flexion, G in 1/arcsec

(h) Semi-log plot of G flexion, G in 1/arcsec

Figure 4.5: Galaxy-Galaxy lensing correlations for each component and arbitrary lens planes

25

Chapter 4: Results and interpretation

We see that indeed, for satellite galaxies, the signal decreases slower than for the central
ones, with some sort of "bump" visible towards the middle scales, due to the effects described
earlier, the satellite galaxies, located in the outer part of the halos, cause signal further away
due to the proximity of the main halo.

4.5 Fitting the SIS and NFW models

4.5.1 Central galaxies

(a) κ lensing signal

(b) γ+ lensing signal

(c) G+ lensing signal

(d) F+ lensing signal

Figure 4.6: Lensing signals for each component, fitted with the SIS and NFW models. Com-

puted by the code in Appendix B.

For information, we fit the models with a certain cutoff radius set arbitrarily, because we want
to restrict the models to the one-halo term. Here, the radius cutoff has been chosen to be
around 800 kpc. We do the same fits for the satellite galaxies.

26

Chapter 4: Results and interpretation

4.5.2 Satellite galaxies

(a) κ lensing signal

(b) γ+ lensing signal

(c) G+ lensing signal

(d) F+ lensing signal

Figure 4.7: Lensing signals for each component, fitted with the SIS and NFW models. Com-

puted by the code in Appendix B.

27

CHAPTER 5

Euclid, SourceExtractor++ and Machine
Learning

As a quick complement before ending this report, this section stands to present the rest of
the work that will be done during the time that lasts between the submission and the end of
this internship. Now that the flexion signal has been studied, measured and understood in the
Horizon-AGN simulated data, we want to apply the measurement process on actual sources
from the Euclid survey. To do so, we will use the SourceExtractor++ software.

5.1 The Euclid Survey

The Euclid satellite is a mission led by the European Space Agency (ESA) with the primary
objective of mapping the geometry of the dark universe. Launched to observe billions of
galaxies over more than a third of the sky, Euclid aims to investigate the nature of dark energy
and dark matter, which together constitute about 95% of the universe’s mass-energy content.
Let us recall that Euclid is equipped with a visible imaging channel (VIS) and a near-infrared
spectrometer and photometer (NISP) thanks to which it will measure the shapes and redshifts
of distant galaxies with unprecedented precision. By doing so, it will provide crucial data for
understanding the expansion history of the universe and the growth of cosmic structures.

Euclid’s data will enable scientists to test theories of gravity on cosmic scales, offering
insights into fundamental questions about the nature of the universe. The mission’s com-
prehensive sky survey and its ability to combine multiple observational techniques make it a
huge project for cosmology and astrophysics in the coming decade.

5.2 SourceExtractor++

SourceExtractor++ (also known as SExtractor++ or SEPP) is a state-of-the-art tool designed
for extracting and analyzing sources such as stars and galaxies from astronomical images.
It enhances the older SourceExtractor (SExtractor) software, providing superior features for
analyzing astronomical data. The software boosts performance through multi-threading and

28

Chapter 5: Euclid, SourceExtractor++ and Machine Learning

parallel processing, making it efficient in handling large datasets. Improved memory manage-
ment further aids in processing large-scale images and complex datasets. SourceExtractor++
features advanced detection algorithms to accurately identify faint and overlapping sources,
while reducing noise and background fluctuations in astronomical images. A key benefit of
SourceExtractor++ is its modularity, which allows users to customize and expand the software
according to their specific requirements. This flexibility is supported by various plug-ins
that increase functionality and adaptability. SourceExtractor++ can output data in multiple
formats, ensuring integration with other software and processing workspaces. The software
produces detailed catalogs of detected sources, covering a wide range of characteristics and
parameters.
It is widely used in significant astronomical surveys to identify and catalog
celestial objects and is applicable to any project requiring source extraction and analysis from
image data. The software delivers highly accurate source extraction and access to photometry
Its high-performance capability allows for the processing of
information on the sources.
large datasets within feasible timeframes, while its modular design and customization options
accommodate diverse research needs. We mentionned that SEPP was an upgrade from SE, via
the use of multithreading for instance, but we are most interested in its compatibility with the
ONNX framework, which will allow us to modify SEPP "from outside", since it is a complicated
software written mainly in C++, we don’t want to do the dirty work of rewriting by hand the
new flexion parameters that need to be fitted.

5.3 ONNX

ONNX, which stands for Open Neural Network Exchange, is an open-source format designed
to represent machine learning models. ONNX provides tools to interchange models between
different deep learning frameworks.
It facilitates model portability and deployment across
various platforms, thereby promoting interoperability in the ecosystem of machine learning
tools. Integrating ONNX with SourceExtractor++ will allow us to enhance the capabilities of
the software, especially in source shape measurement, incorporating flexion terms in SEPP,
without having to modify the source code. Using a deep learning framework like PyTorch
or TensorFlow will allow us to convert a trained model to the ONNX format. Since SEPP is
already compatible with ONNX, the work will lie in our capacity to create the first model and
then validating it with data.

29

CHAPTER 6

Conclusion and Future Work

30

Bibliography

[1] M. Meneghetti, Introduction to Gravitational Lensing With Python Examples, vol. 956.
→ [p4]

2021.

[2] D. J. Bacon, D. M. Goldberg, B. T. P. Rowe, and A. N. Taylor, “Weak gravitational flexion,”
→ [p6], [p7], [p41]

, vol. 365, pp. 414–428, Jan. 2006.

[3] Y. Dubois, C. Pichon, C. Welker, D. Le Borgne, J. Devriendt, C. Laigle, S. Codis,
D. Pogosyan, S. Arnouts, K. Benabed, E. Bertin, J. Blaizot, F. Bouchet, J. F. Cardoso,
S. Colombi, V. de Lapparent, V. Desjacques, R. Gavazzi, S. Kassin, T. Kimm, H. McCracken,
B. Milliard, S. Peirani, S. Prunet, S. Rouberol, J. Silk, A. Slyz, T. Sousbie, R. Teyssier,
L. Tresse, M. Treyer, D. Vibert, and M. Volonteri, “Dancing in the dark: galactic properties
→ [p8]
trace spin swings along the cosmic web,” , vol. 444, pp. 1453–1468, Oct. 2014.

[4] R. Teyssier, “Cosmological hydrodynamics with adaptive mesh refinement. A new high
→ [p8]

resolution code called RAMSES,” , vol. 385, pp. 337–364, Apr. 2002.

[5] E. F. Toro, M. Spruce, and W. Speares, “Restoration of the contact surface in the HLL-
→ [p8]

Riemann solver,” Shock Waves, vol. 4, pp. 25–34, July 1994.

[6] C. Gouin, R. Gavazzi, C. Pichon, Y. Dubois, C. Laigle, N. E. Chisari, S. Codis, J. Devriendt,
and S. Peirani, “Weak lensing in the Horizon-AGN simulation lightcone. Small-scale
→ [p9], [p10], [p41]
baryonic effects,” , vol. 626, p. A72, June 2019.

[7] M. Petkova, R. B. Metcalf, and C. Giocoli, “glamer – II. Multiple-plane gravitational
lensing,” Monthly Notices of the Royal Astronomical Society, vol. 445, pp. 1954–1966, 10
→ [p9]
2014.

[8] J. Binney and S. Tremaine, Galactic dynamics. 1987.

→ [p13]

[9] M. Bartelmann and P. Schneider, “Weak gravitational lensing,” Physics Reports, vol. 340,
→ [p14]

p. 291–472, Jan. 2001.

[10] J. F. Navarro, C. S. Frenk, and S. D. M. White, “The Structure of Cold Dark Matter Halos,”
→ [p14]

, vol. 462, p. 563, May 1996.

[11] B. Jain and U. Seljak, “Cosmological model predictions for weak lensing: Linear and
→ [p17]

nonlinear regimes,” The Astrophysical Journal, vol. 484, p. 560–573, Aug. 1997.

31

CHAPTER 7

Appendix

A Julia code for computing the lensing maps

The following Julia functions are : alpha2jac() used to compute the Jacobian matrix of
the source-image coordinate transformation, and the various jac2...() functions, such as
jac2kappa() that take as input the Jacobian computed in alpha2jac() and transform it into
the lensing maps for each quantity, as well as some extra maps such as the rotational map,
which will allow us to verify that Eq. (3.3) is verified, computed with jac2rot(), we will also
check that our numerical derivatives commute correctly, thanks to jac2crossderivativea1()
and jac2crossderivativea2()

function alpha2jac( alpha::Array{Float32,3} ; scale::Float32=1.0f0 )

nn = size(alpha)
(nn1::Int64,nn2::Int64) = (nn[1],nn[2])
(conv1,conv2) = (nn1/scale, nn2/scale)
jac = Array{Float32,3}(undef,nn1,nn2,12)
d = Array{Float32,1}(undef, nn1)
dd = Array{Float32,1}(undef, nn1)
ddd = Array{Float32,1}(undef, nn1)

## Computing d1a1 and d1d1a1

@inbounds for j=1:nn2

@inbounds @simd for i=1:nn1

d[i] = alpha[i,j,1] * conv1

end
deriv!( d, dd )
deriv!( dd, ddd )
@inbounds @simd for

i=1:nn1

jac[i,j,1] = dd[i]
jac[i,j,5] = ddd[i]

end

end

32

Chapter 7: Appendix

## Computing d1a2 and d1d1a2

@inbounds for

j=1:nn2
@inbounds @simd for i=1:nn1

d[i] = alpha[i,j,2] * conv1

end
deriv!( d, dd )
deriv!( dd, ddd )
@inbounds @simd for i=1:nn1

jac[i,j,2] = dd[i]
jac[i,j,6] = ddd[i]

end

end

d = Array{Float32,1}(undef, nn2)
dd = Array{Float32,1}(undef, nn2)
ddd = Array{Float32,1}(undef, nn2)

## Computing d2a1 and d2d2a1
@inbounds for i=1:nn1

@inbounds @simd for j=1:nn2

d[j] = alpha[i,j,1] * conv2

end
deriv!( d, dd )
deriv!( dd, ddd )
@inbounds @simd for j=1:nn2
jac[i,j,3] = dd[j]
jac[i,j,7] = ddd[j]

end

end

## Computing d2a2 and d2d2a2
@inbounds for i=1:nn1

@inbounds @simd for j=1:nn2

d[j] = alpha[i,j,2] * conv2

end
deriv!( d, dd )
deriv!( dd, ddd )
@inbounds @simd for j=1:nn2
jac[i,j,4] = dd[j]
jac[i,j,8] = ddd[j]

end

end

33

Chapter 7: Appendix

## Computing mixed derivatives

## d2d1a1

@inbounds for i=1:nn1

@inbounds @simd for j=1:nn2
d[j] = jac[i,j,1] /conv1

end
deriv!( d, dd )
@inbounds @simd for j=1:nn2
jac[i,j,9] = dd[j]

end

end

## d2d1a2

@inbounds for i=1:nn1

@inbounds @simd for j=1:nn2
d[j] = jac[i,j,2] /conv1

end
deriv!( d, dd )
@inbounds @simd for j=1:nn2
jac[i,j,10] = dd[j]

end

end

## d1d2a1

@inbounds for j=1:nn2

@inbounds @simd for i=1:nn1

d[i] = jac[i,j,3] /conv2

end
deriv!( d, dd )
@inbounds @simd for i=1:nn1
jac[i,j,11] = dd[j]

end

end

## d1d2a2

@inbounds for j=1:nn2

@inbounds @simd for i=1:nn1

d[i] = jac[i,j,4] /conv2

end
deriv!( d, dd )
@inbounds @simd for i=1:nn1
jac[i,j,12] = dd[j]

end

end

34

Chapter 7: Appendix

return jac

end

function jac2kappa( jac::Array{Float32,3} )

nn = size(jac)
res = Array{Float32,2}( undef, nn[1], nn[2] )
@inbounds for j=1:nn[2]

@inbounds @simd for i=1:nn[1]

res[i,j] = 0.5f0*( jac[i,j,1] + jac[i,j,4] )

end

end
return res

end
##############################################
function jac2gamma1( jac::Array{Float32,3} )

nn = size(jac)
res = Array{Float32,2}( undef, nn[1], nn[2] )
@inbounds for j=1:nn[2]

@inbounds @simd for i=1:nn[1]

res[i,j] = 0.5f0*( jac[i,j,1] - jac[i,j,4] )

end

end
return res

end
##############################################
function jac2gamma2( jac::Array{Float32,3} )

nn = size(jac)
res = Array{Float32,2}( undef, nn[1], nn[2] )
@inbounds for j=1:nn[2]

@inbounds @simd for i=1:nn[1]

res[i,j] = 0.5f0*( jac[i,j,2] + jac[i,j,3] )

end

end
return res

end
##############################################
function jac2F1( jac::Array{Float32,3} )

nn = size(jac)
res = Array{Float32,2}( undef, nn[1], nn[2] )
@inbounds for j=1:nn[2]

@inbounds @simd for i=1:nn[1]

res[i,j] = 0.5f0*(jac[i,j,5]-jac[i,j,12])
+ 0.5f0*(jac[i,j,7]+jac[i,j,10])

35

Chapter 7: Appendix

end

end
return res

end
##############################################
function jac2F2( jac::Array{Float32,3} )

nn = size(jac)
res = Array{Float32,2}( undef, nn[1], nn[2] )
@inbounds for j=1:nn[2]

@inbounds @simd for i=1:nn[1]

res[i,j] = 0.5f0*(jac[i,j,6]+jac[i,j,11])

- 0.5f0*(jac[i,j,9]-jac[i,j,8])

end

end
return res

end
##############################################
function jac2G1( jac::Array{Float32,3} )

nn = size(jac)
res = Array{Float32,2}( undef, nn[1], nn[2] )
@inbounds for j=1:nn[2]

@inbounds @simd for i=1:nn[1]

res[i,j] = 0.5f0*(jac[i,j,5]-jac[i,j,12])

-0.5f0*(jac[i,j,10]+jac[i,j,7])

end

end
return res

end
##############################################
function jac2G2( jac::Array{Float32,3} )

nn = size(jac)
res = Array{Float32,2}( undef, nn[1], nn[2] )
@inbounds for j=1:nn[2]

@inbounds @simd for i=1:nn[1]
res[i,j] = 0.5f0*(jac[i,j,6]+jac[i,j,11])
+ 0.5f0*(jac[i,j,9]-jac[i,j,8])

end

end
return res

end
##############################################
function jac2rot( jac::Array{Float32,3} )

nn = size(jac)

36

Chapter 7: Appendix

res = Array{Float32,2}( undef, nn[1], nn[2] )
@inbounds for j=1:nn[2]

@inbounds @simd for i=1:nn[1]

res[i,j] = 0.5f0*( jac[i,j,2] - jac[i,j,3] )

end

end
return res

end
##############################################
function jac2crossderivativea1( jac::Array{Float32,3} )

nn = size(jac)
res = Array{Float32,2}(undef, nn[1], nn[2] )
@inbounds for j=1:nn[2]

@inbounds @simd for i=1:nn[1]

res[i,j] = 0.5f0*(jac[i,j,11]-jac[i,j,9])

end

end
return res

end
##############################################
function jac2crossderivativea2( jac::Array{Float32,3} )

nn = size(jac)
res = Array{Float32,2}(undef, nn[1], nn[2] )
@inbounds for j=1:nn[2]

@inbounds @simd for i=1:nn[1]

res[i,j] = 0.5f0*(jac[i,j,10]-jac[i,j,12])

end

end
return res

end

B Julia code for computing the galaxy-galaxy lensing

correlation function

function winc( i::Int64, j::Int64, flexf1::Array{Float32,2},
flexf2::Array{Float32,2}, flexg1::Array{Float32,2}, flexg2::Array{Float32,2},
g1::Array{Float32,2}, g2::Array{Float32,2}, kappa::Array{Float32,2},
bb, cd2, acd1, rmin2, lrmin, dlr, nr, w::Vector{Float32}, ir::Vector{Int64},
kp::Vector{Bool} )

w[1] = (j-bb[6])*cd2 #y
w[2] = w[1]^2

#y2

37

Chapter 7: Appendix

#x2

w[3] = (i-bb[5])*acd1 #x
w[4] = w[3]^2
w[5] = w[4] + w[2] #r2
r = sqrt(w[5])
r3 = r^3
x3 = w[3]^3
y3 = w[1]^3
kp[1] = w[5]>rmin2
kp[1] || return
w[6] = ( -(w[4]-w[2])*g1[i,j] - 2.f0*w[1]*w[3]*g2[i,j] ) / w[5]
w[7] = kappa[i,j]
w[8] = (-flexf1[i,j]*w[3]-w[1]*flexf2[i,j])/r
#w[8] = (flexf1[i,j]*w[1]-w[3]*flexf2[i,j])/r
w[9] = -flexg1[i,j]*(4*(x3/r3)-3*w[3]/r)-flexg2[i,j]*(-4*(y3/r3)+3*w[1]/r)
#-Re(G)
#w[9] = -flexg2[i,j]*(4*(x3/r3)-3*w[3]/r)+flexg1[i,j]*(-4*(y3/r3)+3*w[1]/r)
#-Im(G)
ir[1] = ceil( Int64, (0.5f0*log(w[5])-lrmin)/dlr )
kp[1] = (ir[1]>=1) & (ir[1]<=nr)

#-Re(F)

#-Im(F)

end

##############################################
function comp_gs_corr( ra, dec, fn::AbstractString; rf=4,
dx=1.0f0, rmin=0.05f0, rmax=5f0, nr=50 )

ng=length(ra)
@assert length(dec)==ng

######### Deflection and derived information
@time (alpha,zz) = read_bin( fn )
nn=size(alpha)
@info "done reading deflection map!"
#@time jac=alpha2jac( alpha, scale=dx )
#@info "done computing 1st order derivatives!"
#@time jac=jac2der(jac)
@time jac = alpha2jac_flexion_w_cross_terms(alpha,scale=dx)
alpha=nothing ; GC.gc(true)
@info "done computing jacobian matrix map"
@time begin

kappa = jac2kappa(jac)
gamma1 = jac2gamma1(jac)
gamma2 = jac2gamma2(jac)

38

Chapter 7: Appendix

F1 = jac2F1(jac)
F2 = jac2F2(jac)
G1 = jac2G1(jac)
G2 = jac2G2(jac)

end
jac=nothing ; GC.gc(true)
@info "done computing raw convergence, shear and flexion maps"
nn=(nn[1]/rf, nn[2]/rf)
@time begin

kappa=rebin(kappa,factor=rf)
gamma1=rebin(gamma1,factor=rf)
gamma2=rebin(gamma2,factor=rf)
F1=rebin(F1,factor=rf)
F2=rebin(F2,factor=rf)
G1=rebin(G1,factor=rf)
G2=rebin(G2,factor=rf)

end
@info "done computing rebinned convergence & shear maps"
cd1 = -dx/nn[1] ; acd1=abs(cd1)
cr1 = (nn[1]+1)/2f0
cd2 = dx/nn[2]
cr2 = (nn[2]+1)/2f0

######### Rad bins
lrmin=log(rmin/60.f0) ; lrmax=log(rmax/60.f0) ;rmin2=(rmin/60.f0)^2
dlr=(lrmax-lrmin)/(nr)
W=ceil(Int64,rmax/abs(cd2)/60.f0)
rbin = exp.( range( lrmin, stop=lrmax,length=nr+1 ) ) ;

Sres=zeros(5,nr,ng)

#######
Threads.@threads for ig in eachindex(ra)

bb=get_bbox2( ra[ig], dec[ig], W, cr1, cr2, cd1, cd2, nn )
#afficher bb 1,2,3,4
w=zeros(Float32,9)
lres=zeros(Float64,5,nr)
ir=[1]
kp=[false]

@inbounds for j=bb[3]:bb[4]

@inbounds for i=bb[1]:bb[2]

winc( i, j, F1, F2, G1, G2, gamma1, gamma2, kappa,

39

Chapter 7: Appendix

bb, cd2, acd1, rmin2, lrmin, dlr, nr, w, ir, kp )
if kp[1]

lres[1,ir[1]] += Float64(w[9]) #storing G+
lres[2,ir[1]] += Float64(w[8]) #storing F+
lres[3,ir[1]] += Float64(w[7]) #storing kappa
lres[4,ir[1]] += Float64(w[6]) #storing gamma+
lres[5,ir[1]] += 1.0

#counting pixels

end

end

end
Sres[:,:,ig] .= lres

end
#
#
return Sres, rbin, kappa

Sres = pmap( ig->loc_calc(ra[ig],dec[ig]) , 1:ng )
res = sdata(Sres)

end

40

.

.

.

.

.

.

5

7

9

11

12

12

List of Figures

2.1 Definition of lensing angles and distances. Figure by Michael SACHS .
2.2

.
Effects of the different lensing fields on a Gaussian galaxy of radius 1 arcsec.
10% convergence/shear and 0.28 arcsec−1 flexion (which is a very large value
for this quantity, chosen only for visualizing) are applied. Figure extracted
.
.
from [2] .

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

3.2

3.1 Two dimensional graph of the past lightcone generated with Horizon-AGN.
Each blue cell is a replicate of the H-AGN simulation box. The orange vertical
line marks z = 1 and the tiling repeats up to z = 8. The black lines represents
the aperture angle of the cone which goes from 2.25 arcseconds to 1.00 arcsec-
ond after z = 1. Figure extracted from [6] .
.
Expressions as derivatives of the two components of the deflection fields of
.
the seven lensing quantities defined earlier and used .
3.3 Tangential projection of the lensing quantities, generated by gravitational
.
lensing, used to effectively compute and measure the lensing quantities.
3.4 Cross projection of the lensing quantities, which should yield no signal and
.
3.5 Tangential projection of the lensing quantities, expressed in flat sky plane
.
.

thus can be used as a robust way of checking for systematical errors.

coordinate.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

12
.
3.6 Cross projection of the lensing quantities, expressed in flat sky plane coordinates. 13
3.7 The process of measuring the galaxy-galaxy lensing correlation statistic. Fig-
.
.

ure by Masamune Oguri

18

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.
.

.
.

.
.

.
.

.
.

Appendix A.

of derivatives of the two deflection field maps.
.

4.1 Convergence map for a redshift zplane ≈ 1.016, computed from a combination
.
.
4.2 The two different shear fields γ1 and γ2
4.3 The four different flexion fields F1, F2 and G1, G2. Computed by the code in
.
.
.

22
.
4.4 Galaxy-Galaxy lensing correlations for each component and arbitrary lens planes 23
4.5 Galaxy-Galaxy lensing correlations for each component and arbitrary lens planes 25
Lensing signals for each component, fitted with the SIS and NFW models.
4.6
Computed by the code in Appendix B. .
.
Lensing signals for each component, fitted with the SIS and NFW models.
.
Computed by the code in Appendix B. .

.
.

21
21

4.7

27

26

.
.

.
.

.
.

.
.

.
.

.
.

.
.

.
.

.
.

.
.

.
.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

41

List of Tables

3.1

List of parameters used for the flat ΛCDM cosmology defined in the Horizon-
AGN simulation and their values, all parameters are unitless, except for the
.
Hubble parameter H0, expressed in km s−1Mpc−1 .
3.2 Available information in the header of the galaxy catalog within the Horizon-
.
.

AGN simulated light cone.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

.

8

16

42

List of Listings

43

