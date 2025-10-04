DistributeLib
=============

A small library that contains:

- A number of 32-bit sequence generators -- various classic RNGs, plus a set of
  stratified progressive sequences.

- Functions for converting 32-bit integer inputs from such generators into
  various 1D, 2D, and 3D distributions and types. The intended use is to convert
  the raw output of any generator into various useful distributions.

Features:

* Generation of int, float, 2D and 3D samples over various shapes
* Cheap routines that only take one input number, but may be less accurate or
  controllable
* More accurate routines that take a full number of input samples, e.g., three
  for a cube
* Pre-modulation of inputs into triangular, gaussian-like, and weighted forms
* Gaussian (normal) distribution as both a cheap approximation and via full
  Box-Muller transform
* Generators supplying LCG/PCG/XorShift/Hash/Halton/Sobol/Golden/'R' sequences.

To build and run the test app:

	make
    ./distribute 1000  # also see distribute.svg after running

Or add Generate.* and Distribute.* to your favourite IDE.

There is an [online demo](https://andrewwillmott.github.io/app/GenDistTest.html)
of most of the included generators, modifiers, and distribution functions.

Examples
--------

Assuming 'rng' returns the next sample via operator uint32_t(), usage can be as
simple as that below. Otherwise, substitute MyRNG() or whatever other syntax is
appropriate for your sample generator.

1D:

	float score     = ToFloat(rng, 1.0f, 100.0f);
	int   modifier  = ToInt32Signed(rng, 5);
	int   dayOfYear = ToInt32Inclusive(rng, 1, 365);
	float weightKG  = ToFloat(ModGaussLike(rng), 50.0f, 130.0f);

2D/3D:

    Vec2f discLoc       = ToDisc(generator);
	Vec2f pixTentSample = ToSquare(ModTriangle(rng), ModTriangle(rng));
	Vec3f rayDir        = ToDirection3(rng);

Warning: do not use rand() to feed these functions, as, in addition to its low
quality, its range is not guaranteed to be the full 32 bits.


Output
------

* 1D: float [-1, 1], integer [-10, 10], with a simple LCG input.

	![](images/float.png "Float [-1, 1]")
	![](images/int.png "Integer [-10, 10]")

* 2D: square, ring, triangle

	![](images/square.png "square")
	![](images/ring.png "ring")
	![](images/triangle.png "triangle")


* Triangle Modifier: integer, float, square in both dimensions

	![](images/int-triangle.png "Integer triangular distribution")
	![](images/float-triangle.png "Float triangular distribution")
	![](images/square-triangle.png "Square triangular distribution")

* Gauss-like Modifier: integer, square, and full gaussian for comparison.

	![](images/int-gauss-like.png "Gauss-like integer -10 -- 10")
	![](images/square-gauss-like.png "Gauss-like over the square")
	![](images/gaussian.png "Full Gaussian")

* Weighted Modifier: float and circle weighted by [1, 8, 0, 4, 1]

	![](images/float-weighted.png "Float weighted by 1, 8, 0, 4, 1")
	![](images/circle-radially-weighted.png "Circle radially weighted")

* Non-random inputs: Halton 100, 1000, 5000 samples

	![](images/halton-100.png "Halton 2D 100 points")
	![](images/halton-1000.png "Halton 2D 1000 points")
	![](images/halton-5000.png "Halton 2D 5000 points")

* Non-random inputs: Golden ratio square, circle, and sphere surface

	![](images/square-golden.png "Golden Square")
	![](images/circle-golden.png "Golden Circle")
	![](images/sphere-surface-golden.png "Golden Sphere Surface")

* Progressively adding samples from the PCG, Halton, and Golden sequences

    <a href="images/circle-pcg-anim.gif"   ><img src="images/circle-pcg-anim.gif"    alt="PCG Circle"             width="256"/></a>
    <a href="images/circle-halton-anim.gif"><img src="images/circle-halton-anim.gif" alt="Halton Circle"          width="256"/></a>
    <a href="images/circle-golden-anim.gif"><img src="images/circle-golden-anim.gif" alt="Golden (Spiral) Circle" width="256"/></a>

* Progressively adding samples from the Rd sequence: square, circle, cube

    <a href="images/square-rd-anim.gif"><img src="images/square-rd-anim.gif" alt="Rd Square" width="256"/></a>
    <a href="images/circle-rd-anim.gif"><img src="images/circle-rd-anim.gif" alt="Rd Circle" width="256"/></a>
    <a href="images/cube-rd-anim.gif"  ><img src="images/cube-rd-anim.gif"   alt="Rd Cube"   width="256"/></a>
