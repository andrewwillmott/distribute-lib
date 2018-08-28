DistributeLib
=============

A single-source-file library for converting 32-bit integer inputs into various 1D, 2D, and 
3D distributions and types. The intended use is to convert the output of an RNG into various 
more useful forms.

Features:

* Generation of int, float, 2D and 3D samples over various shapes
* Cheap routines that only take one input number, but may be less accurate or controllable
* More accurate routines that take a full number of input samples, e.g., three for a cube
* Pre-modulation of inputs into triangular, gaussian-like, and weighted forms
* Gaussian (normal) distribution as both a cheap approximation and via full Box-Muller transform


To build and run the test app:

    c++ --std=c++11 Distribute.cpp DistributeTest.cpp -o distribute && ./distribute 200


Examples
--------

Assuming 'rng' returns the next sample via operator uint32_t(), usage can be as simple as 
below, otherwise substitute rng() or rng.next() or whatever else is appropriate.

1D:

	float score     = ToFloat(rng, 1.0f, 100.0f);
	int   modifier  = ToInt32Signed(rng, 5);
	int   dayOfYear = ToInt32Inclusive(rng, 1, 365);
	float weightKG  = ToFloat(ToGaussLike(rng), 50.0f, 130.0f);

2D/3D:

	Vec2f circleLoc     = ToCircle(rng);
	Vec2f pixTentSample = ToSquare(ModTriangle(rng), ModTriangle(rng));
	Vec3f rayDir        = ToDirection3(rng);


Warning: do not simply use rand() to feed these functions, as, in addition to its low 
quality, its range is not guaranteed to be a full 32 bits. Consider a simple LCG or 
something like http://www.pcg-random.org instead.


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

* Gauss-like Modifier: integer, square, and full gaussian for comparison. (Note long tail.)

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
