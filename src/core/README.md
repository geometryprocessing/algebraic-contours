# Core

Core geometric methods and data structures. These are simple and general geometric objects that are used for downstream spline construction and contour computation. Methods act on general data structures representing some geometric object, and data structures encapsulate more complex data types.

Included geometric objects and methods:

- `AbstractCurveNetwork`: Abstract network of transversally intersecting curves.
- `AffineManifold`: Triangulated manifold with affine transition maps between charts.
- `Bivariate Quadratic Functions`: Methods to operate on bivariate quadratics represented by coefficient vectors.
- `Conic`: Rational function of degree 2
- `ConvexPolygon`: Convex polygon formed by intersecting half planes.
- `Halfedge`: Boilerplate halfedge mesh representation.
- `Interval`: Connected open, closed, bounded, and unbounded intervals.
- `LineSegment`: Conic of degree 1
- `Polynomial Functions`: Methods to operate on polynomials represented by coefficient vectors.
- `RationalFunction`: Quotient of scalar or vector valued polynomial functions over an interval.
