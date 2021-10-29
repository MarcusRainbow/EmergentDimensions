# Emergent Dimensions

 Exploring the emergence of finite dimensionality in models such as Quantum Loop Gravity.

 We start with a toy universe, made up of points that are not specified in some cartesian geometry (that would be begging the question). Instead, the points are specified purely in terms of the distances between them.

 Assuming Euclidean geometry (or near Euclidean, anyway), we count the minimum number of dimensions required to reproduce the distances between the points.

 Normally, I'd expect the number of dimensions to be one less than the number of points, unless we are lucky with the distances.

 We now experiment with ways to evolve the points. Is there a natural evolution or dynamic that results in a small finite number of dimensions?

## A Creation Myth

 In the beginning, everything was without form and void. There were no spatial dimensions, but a very large number of points separated by widely varying distances. In fact, it would not be possible to define spatial dimensions to accomodate the points, because the distances between them was inconsistent -- for example, it might be a shorter distance to go from one point to another via a third point than to go direct.

 Within this chaos, there was a force. Maybe it was attractive, maybe repulsive -- let's call it inflation for now. In later eras of just three spatial dimensions, the force might follow an inverse square law, but for now, it is extremely sensitive to distance, because if the world is dimensional at all, it has a very large number of dimensions.

 This early universe was rather like a many dimensional onion, with each layer of the onion being made up of points roughly the same distance from each other. The effect of the force was to blow the layers apart, and as they grew further apart the difference in force grew without limit.

 At the end of the first day, all that was left in our universe was a collection of points the same distance apart -- mathematically we call this a _simplex_. There was an onion shell larger than ours, but so much larger and slower that its smallest division of time would be greater than an eternity for us. There was once an onion shell smaller than ours, but it had gone through its entire lifecycle. Galaxies had appeared. On a few of them, intelligent civilisations had come and gone. The galaxies had moved further apart and grown colder until there was nothing left but black holes which had evaporated. And all in a flicker of an instant for the universe that is our onion shell. And the onion shells inside and outside of these are inconceivable to us.

 The second day was a day of rest. The whole of our universe was made up of points that were exactly the same distance apart, and defining a spatial universe with as many dimensions as points (actually one less). Nothing happened for a long time. (Okay, a long time compared with the first day, but a vanishingly small time by our standards today.)

 On the third day, tiny differences between the distances appeared. That is where our simulation appears, with a universe that is initially very nearly an equilateral simplex but with tiny differences in momenta or distances.

## Some Implementation Details

## Experiments

### Creating a consistent dimensional universe

One of the first problems I hit was creating a universe that is consistent with a dimensional view of the world. There are very many ways for the universe to be inconsistent. For example, in one dimension we require for every three points A, B and C that AB + BC >= AC. There are similar though more complicated rules for any number of dimensions.

One way to solve the problem is to initially create a universe within a dimensional space. For example you can use N cartesian coordinates and randomly pick points within this space. However, this feels like begging the question. We are interested in the emergence of dimensionality, so we should not start with an assumption of dimensionality.

A method that I found works is to create a number of points that form a roughly equilateral simplex. An equilateral simplex has the same distance between every pair of points. For example, a two dimensional simplex is an equilateral triangle and a three dimensional one is a perfect tetrahedron. The distances can be randomised slightly without causing the failures described in the previous paragraph. A near equilateral simplex can be justified as the result of a creation event, or by the arguments put forward in the creation myth above.

### Verifying the dimensionality of the universe

My method is to take a matrix of the distances between points, which must be symmetric and have a zero main diagonal, as every point has zero distance from itself. I take the first point to be the origin, and the second point to be along the first axis (x axis).

The next point is defined by its distance from the first two points. It is possible that it lies along the x axis. This can be calculated by the two simultaneous equations defined by the pythagorean distances from the first two points, and assuming a coordinate along some new y axis. The triangular nature of the simultaneous equations makes this really easy to do, as the first equation is of the form:

(x - 0)^2 + (y - 0)^2 = s0^2

and the second equation is of the form:

(x - x0)^2 + (y - 0)^2 = s1^2

And so we continue. Adding each new point means working through the existing axes to find the coordinates on those axes, then adding a new axis, which may have a zero coordinate.

This algorith could be improved by a method similar to pivoting when solving a set of linear simultaneous equations. Rather than simply iterating through the points in the supplied order, we could choose an ordering that prioritises the best-defined dimensions.

Dimensions with an extent less than some absolute minimum are ignored, so that numerical noise does not lead to spurious extra dimensions.

### Dynamical evolution of the early universe

Again, to avoid begging the question, I wanted to express the state of the universe in terms of distances and rates of change of distances, rather than assuming that which we wish to prove -- a rectangular coordinate system.

We want to calculate forces aligned along the distances between pairs of points. If we assume some gravitational or dark-energy/inflation force, the force directly between two points is aligned between them and is some function of the distance between the points. However, there will generally be an indirect force on the distance that is caused by other points in the vicinity. This can be calculated as the dot product of two vectors: from the point experiencing the force towards the point causing the force; and the vector aligned along the distance being affected by the force.

Calculating the dot product requires us to use a coordinate system -- in my experiments I used the rectangular cartesian coordinates calculated when validating the dimensionality of the system. Perhaps this is begging the question, but it is only done to calculate the force, and some assumption about dimensionality must be made in order to calculate the dot product.

Another assumption that must be made is the nature of the force acting on any pair of points. I made the assumption that the direction of the force was in line with the two points, and that the magnitude is inversely proportional to some power of the distance. In our 3d world, we are used to inverse square law forces, but this is a direct consequence of the three dimensionality. In my simulation, I assumed a power equal to the number of dimensions minus one, where the number of dimensions is calculated as above.

I assume a simple newtonian mechanics, so each point moves according to:

dv_ij = g F_ij dt
ds_ij = v_ij dt

Where s_ij is the distance between i and j, and F_ij is the force acting on that distance.

Although this formula makes sense to first order, it ignores the second order effect that distances cannot change at a constant rate and independently from other distances. We end up needing to add a correction of the form:

dv_ij/dt = - v_ij^2 (V_ij . S_ij) / s_ij^3

Where v_ij and s_ij are as before, and V_ij and S_ij are the corresponding vectors.

The good news is that indeed a multidimensional system remains multidimensional, rather than ceasing to be dimensional at all. The bad news, in my tests at least, is that dimensions do not disappear. An N dimensional system continues to be N dimensional. The only exception to this is when points get so close to each other that they and their dimensions cease to be important. However, this is not what we want. We do not want a three dimensional universe with only four points in it.

### A better test

How can we improve all this? Let's consider a thought experiment with a really simple world, with only five points. Four points define three dimensions (except for some very special cases, which we ignore here). Let us add a fifth point, and for the sake of argument, make it equidistant from the other four points. There is one distance that works in terms of making that fifth point fit into the same three dimensions as the others. If we make that distance shorter, the fifth point cannot fit into any real universe, of any dimensionality. Making the distance longer means that we need a fourth spatial dimension.

For three dimensions to emerge, we need either of two things to happen: either we need a force that acts to shrink these extra distances so that all points tend to three dimensions; or we need an explosion of scale in just three directions, so that the remaining dimensions are insignificant in scale.

Moreover, this force or scale-explosion must be global. Imagine a universe which is locally three dimensional but where the three dimensional areas are differently rotated relative to each other in some further dimension or dimensions. This is not a description of the universe that we know.
 