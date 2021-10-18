import Data.Vector.Unboxed hiding (foldr1, foldr, force, sum, length, last, zipWith, map, (++))
import qualified Data.Vector.Unboxed as V
import qualified Data.Vector as VV
import Data.Maybe
import Data.List
import Data.Ord
import Data.Either
import Numeric.IEEE
import System.Random
import System.Random.Stateful
import Control.Exception
import Debug.Trace
import System.Environment
import System.Random.MWC as MWC
import System.Random.MWC.Distributions as MWC
import Coords as C

-- |A force is something that causes a change to a distance. The
-- |term is used pretty loosely.
type Force = Double

-- |What is an appropriate force power law for a multi-dimensional world?
-- |In our universe, first-order forces with zero-mass bosons follow a 
-- |one-over-r-squared law, but this is because our universe is 3d.
-- |
-- |Find a power-law function, given a set of coords. This initial version
-- |assumes the same power law for every distance. Arguably, the cutoff when
-- |counting dimensions depends on the distance scale, so we could try a power
-- |law that has more dimensions for smaller distances.
power_law :: VV.Vector (Vector C.Coord) -> (Distance -> Double)
power_law c = let count = count_coord_dimensions c in
    (\s -> 1.0 / (s^(count - 1))) -- count - 1, because we differentiate

-- |Given a vector of coordinates, return a vector of forces applied to
-- |each distance between two points
forces :: VV.Vector (Vector C.Distance) -> VV.Vector (Vector C.Coord) -> VV.Vector (Vector Force)
forces c = forces_given_power_law (power_law c) c

-- |Given a vector of coordinates and a power law,
-- |return a vector of forces applied to
-- |each distance between two points
forces_given_power_law :: (Distance -> Double) -> 
    VV.Vector (Vector C.Coord) -> VV.Vector (Vector C.Distance) -> 
    VV.Vector (Vector Force)
forces_given_power_law f c s = VV.generate (VV.length c) generator where
    generator i = V.generate (V.length (c VV.! i)) (generator' i)
    generator' i = force f c s i

-- |Find the force between two points i and j in a vector of coordinates.
-- |We need to add up the dot products of each vector ik with the vector ij
force :: (Distance -> Double) -> 
    VV.Vector (Vector C.Coord) -> VV.Vector (Vector C.Distance) ->
    Int -> Int -> Force
force f c s i j = foldr (add_both_components f c s i j) 0.0 [0..n] where
    n = (VV.length c) - 1

-- |Adds two vector component to the overall force reducing the distance
-- |from point i to j. Adds the effect of the point k relative to i, and
-- |relative to j.
add_both_components :: (Distance -> Double) -> 
    VV.Vector (Vector C.Coord) -> VV.Vector (Vector C.Distance) ->
    Int -> Int -> Int -> Force -> Force
add_both_components f c s i j k acc =
    -- trace ("add_both_components i=" ++ show i ++ " j=" ++ show j ++ " k=" ++ show k)
    add_component f c s j i k $ add_component f c s i j k acc

-- |Adds one vector component to the overall force in the direction of the distance
-- |from point i to j. The component added is due to point k relative to point i
add_component :: (Distance -> Double) -> 
    VV.Vector (Vector C.Coord) -> VV.Vector (Vector C.Distance) ->
    Int -> Int -> Int -> Force -> Force
add_component f c s i j k acc = 
    if i == k then 
        -- trace ("add_component i=" ++ show i ++ " j=" ++ show j ++ " k=" ++ show k ++ " skip")
        acc     -- a point exerts no force on itself
    else
        let
            dot = dot_product c i j k
            dist_ij = distance s i j
            dist_ik = distance s i k
            factor = f dist_ik
        in
            -- trace ("add_component i=" ++ show i ++ " j=" ++ show j ++ " k=" ++ show k ++ " dot=" ++ show dot)
            -- trace ("    dist_ij=" ++ show dist_ij ++ " dist_ik=" ++ show dist_ik ++ " factor" ++ show factor)
            acc + dot * factor / (dist_ij * dist_ik)

-- |Calculate the dot product of vector ik and vector ij within coordinate set c_ij
dot_product :: VV.Vector (Vector Double) -> Int -> Int -> Int -> Double
dot_product c i j k = 
    let
        ci = c VV.! i
        cj = c VV.! j
        ck = c VV.! k
    in
        V.sum $ izipWith (\l x y -> let z = (get_or 0.0 ci l) in (x - z) * (y - z)) cj ck

-- |Fetch the distance between two points i and j within a distance set s_ij
distance :: VV.Vector (Vector C.Distance) -> Int -> Int -> Distance
distance s i j =
    if i > j then
        (s VV.! (i-1)) ! j
    else
        (s VV.! (j-1)) ! i

-- |Gets a value from a vector or a default if off the end
get_or :: Unbox a => a -> Vector a -> Int -> a
get_or def c i = case c!?i of
    Just x -> x
    Nothing -> def

-- |Tests the get_or method, before and after the end of the vector
test_get_or :: (Int, Int)
test_get_or = 
    let
        v = V.fromList [1, 2]
        a = get_or 0 v 1
        b = get_or 0 v 2
    in
        assert (a == 2 && b == 0) 
        (a, b)     

-- |Tests the distance method, with the end points in either order
test_distances :: (Distance, Distance, Distance)
test_distances =
    let
        s = fromLists [[1.0], [1.1, 2.1], [1.2, 2.2, 3.2]]
        s01 = distance s 0 1
        s13 = distance s 1 3
        s30 = distance s 3 0
    in
        assert (s01 == 1.0 && s13 == 2.2 && s30 == 1.2)
        (s01, s13, s30)

-- |Tests the dot product function
test_dot_product :: Double
test_dot_product =
    let
        c = fromLists [[], [1.0], [1.1, 2.1], [1.2, 2.2, 3.2]]
        (i, j, k) = (1, 2, 3)
        -- vector ij = (1.1 - 1.0, 2.1 - 0.0) = (0.1, 2.1) 
        -- vector ik = (1.2 - 1.0, 2.2 - 0.0, 3.2 - 0.0) = (0.2, 2.2, 3.2)
        -- ij.ik = 0.1 * 0.2 + 2.1 * 2.2 = 4.64
        result = dot_product c i j k
    in
        assert (C.approx_equal 1e-6 result 4.64)
        result

-- |Tests the force function, assuming a simple factor function
-- |that always returns 1.0 operating on a simplex of points
test_force :: Force
test_force =
    let
        dists = fromLists [
            [1.0],
            [1.0, 1.0],
            [1.0, 1.0, 1.0]]
        coords = fromLists [
            [],
            [1.0],
            [0.5, 0.8660254],
            [0.5,0.28867513,0.8164966]]
        result = force (\_ -> 1.0) coords dists 1 3
    in
        assert (C.approx_equal 1e-6 result 1.75)
        result

-- |Tests the forces_given_power_law function, assuming a simple factor function
-- |that always returns 1.0 operating on a simplex of points
test_forces_given_power_law :: VV.Vector (Vector Force)
test_forces_given_power_law =
    let
        dists = fromLists [
            [1.0],
            [1.0, 1.0],
            [1.0, 1.0, 1.0]]
        coords = fromLists [
            [],
            [1.0],
            [0.5, 0.8660254],
            [0.5,0.28867513,0.8164966]]
        expected = fromLists [
            [],
            [2.0],
            [2.0, 1.75],
            [2.0, 1.75, 1.33333333333]]
        result = forces_given_power_law (\_ -> 1.0) coords dists 
    in
        assert (C.approx_equal_vv 1e-6 result expected)
        result

-- |Tests the forces function on a random simplex of dimension n
test_forces_simplex :: StatefulGen g m => g -> Int -> m (Either String (VV.Vector (Vector Force)))
test_forces_simplex rnd n = do
    simplex <- C.create_random_triangular_matrix rnd (0.99, 1.01) n
    return (test_forces simplex)

-- |Tests the forces function on a given set of distances
test_forces :: VV.Vector (Vector Distance) -> Either String (VV.Vector (Vector Force))
test_forces s = do
    c <- C.validated_coords s
    return (forces c s)

main :: IO ()
main = do

    putStrLn ""
    putStrLn "Running Forces tests..."

    putStrLn $ "test_get_or: " ++ (show test_get_or)
    putStrLn $ "test_distances: " ++ (show test_distances)
    putStrLn $ "test_dot_product: " ++ (show test_dot_product)
    putStrLn $ "test_force: " ++ (show test_force)
    putStrLn $ "test_forces_given_power_law: " ++ (show test_forces_given_power_law)

    -- One generator shared by all tests that need random numbers
    rnd <- createSystemRandom
    forces_simplex <- test_forces_simplex rnd 10
    putStrLn $ "test_forces_simplex: " ++ (show forces_simplex)
        
    putStrLn "...all tests passed"
    putStrLn ""
