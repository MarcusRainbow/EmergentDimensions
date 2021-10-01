import Data.Vector.Unboxed hiding (foldr1, sum, length, last, zipWith, map, (++))
import qualified Data.Vector.Unboxed as V
import qualified Data.Vector as VV
import Control.Exception
import Debug.Trace

test_vector_from_list :: Vector Float
test_vector_from_list = fromList [1.0, 2.0, 3.0]

-- |Given a set of distances between points, calculate
-- |a consistent set of coordinates for each point. The
-- |distances are supplied as a vector of vectors as
-- |an encoding of a symmetric matrix with diagonal zero.
-- |The coordinates are returned as vectors where zero
-- |coordinates may be zero-suppressed on the right hand
-- |side. Thus the origin is empty vector. This allows
-- |the implementation to be agnostic of the number of
-- |dimensions.
coords :: VV.Vector (Vector Float) -> VV.Vector (Vector Float)
coords s = let n = (VV.length s + 1) in
    VV.constructN n (coords_constructor s)

-- |Construct one of the coordinate sets x_j given the distances
-- |s_ij and the previous coordinate sets p_ij
coords_constructor ::  VV.Vector (Vector Float) ->  VV.Vector (Vector Float) ->  Vector Float
coords_constructor s p = let n = VV.length p in
    if n == 0 then
        empty  -- first coordinate is always the origin (empty vector)
    else
        calc_coord (s VV.! (n - 1)) p

-- |Similar to Gram-Schmidt orthogonalisation but where
-- |we supply a vector of distances. This call constructs
-- |the cartesian coordinates of one point defined by
-- |its distances from the preceding points, vector s_i, and
-- |the cartesian coordinates of the preceding points p_ij
calc_coord :: Vector Float -> VV.Vector (Vector Float) -> Vector Float
calc_coord s p = case (length p) of
    0 -> empty           -- first coord is the origin (empty vector)
    1 -> fromList [s!0]  -- second coord is (s) where s is distance 
    otherwise -> constructN (V.length (VV.last p) + 1) (coord_constructor s p)

-- |Construct one of the coordinate dimensions x_n given the
-- |preceding dimensions x_i and the preceding coordinates p_ij
coord_constructor :: Vector Float -> VV.Vector (Vector Float) -> Vector Float -> Float
coord_constructor s p x = 
    let 
        -- b = trace ("coord_constructor: s=" ++ (show s) ++ " p=" ++ (show p) ++ " x=" ++ (show x)) (VV.last p)
        b = VV.last p
        n = V.length b
    in
        if (V.length x) < n then
            (s!0^2 - s!n^2 + sum2 b - 2 * (sum_pair b x)) / (2 * b!(n-1))
        else
            sqrt (s!n^2 - sum2_diff x b)

-- |Sum of the squares of the items in a vector
sum2 :: Vector Float -> Float
sum2 = V.foldr (\ x y -> x^2 + y) 0.0

-- |Sum of pairwise products of members in two vectors
sum_pair :: Vector Float -> Vector Float -> Float
-- todo is there a way of doing this without an intermediate vector?
sum_pair v1 v2 = V.sum $ V.zipWith (*) v1 v2

-- |Sum of squares of differences of members in two vectors
sum2_diff :: Vector Float -> Vector Float -> Float
-- todo is there a way of doing this without an intermediate vector?
sum2_diff v1 v2 = sum2 $ V.zipWith (-) v1 v2

-- |Test sum of squares
test_sum2 :: Float
test_sum2 = let
    result = sum2 (fromList [1.0, 2.0, 3.0]) in
    assert (result == 14.0) result

-- |Test sum of pairwise products
test_sum_pair :: Float
test_sum_pair = let 
    result = sum_pair (fromList [1.0, 2.0, 3.0]) (fromList [10.0, 100.0, 1000.0]) in
    assert (result == 3210.0) result

-- |Test sum of pairwise products where vectors are different lengths
test_sum_pair_mismatch :: Float
test_sum_pair_mismatch = let 
    result = sum_pair (fromList [1.0, 2.0, 3.0]) (fromList [10.0, 100.0]) in
    assert (result == 210.0) result

-- |Test sum of pairwise products where vectors are different lengths
test_sum_pair_mismatch2 :: Float
test_sum_pair_mismatch2 = let 
    result = sum_pair (fromList [1.0, 2.0]) (fromList [10.0, 100.0, 1000.0]) in
    assert (result == 210.0) result

-- |Test sum of squares of diffs
test_sum2_diff :: Float
test_sum2_diff = let 
    result = sum2_diff (fromList [11.0, 22.0, 33.0]) (fromList [10.0, 20.0, 30.0]) in
    assert (result == 14.0) result

-- |Test sum of squares of diffs where vectors are different lengths
test_sum2_diff_mismatch :: Float
test_sum2_diff_mismatch = let 
    result = sum2_diff (fromList [11.0, 22.0]) (fromList [10.0, 20.0, 30.0]) in
    assert (result == 5.0) result

-- |Test coord constructor for intermediate step
test_coord_constructor :: Float
test_coord_constructor = let
    root_2 = sqrt 2.0
    sample_dists = fromList [1.0, root_2, root_2, root_2]
    sample_prev = VV.fromList [empty, fromList [1.0], fromList [0.0, 1.0]]
    sample_coords = fromList [0.0]
    result = coord_constructor sample_dists sample_prev sample_coords in
    assert (approx_equal 1e-6 result 0.0) result
    
-- |Test coord constructor for final step
test_coord_constructor_final :: Float
test_coord_constructor_final = let
    root_2 = sqrt 2.0
    sample_dists = fromList [1.0, root_2, root_2, root_2]
    sample_prev = VV.fromList [empty, fromList [1.0], fromList [0.0, 1.0]]
    sample_coords = fromList [0.0, 0.0]
    result = coord_constructor sample_dists sample_prev sample_coords in
    assert (approx_equal 1e-6 result 1.0) result
    
-- |Test coord constructor in case that gave out of bounds error
test_coord_constructor_problem :: Float
test_coord_constructor_problem = let
    sample_dists = fromList [1.0, 1.0, 1.0]
    sample_prev = VV.fromList [empty, fromList [1.0]]
    sample_coords = fromList [0.5]
    result = coord_constructor sample_dists sample_prev sample_coords in
    -- assert (approx_equal 1e-6 result 1.0) 
    result

-- |Test coord generator for a partial hypercube (missing far corners)
test_calc_coord_hypercube :: Vector Float
test_calc_coord_hypercube = let
    root_2 = sqrt 2.0
    sample_dists = fromList [1.0, root_2, root_2, root_2]
    sample_prev = VV.fromList [empty, fromList [1.0], fromList [0.0, 1.0]]
    result = calc_coord sample_dists sample_prev in
    assert (approx_equal_vector 1e-6 result (fromList [0.0, 0.0, 1.0]))
    result

-- |Test coord generator for a 3d simplex
test_calc_coord_simplex :: Vector Float
test_calc_coord_simplex = let
    sample_dists = fromList [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    sample_prev = VV.fromList [empty, fromList [1.0]]
    result = calc_coord sample_dists sample_prev in
    assert (approx_equal_vector 1e-6 result (fromList [0.5, 0.5 * sqrt 3]))
    result

-- |Test coord generator for the first point
test_calc_coord_first :: Vector Float
test_calc_coord_first = let
    sample_dists = fromList [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    sample_prev = VV.empty
    result = calc_coord sample_dists sample_prev in
    assert (V.null result)
    result

-- |Test coord generator for the second point
test_calc_coord_second :: Vector Float
test_calc_coord_second = let
    sample_dists = fromList [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    sample_prev = VV.fromList [empty]
    result = calc_coord sample_dists sample_prev in
    assert (approx_equal_vector 1e-6 result (fromList [1.0]))
    result

-- |Overall test, given a collection of 3 points all one apart
test_coords_simplex :: VV.Vector (Vector Float)
test_coords_simplex = let
    sample_dists = fromLists [
        [1.0],
        [1.0, 1.0],
        [1.0, 1.0, 1.0]]
    result = coords sample_dists in
    assert (approx_equal_vv 1e-6 result (fromLists [
        [],
        [1.0],
        [0.5, 0.8660254],
        [0.57735026,0.24401695,0.7791806]]))
    result

-- |Overall test, given a collection of 5 points all one apart
test_coords_4d_simplex :: VV.Vector (Vector Float)
test_coords_4d_simplex = let
    sample_dists = fromLists [
        [1.0],
        [1.0, 1.0],
        [1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0]]
    result = coords sample_dists in
    assert (approx_equal_vv 1e-6 result (fromLists [
        [],
        [1.0],
        [0.5, 0.8660254],
        [0.57735026,0.24401695,0.7791806],
        [0.64169973,0.16621879,0.11416383,0.7399725]]))
    result

-- |Test of a set of points that form a 3, 4, 5 triangle
test_coords_345 :: VV.Vector (Vector Float)
test_coords_345 = let
    sample_dists = fromLists [
        [3.0],
        [4.0, 5.0]]
    result = coords sample_dists in
    assert (approx_equal_vv 1e-6 result (fromLists [
        [],
        [3.0],
        [0.0, 4.0]]))
    result

-- |Test of a set of points with one degenerate dimension
test_coords_degenerate_345 :: VV.Vector (Vector Float)
test_coords_degenerate_345 = let
    sample_dists = fromLists [
        [3.0],
        [4.0, 5.0],
        [5.0, 4.0, 3.0]]
    result = coords sample_dists in
    -- assert (approx_equal_vv 1e-6 result (fromLists [
    --     [],
    --     [1.0],
    --     [0.5, 0.5 * sqrt 3]]))
    result

-- |Approx absolute comparison of floats
approx_equal :: Float -> Float -> Float -> Bool
approx_equal tol x y = abs (x - y) < tol

-- |Approx absolute comparison of vectors of floats
approx_equal_vector :: Float -> Vector Float -> Vector Float -> Bool
approx_equal_vector tol x y =
    V.length x == V.length y &&
    V.foldr (\ x y -> (approx_equal tol x 0.0) && y) True (V.zipWith (-) x y)

-- |Construct a vector of vectors from a list of lists
fromLists :: [[Float]] -> VV.Vector (Vector Float)
fromLists v = VV.fromList (map fromList v)

-- |Approx absolute comparison of vectors of vectors of floats
approx_equal_vv :: Float -> VV.Vector (Vector Float) -> VV.Vector (Vector Float) -> Bool
approx_equal_vv tol x y = 
    VV.length x == VV.length y && 
    VV.and ((VV.zipWith (approx_equal_vector tol)) x y)

main :: IO ()
main = do
    putStrLn ""
    putStrLn "Running tests..."
    putStrLn $ "test_vector_from_list: " ++ (show test_vector_from_list)
    putStrLn $ "test_sum2: " ++ (show test_sum2)
    putStrLn $ "test_sum_pair: " ++ (show test_sum_pair)
    putStrLn $ "test_sum_pair_mismatch: " ++ (show test_sum_pair_mismatch)
    putStrLn $ "test_sum_pair_mismatch2: " ++ (show test_sum_pair_mismatch2)
    putStrLn $ "test_sum2_diff: " ++ (show test_sum2_diff)
    putStrLn $ "test_sum2_diff_mismatch: " ++ (show test_sum2_diff_mismatch)
    putStrLn $ "test_coord_constructor: " ++ (show test_coord_constructor)
    putStrLn $ "test_coord_constructor_final: " ++ (show test_coord_constructor_final)
    putStrLn $ "test_coord_constructor_problem: " ++ (show test_coord_constructor_problem)
    putStrLn $ "test_calc_coord_hypercube: " ++ (show test_calc_coord_hypercube)
    putStrLn $ "test_calc_coord_simplex: " ++ (show test_calc_coord_simplex)
    putStrLn $ "test_calc_coord_first: " ++ (show test_calc_coord_first)
    putStrLn $ "test_calc_coord_second: " ++ (show test_calc_coord_second)
    putStrLn $ "test_coords_simplex: " ++ (show test_coords_simplex)
    putStrLn $ "test_coords_4d_simplex: " ++ (show test_coords_4d_simplex)
    putStrLn $ "test_coords_345: " ++ (show test_coords_345)
    putStrLn $ "test_coords_degenerate_345: " ++ (show test_coords_degenerate_345)
    putStrLn "...all tests passed"
    putStrLn ""
