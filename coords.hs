import Data.Vector.Unboxed hiding (foldr1, sum, length, last, zipWith, map, (++))
import qualified Data.Vector.Unboxed as V
import qualified Data.Vector as VV
import Data.Maybe
import Data.Either
import Numeric.IEEE
import System.Random
import System.Random.Stateful
import Control.Exception
import Debug.Trace
import System.Environment
import System.Random.MWC as MWC
import System.Random.MWC.Distributions as MWC

-- |Tolerance used for the absolute comparison of coordinate generation etc.
tolerance :: Float
tolerance = 1e-6

-- |Given a set of distances between points, count
-- |the dimensions. If there is an error, such as a set of
-- |distances that do not fit into any possible space, a
-- |explanatory string is returned. Otherwise the count. 
count_dimensions :: VV.Vector (Vector Float) -> Either String Int
count_dimensions s = do
    c <- coords s
    t <- uncoords c
    if approx_equal_vv tolerance s t then
        return (if VV.null c then 0 else V.length (VV.last c))
    else
        Left "Not all distances round-trip successfully in Euclidean space"

-- |Given a set of distances between points, calculate
-- |a consistent set of coordinates for each point. The
-- |distances are supplied as a vector of vectors as
-- |an encoding of a symmetric matrix with diagonal zero.
-- |The coordinates are returned as vectors where zero
-- |coordinates may be zero-suppressed on the right hand
-- |side. Thus the origin is empty vector. This allows
-- |the implementation to be agnostic of the number of
-- |dimensions. If the coordinates cannot be represented,
-- |an error string is returned.
coords :: VV.Vector (Vector Float) -> Either String (VV.Vector (Vector Float))
coords s = validate_non_nan $ VV.constructN (VV.length s + 1) (coords_constructor s)

-- |Construct one of the coordinate sets x_j given the distances
-- |s_ij and the previous coordinate sets p_ij
coords_constructor ::  VV.Vector (Vector Float) ->  VV.Vector (Vector Float) -> Vector Float
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
    otherwise -> let n = V.length (VV.last p) + 1 in 
        trim1 tolerance (constructN n (coord_constructor s p))

-- |Validates a vector of vectors, checking the values for NaN
validate_non_nan :: VV.Vector (Vector Float) -> Either String (VV.Vector (Vector Float))
validate_non_nan v = 
    if VV.any (V.any isNaN) v then
        Left $ "At least one coordinate is NaN in " ++ show v
    else
        Right v

-- |Construct one of the coordinate dimensions x_n given the
-- |preceding dimensions x_i and the preceding coordinates p_ij
coord_constructor :: Vector Float -> VV.Vector (Vector Float) -> Vector Float -> Float
coord_constructor s p x = 
    let 
        -- n = trace ("coord_constructor: s=" ++ (show s) ++ " p=" ++ (show p) ++ " x=" ++ (show x)) (V.length x) + 1
        n = 1 + V.length x
        m = if (VV.null p) then 0 else V.length (VV.last p)
    in
        -- we do not validate against points we are skipping, but that's okay
        -- because we can finally validate by using uncoord
        if n <= m then let b = first_coord_of_dim n p in
            (s!0^2 - s!n^2 + sum2 b - 2 * (sum_pair b x)) / (2 * V.last b)
        else
            case safe_sqrt ((V.last s)^2 - sum2_diff x (VV.last p)) of
                Just x  -> x
                Nothing -> nan -- sadly there's no way to pass monads through constructN

-- |Scans through the vector of previous points to find the first with
-- |the requisite dimension. Must exist (throws error if not)
first_coord_of_dim :: Int -> VV.Vector (Vector Float) -> Vector Float
first_coord_of_dim n p = fromJust (VV.find (\ v -> V.length v == n) p)

-- |Trim zeroes off the right end of a vector. This call trims no
-- |more than one zero. (Useful because we want the vector of
-- |coordinates to be monotonic increasing in length.)
trim1 :: Float -> Vector Float -> Vector Float
trim1 tol v = let l = V.length v in
    if l > 0 && approx_equal tol (V.last v) 0.0 then
        unsafeTake (l - 1) v
    else
        v

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

-- |Reverse the action of coords. Given a set of coordinates, calculate
-- |the distances between them. If there is an error, returns an error
-- |string.
uncoords :: VV.Vector (Vector Float) -> Either String (VV.Vector (Vector Float))
uncoords c =
    VV.generateM (VV.length c - 1) (uncoords_generator c)

-- |Construct a set of distances from one of a set of coordinates
uncoords_generator :: VV.Vector (Vector Float) -> Int -> Either String (Vector Float)
uncoords_generator c i =
    generateM (i + 1) (find_distance c (i + 1))

-- |Find the distance between two points given the set of coordinates
find_distance :: VV.Vector (Vector Float) -> Int -> Int -> Either String Float
find_distance c i j = case calc_distance (c VV.! i) (c VV.! j) of
    Just x  -> Right x
    Nothing -> Left $ "sqrt negative distance between points " ++ show i ++ " and " ++ show j 

-- |Calculate the distance between two coordinates
calc_distance :: Vector Float -> Vector Float -> Maybe Float
calc_distance x y = safe_sqrt ((sum2_diff x y) + (sum2_rem x y))

-- |Given two vectors of different lengths, sum the squares of the
-- |tail elements in the longer vector
sum2_rem :: Vector Float -> Vector Float -> Float
sum2_rem x y =
    let
        lx = V.length x
        ly = V.length y
        d = lx - ly
    in case compare d 0 of
        LT -> sum2 (unsafeDrop lx y)
        EQ -> 0.0
        GT -> sum2 (unsafeDrop ly x)

-- |Safer version of sqrt, which will not give NaN for slightly negative
safe_sqrt :: Float -> Maybe Float
safe_sqrt x = 
    if x >= 0.0 then
        Just (sqrt x)
    else 
        if x > -tolerance then
            Just 0.0
        else
            Nothing

-- |Approx absolute comparison of floats
approx_equal :: Float -> Float -> Float -> Bool
approx_equal tol x y = abs (x - y) < tol

-- |Approx absolute comparison of vectors of floats
approx_equal_vector :: Float -> Vector Float -> Vector Float -> Bool
approx_equal_vector tol x y =
    V.length x == V.length y &&
    V.foldr (\ x y -> (approx_equal tol x 0.0) && y) True (V.zipWith (-) x y)

-- |Approx absolute comparison of vectors of vectors of floats
approx_equal_vv :: Float -> VV.Vector (Vector Float) -> VV.Vector (Vector Float) -> Bool
approx_equal_vv tol x y = 
    VV.length x == VV.length y && 
    VV.and ((VV.zipWith (approx_equal_vector tol)) x y)

match_uncoords :: Float -> VV.Vector (Vector Float) -> VV.Vector (Vector Float) -> Bool
match_uncoords tol original result = case uncoords result of
    Left _ -> False
    Right roundtrip -> approx_equal_vv tol original roundtrip

-- |Construct a vector of vectors from a list of lists
fromLists :: [[Float]] -> VV.Vector (Vector Float)
fromLists v = VV.fromList (map fromList v)

-- |Construct a vector of vectors as a triangle, given Monadic random generator,
-- |a range, and a size.
create_random_triangular_matrix :: StatefulGen g m => g -> (Float, Float) -> Int -> m (VV.Vector (Vector Float))
create_random_triangular_matrix rnd range n = 
    VV.generateM n (\ i -> create_random_vector rnd range (i + 1))

-- |Construct a vector of random numbers uniform in the range 0..1, given a Monadic
-- |random generator and a size.
create_random_vector :: StatefulGen g m => g -> (Float, Float) -> Int -> m (Vector Float)
create_random_vector rnd range n =
    V.replicateM n (uniformRM range rnd)

-- |Various tests of safe_sqrt
test_safe_sqrt :: Maybe Float
test_safe_sqrt = do
    test1 <- safe_sqrt 4.0
    test2 <- safe_sqrt (-1e-7)
    let result = test1 + test2
    return (assert (result == 2.0) result)

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

-- |Test sum2_rem where second vector is empty
test_sum2_rem :: Float
test_sum2_rem = let 
    result = sum2_rem (fromList [3.0]) empty in
    assert (result == 9.0) 
    result

-- |Test sum2_rem where first vector is longer
test_sum2_rem_gt :: Float
test_sum2_rem_gt = let 
    result = sum2_rem (fromList [1.0, 2.0, 3.0, 4.0]) (fromList [10.0, 20.0]) in
    assert (result == 25.0) result

-- |Test sum2_rem where second vector is longer
test_sum2_rem_lt :: Float
test_sum2_rem_lt = let 
    result = sum2_rem (fromList [1.0, 2.0]) (fromList [10.0, 20.0, 30.0, 40.0]) in
    assert (result == 2500.0) result

-- Test calc_distance for two vectors of different lengths
test_calc_distance :: Maybe Float
test_calc_distance = do
    result <- calc_distance (fromList [3.0]) (fromList [])
    return (assert (result == 3.0) result) 

-- |Test coord constructor for intermediate step
test_coord_constructor :: Float
test_coord_constructor = let
    root_2 = sqrt 2.0
    sample_dists = fromList [1.0, root_2, root_2, root_2]
    sample_prev = fromLists [[], [1.0], [0.0, 1.0]]
    sample_coords = fromList [0.0]
    result = coord_constructor sample_dists sample_prev sample_coords in
    assert (approx_equal 1e-6 result 0.0) result
    
-- |Test coord constructor for final step
test_coord_constructor_final :: Float
test_coord_constructor_final = let
    root_2 = sqrt 2.0
    sample_dists = fromList [1.0, root_2, root_2, root_2]
    sample_prev = fromLists [[], [1.0], [0.0, 1.0]]
    sample_coords = fromList [0.0, 0.0]
    result = coord_constructor sample_dists sample_prev sample_coords in
    assert (approx_equal 1e-6 result 1.0) result
    
-- |Test coord constructor in case that gave out of bounds error
test_coord_constructor_problem :: Float
test_coord_constructor_problem = let
    sample_dists = fromList [1.0, 1.0, 1.0]
    sample_prev = fromLists [[], [1.0]]
    sample_coords = fromList [0.5]
    result = coord_constructor sample_dists sample_prev sample_coords in
    assert (approx_equal 1e-6 result 0.8660254) 
    result
    
-- |Test coord constructor in case that gave NaN
test_coord_constructor_degenerate1 :: Float
test_coord_constructor_degenerate1 = let
    sample_dists = fromList [5.0, 4.0, 3.0]
    sample_prev = fromLists [[], [3.0], [0.0, 4.0]]
    sample_coords = fromList []
    result = coord_constructor sample_dists sample_prev sample_coords in
    assert (approx_equal 1e-6 result 3.0) 
    result

-- |Test coord constructor in case that gave NaN
test_coord_constructor_degenerate2 :: Float
test_coord_constructor_degenerate2 = let
    sample_dists = fromList [5.0, 4.0, 3.0]
    sample_prev = fromLists [[], [3.0], [0.0, 4.0]]
    sample_coords = fromList [3.0]
    result = coord_constructor sample_dists sample_prev sample_coords in
    assert (approx_equal 1e-6 result 4.0) 
    result

-- |Test coord constructor in case that gave invalid final coordinate
test_coord_constructor_degenerate_plus :: Float
test_coord_constructor_degenerate_plus = let
    sample_dists = fromList [12.0, 12.369317, 12.649111, 13.0]
    sample_prev = fromLists [[], [3.0], [0.0, 4.0], [3.0, 4.0]]
    sample_coords = fromList [0.0, 0.0]
    result = coord_constructor sample_dists sample_prev sample_coords in
    assert (approx_equal 1e-6 result 12.0) 
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

-- |Test coord generator for a simplex
test_calc_coord_simplex :: Vector Float
test_calc_coord_simplex = let
    sample_dists = fromList [1.0, 1.0]
    sample_prev = fromLists [[], [1.0]]
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
test_coords_simplex :: Either String (VV.Vector (Vector Float))
test_coords_simplex = let
    sample_dists = fromLists [
        [1.0],
        [1.0, 1.0],
        [1.0, 1.0, 1.0]]
    expected = fromLists [
        [],
        [1.0],
        [0.5, 0.8660254],
        [0.5,0.28867513,0.8164966]]
    in test_coords sample_dists expected

-- |Helper to run any overall test, checking coords and uncoords
test_coords :: VV.Vector (Vector Float) -> VV.Vector (Vector Float) -> Either String (VV.Vector (Vector Float))
test_coords sample_dists expected =
    case coords sample_dists of
        Left result ->
            error result
            Left result 
        Right result -> 
            assert (approx_equal_vv 1e-6 result expected)
            assert (match_uncoords 1e-6 sample_dists result)
            Right result

-- |Overall test, given a collection of 5 points all one apart
test_coords_4d_simplex :: Either String (VV.Vector (Vector Float))
test_coords_4d_simplex = let
    sample_dists = fromLists [
        [1.0],
        [1.0, 1.0],
        [1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0]]
    expected = fromLists [
        [],
        [1.0],
        [0.5, 0.8660254],
        [0.5,0.28867513,0.8164966],
        [0.5,0.28867513,0.20412417,0.7905694]]
    in test_coords sample_dists expected

-- |Round trip test of a near-equilateral 4d simplex
test_roundtrip_4d_simplex :: Either String (VV.Vector (Vector Float))
test_roundtrip_4d_simplex = let
    sample_dists = fromLists [
        [1.01],
        [1.02, 1.03],
        [1.04, 1.05, 1.06],
        [1.07, 1.08, 1.09, 1.0]]
    expected = fromLists [
        [],
        [1.01],
        [0.49485147,0.8919204],
        [0.49465352,0.28524965,0.8692241],
        [0.49435645,0.2847417,0.33074602,0.8426393]]
    in test_coords sample_dists expected

-- |Test of a set of points that form a 3, 4, 5 triangle
test_coords_345 :: Either String (VV.Vector (Vector Float))
test_coords_345 = let
    sample_dists = fromLists [
        [3.0],
        [4.0, 5.0]]
    expected = fromLists [
        [],
        [3.0],
        [0.0, 4.0]]
    in test_coords sample_dists expected

-- |Test of a set of points that form a 3, 4, 5 triangle plus a third dimension
test_coords_345_plus :: Either String (VV.Vector (Vector Float))
test_coords_345_plus = let
    sample_dists = fromLists [
        [3.0],
        [4.0, 5.0],
        [12.0, (sqrt 153.0), (sqrt 160.0)]]
    expected = fromLists [
        [],
        [3.0],
        [0.0, 4.0],
        [0.0, 0.0, 12.0]]
    in test_coords sample_dists expected

-- |Test of a set of points with one degenerate dimension
test_coords_degenerate_345 :: Either String (VV.Vector (Vector Float))
test_coords_degenerate_345 = let
    sample_dists = fromLists [
        [3.0],
        [4.0, 5.0],
        [5.0, 4.0, 3.0]]
    expected = fromLists [
        [],
        [3.0],
        [0.0, 4.0],
        [3.0, 4.0]]
    in test_coords sample_dists expected

-- |Test of a set of points with one degenerate dimension plus another point
test_coords_degenerate_345_plus :: Either String (VV.Vector (Vector Float))
test_coords_degenerate_345_plus = let
    sample_dists = fromLists [
        [3.0],
        [4.0, 5.0],
        [5.0, 4.0, 3.0],
        [12.0, (sqrt 153.0), (sqrt 160.0), 13.0]]
    expected = fromLists [
        [],
        [3.0],
        [0.0, 4.0],
        [3.0, 4.0],
        [0.0, 0.0, 12.0]]
    in test_coords sample_dists expected

-- |Test of a set of points with one degenerate dimension plus another point
test_count_dimensions_345_plus :: Either String Int
test_count_dimensions_345_plus = let
    sample_dists = fromLists [
        [3.0],
        [4.0, 5.0],
        [5.0, 4.0, 3.0],
        [12.0, (sqrt 153.0), (sqrt 160.0), 13.0]]
    in test_count_dimensions sample_dists 3

-- |Run any test of count_dimensions
test_count_dimensions :: VV.Vector (Vector Float) -> Int -> Either String Int
test_count_dimensions sample n = case count_dimensions sample of
    Left e ->
        error e
        Left e
    Right result ->
        assert (result == n)
        Right result

test_count_dimensions_simplex :: StatefulGen g m => g -> Int -> m (Either String Int)
test_count_dimensions_simplex rnd n = do
    simplex <- create_random_triangular_matrix rnd (0.99, 1.01) n
    return (test_count_dimensions simplex n)

-- |Test the reverse calculation of distances from coords
test_uncoords :: Either String (VV.Vector (Vector Float))
test_uncoords = let
    sample_coords = fromLists [
        [],
        [3.0],
        [0.0, 4.0],
        [3.0, 4.0, 0.0]]
    in case uncoords sample_coords of
        Left e ->
            error e
            Left e
        Right result ->
            assert (approx_equal_vv 1e-6 result (fromLists [
                [3.0],
                [4.0, 5.0],
                [5.0, 4.0, 3.0]]))
            Right result

-- |Test list list constructor
test_from_lists :: VV.Vector (Vector Float)
test_from_lists = fromLists [[], [1.0], [1.0, 2.0]]

-- |Ultra-simplistic argument parsing. Either returns 0
-- |or a number of points to test
parse_args :: [String] -> Int
parse_args [] = 0
parse_args (arg:[]) = read arg
parse_args _ = -1

main :: IO ()
main = do

    args <- getArgs
    let arg = parse_args args
    let n = if arg > 0 then arg else 30

    putStrLn ""
    putStrLn "Running tests..."

    putStrLn $ "test_from_lists: " ++ (show test_from_lists)
    putStrLn $ "test_safe_sqrt: " ++ (show test_safe_sqrt)
    putStrLn $ "test_sum2: " ++ (show test_sum2)
    putStrLn $ "test_sum_pair: " ++ (show test_sum_pair)
    putStrLn $ "test_sum_pair_mismatch: " ++ (show test_sum_pair_mismatch)
    putStrLn $ "test_sum_pair_mismatch2: " ++ (show test_sum_pair_mismatch2)
    putStrLn $ "test_sum2_diff: " ++ (show test_sum2_diff)
    putStrLn $ "test_sum2_diff_mismatch: " ++ (show test_sum2_diff_mismatch)
    putStrLn $ "test_sum2_rem: " ++ (show test_sum2_rem)
    putStrLn $ "test_sum2_rem_lt: " ++ (show test_sum2_rem_lt)
    putStrLn $ "test_sum2_rem_gt: " ++ (show test_sum2_rem_gt)
    putStrLn $ "test_calc_distance: " ++ (show test_calc_distance)
    putStrLn $ "test_coord_constructor: " ++ (show test_coord_constructor)
    putStrLn $ "test_coord_constructor_final: " ++ (show test_coord_constructor_final)
    putStrLn $ "test_coord_constructor_problem: " ++ (show test_coord_constructor_problem)
    putStrLn $ "test_coord_constructor_degenerate1: " ++ (show test_coord_constructor_degenerate1)
    putStrLn $ "test_coord_constructor_degenerate2: " ++ (show test_coord_constructor_degenerate2)
    putStrLn $ "test_coord_constructor_degenerate_plus: " ++ (show test_coord_constructor_degenerate_plus)
    putStrLn $ "test_calc_coord_hypercube: " ++ (show test_calc_coord_hypercube)
    putStrLn $ "test_calc_coord_simplex: " ++ (show test_calc_coord_simplex)
    putStrLn $ "test_calc_coord_first: " ++ (show test_calc_coord_first)
    putStrLn $ "test_calc_coord_second: " ++ (show test_calc_coord_second)
    putStrLn $ "test_coords_simplex: " ++ (show test_coords_simplex)
    putStrLn $ "test_coords_4d_simplex: " ++ (show test_coords_4d_simplex)
    putStrLn $ "test_roundtrip_4d_simplex: " ++ (show test_roundtrip_4d_simplex)
    putStrLn $ "test_coords_345: " ++ (show test_coords_345)
    putStrLn $ "test_coords_345_plus: " ++ (show test_coords_345_plus)
    putStrLn $ "test_coords_degenerate_345: " ++ (show test_coords_degenerate_345)
    putStrLn $ "test_coords_degenerate_345_plus: " ++ (show test_coords_degenerate_345_plus)
    putStrLn $ "test_uncoords: " ++ (show test_uncoords)
    putStrLn $ "test_count_dimensions_345_plus: " ++ (show test_count_dimensions_345_plus)

    -- One generator shared by all tests that need random numbers
    rnd <- createSystemRandom
    count_dimensions_simplex <- test_count_dimensions_simplex rnd n
    putStrLn $ "test_count_large_random_simplex: " ++ (show count_dimensions_simplex)

    putStrLn "...all tests passed"
    putStrLn ""
