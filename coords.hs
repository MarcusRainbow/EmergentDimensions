import Data.Vector.Unboxed hiding (foldr1, sum, length, last, zipWith, (++))
import qualified Data.Vector.Unboxed as V
import qualified Data.Vector as VV
import Control.Exception
import Debug.Trace

test_vector_from_list :: Vector Float
test_vector_from_list = fromList [1.0, 2.0, 3.0]

-- |Similar to Gram-Schmidt orthogonalisation but where
-- |we supply a vector of distances. This call constructs
-- |the cartesian coordinates of one point defined by
-- |its distances from the preceding points, vector s_i, and
-- |the cartesian coordinates of the preceding points p_ij
calc_coord :: Vector Float -> VV.Vector (Vector Float) -> Vector Float
calc_coord s p = constructN (V.length (VV.last p) + 1) (coord_constructor s p)

-- |Construct one of the coordinate dimensions x_n given the
-- |preceding dimensions x_i and the preceding coordinates p_ij
coord_constructor :: Vector Float -> VV.Vector (Vector Float) -> Vector Float -> Float
coord_constructor s p x = 
    let 
        b = trace ("coord_constructor s=" ++ (show s) ++ " p=" ++ (show p) ++ " x=" ++ (show x))(VV.last p)
        m = (V.length b) - 1
    in
        if (V.length x) <= m then
            (s!0^2 - s!m^2 + sum2 b - 2 * (sum_pair b x)) / (2 * b!m)
        else
            sqrt (s!m^2 - sum2_diff x b)

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

-- |Test sum of squares of diffs
test_sum2_diff :: Float
test_sum2_diff = let 
    result = sum2_diff (fromList [11.0, 22.0, 33.0]) (fromList [10.0, 20.0, 30.0]) in
    assert (result == 14.0) result

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

-- |Test coord generator for a partial hypercube (missing far corners)
test_calc_coord_hypercube :: Vector Float
test_calc_coord_hypercube = let
    root_2 = sqrt 2.0
    sample_dists = fromList [1.0, root_2, root_2, root_2]
    sample_prev = VV.fromList [empty, fromList [1.0], fromList [0.0, 1.0]]
    result = calc_coord sample_dists sample_prev in
    assert (approx_equal_vector 1e-6 result (fromList [0.0, 0.0, 1.0]))
    result

-- |Approx absolute comparison of floats
approx_equal :: Float -> Float -> Float -> Bool
approx_equal tol x y = abs (x - y) < tol

-- |Approx absolute comparison of vectors of floats
approx_equal_vector :: Float -> Vector Float -> Vector Float -> Bool
approx_equal_vector tol x y = V.foldr (\ x y -> (approx_equal tol x 0.0) && y) True (V.zipWith (-) x y)

-- |Testing some basic vector indexing
test_indexing :: Float
test_indexing =
    s!0 + s!1 + s!2^2
    where
        s = fromList [0.0, 1.0, 2.0]

main :: IO ()
main = do
    putStrLn ""
    putStrLn "Running tests..."
    putStrLn $ "test_vector_from_list: " ++ (show test_vector_from_list)
    putStrLn $ "test_sum2: " ++ (show test_sum2)
    putStrLn $ "test_sum_pair: " ++ (show test_sum_pair)
    putStrLn $ "test_sum2_diff: " ++ (show test_sum2_diff)
    putStrLn $ "test_indexing: " ++ (show test_indexing)
    putStrLn $ "test_coord_constructor: " ++ (show test_coord_constructor)
    putStrLn $ "test_coord_constructor_final: " ++ (show test_coord_constructor_final)
    putStrLn $ "test_calc_coord_hypercube: " ++ (show test_calc_coord_hypercube)
    putStrLn "...all tests passed"
    putStrLn ""
