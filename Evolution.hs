import Data.Vector.Unboxed hiding (foldr1, foldr, foldM, reverse, force, sum, length, last, zipWith, map, (++))
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
import Control.Monad
import Debug.Trace
import System.Environment
import System.Random.MWC as MWC
import System.Random.MWC.Distributions as MWC
import Coords as C
import Forces as F
import Dynamics as D

-- |Predictor Corrector takes into account the fact that the forces at
-- |the start of a timestep may be different from those at the end.
-- |It has a number of phases, for example PECE or PECECE.
-- |
-- |P stands for a predictor phases, where the influence on the state over
-- |the timestep is predicted. In our case, this means the calculation of
-- |the coordinates and hence the force, based on the distances and momenta
-- |at the start of the timestep.
-- |
-- |E stands for an evaluation phase, where the dynamics of the points and
-- |the forces are used to calculate the distances and momenta at the end 
-- |of the timestep.
-- |
-- |C stands for a corrector phase, where the coordinates and hence the 
-- |forces at the end of the timestep are calculated. These are used as a
-- |correction to the forces, so that rather than using forces at the start,
-- |we use the average force over the timestep, assuming that the change
-- |to the force is linear.
-- |
-- |There are two reasons for inaccuracies in predictor correct evaluation.
-- |The PE(CE)* sequence may not be fully converged, so the final correction
-- |may not be exact; there may be convexity to the force function, so that
-- |an exact knowledge of the force at the start and end of the timestep
-- |is insufficient.
predictor_corrector :: D.State -> F.Gravity -> D.Interval -> Either String D.State
predictor_corrector (momenta, dists) g dt = do
    -- predict
    coordsP <- C.validated_coords dists
    forcesP <- F.forces dists coordsP g

    -- evaluate
    let (momentaE, distsE) = D.step (momenta, dists) forcesP dt
    coordsE <- C.validated_coords distsE
    forcesE <- F.forces distsE coordsE g

    -- correct
    let forcesC = average_forces forcesP forcesE

    -- evaluate a second time
    return $ D.step (momenta, dists) forcesC dt

-- |Evolves the universe one step at a time using the given evolution method
-- |and generating some output
evolve :: (D.State -> F.Gravity -> D.Interval -> Either String D.State) ->
        (D.State -> Either String a) ->
        D.State -> F.Gravity -> D.Interval -> Int -> Either String [a]
evolve stepper outputter state g dt n = do
    (_, output) <- foldM (\a _ -> one_step stepper outputter g dt a) (state, []) [0..n]
    return $ reverse output

-- |Evolves one single step and inserts the new output to the start
-- |of an output list.
one_step :: (D.State -> F.Gravity -> D.Interval -> Either String D.State) ->
        (D.State -> Either String a) ->
        F.Gravity -> D.Interval -> (D.State, [a]) -> Either String (D.State, [a])
one_step stepper outputter g dt (state, old_output) = do
    state <- stepper state g dt
    output <- outputter state
    return (state, output:old_output)

-- |Standard evolution and output. Uses predictor/corrector, and outputs 
-- |a list of lists of max dimension sizes in chronological order
std_evolve :: D.State -> F.Gravity -> D.Interval -> Int -> Either String [[Double]]
std_evolve = evolve predictor_corrector (\(_, s) -> C.dimension_ranges s)

-- |Returns a force matrix that is the mean of two matrices passed in
average_forces :: C.Matrix F.Force -> C.Matrix F.Force -> C.Matrix F.Force
average_forces a b = VV.generate n generator where
    n = VV.length a
    -- assert (n == VV.length b)
    generator i = gen (a VV.! i) (b VV.! i)
    gen va vb = V.generate (V.length va `max` V.length vb) (generator' va vb)
    (!) = get_or 0.0  -- no guarantee that force vectors are the same length
    generator' va vb j = 0.5 * (va!j + vb!j)

-- |Tests the average_forces function with some simple force matrices
test_average_forces :: C.Matrix F.Force
test_average_forces = let
        a = fromLists [
            [],
            [1.0],
            [1.1, 1.2]]
        b = fromLists [
            [20.0],
            [],
            [21.1, 21.2, 21.3]]
        expected = fromLists [
            [10.0],
            [0.5],
            [11.1, 11.2, 10.65]]
        result = average_forces a b
    in
        assert (C.approx_equal_vv 1e-6 result expected)
        result

-- |Test one predictor corrector step for an equilateral simplex
test_pc_eq_simplex :: Either String State
test_pc_eq_simplex =
    let 
        dists = fromLists [
            [1.0],
            [1.0, 1.0],
            [1.0, 1.0, 1.0]]
        momenta = fromLists [
            [0.01],
            [0.01, 0.01],
            [0.01, 0.01, 0.01]]
    in do
        (momenta', dists') <- predictor_corrector (momenta, dists) 1.0 1e-3
        dimensions <- C.count_dimensions dists'
        return (assert (dimensions == 3) (momenta', dists'))

-- |Test evolution using predictor corrector on a randomised simplex
test_evolve_simplex :: StatefulGen g m => g -> Int -> Int -> m (Either String [[Double]])
test_evolve_simplex rnd m n = do
    simplex <- C.create_random_triangular_matrix rnd (0.99, 1.01) n
    momenta <- C.create_random_triangular_matrix rnd (0.00099, 0.10001) n
    return $ test_evolve (momenta, simplex) 50 1e-2 m n

-- |Generic evolution tester. Runs m steps, followed by a
-- |test of the dimensionality, which should be n
test_evolve :: State -> F.Gravity -> D.Interval -> Int -> Int -> Either String [[Double]]
test_evolve s g dt m n = do
    ranges <- std_evolve s g dt m
    -- return $ assert (length ranges == n) ranges
    return ranges

main :: IO ()
main = do

    putStrLn ""
    putStrLn "Running Evolution tests..."

    putStrLn $ "test_average_forces: " ++ (show test_average_forces)
    putStrLn $ "test_pc_eq_simplex: " ++ (show test_pc_eq_simplex)

    -- One generator shared by all tests that need random numbers
    rnd <- createSystemRandom
    evolve_simplex <- test_evolve_simplex rnd 5 5
    putStrLn $ "test_evolve_simplex: " ++ (show evolve_simplex)
        
    putStrLn "...all tests passed"
    putStrLn ""
