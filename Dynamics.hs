module Dynamics (
    Momentum,
    State,
    Interval,
    step,
) where

import Data.Vector.Unboxed hiding (foldr1, foldr, foldM, force, sum, length, last, zipWith, map, (++))
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

-- |Dynamics defines a momentum and a position for each point
type Momentum = Double

-- |A point in state space is defined by the position and momentum of
-- |all the points.
type State = (Matrix Momentum, Matrix Coord)

-- |An interval in time
type Interval = Double

-- |Given a sets of points and momenta, and a matching set of forces,
-- |Calculate the next set of points and momenta after a small time.
step :: State -> Matrix Force -> Interval -> State
step (p, s) f dt = (p', s') where
    n = VV.length s
    -- assert (n == VV.length p && n == VV.length f)
    p' = VV.generate n gen_p
    s' = VV.generate n gen_s
    gen_p i = dynavec_p (f VV.! i) (p VV.! i) dt
    gen_s i = dynavec_s (p VV.! i) (p' VV.! i) (s VV.! i) dt

-- |Step forward the momentum of a vector of points, 
-- |passing in a vector of forces and momenta, and the incremental time
dynavec_p :: Vector Force -> Vector Momentum -> Interval -> Vector Momentum
dynavec_p f p dt = V.generate n gen where
    n = V.length p
    gen i = dynamics_p (get_or 0.0 f i) (get_or 0.0 p i) dt

-- |Step forward the distances of a vector of points, 
-- |passing in a vector of forces, momenta, new momenta, distances, and the incremental time
dynavec_s :: Vector Momentum -> Vector Momentum -> Vector Distance -> Interval -> Vector Distance
dynavec_s p p' s dt = V.generate n gen where
    n = max (V.length p) (V.length s)
    (!) = C.get_or 0.0  -- redefine ! to not walk off end of vector
    gen i = dynamics_s (p!i) (p'!i) (s!i) dt

-- |Step forward the relative momentum of one point pair, given the force. We assume that
-- |units are chosen such that we don't have to worry about mass
-- |or any other conversions
dynamics_p :: Force -> Momentum -> Interval -> Momentum
dynamics_p f p dt = p + f * dt

-- |Step forward the distance of one point pair, given the force. We assume that
-- |units are chosen such that we don't have to worry about mass
-- |or any other conversions. Note that we currently use the average
-- |momentum over the interval to calculate the change in distance.
-- |Also note that we allow distances to go negative even though this
-- |seems counter intuitive. We do not want momentum and distance to change
-- |sign relative to each other as we pass the origin. Moreover, distance
-- |is always squared in use.
dynamics_s :: Momentum -> Momentum -> Distance -> Interval -> Distance
dynamics_s p p' s dt = s + (p + p') * 0.5 * dt

-- |Test a single step for a equilateral simplex
test_step_eq_simplex :: Either String State
test_step_eq_simplex =
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
        (momenta', dists') <- test_step (momenta, dists) 1e-3
        dimensions <- C.count_dimensions dists'
        return (assert (dimensions == 3) (momenta', dists'))

-- |Test run one or more dynamics steps
test_step :: State -> Interval -> Either String State
test_step (momenta, dists) dt = do
    coords <- C.validated_coords dists
    forces <- F.forces dists coords 1.0
    return $ step (momenta, dists) forces dt

-- |Experiment with a slightly randomised equilateral simplex
test_step_simplex :: StatefulGen g m => g -> Int -> m (Either String State)
test_step_simplex rnd n = do
    simplex <- C.create_random_triangular_matrix rnd (0.99, 1.01) n
    momenta <- C.create_random_triangular_matrix rnd (0.09999, 0.10001) n
    return $ test_step (momenta, simplex) 1e-3

-- |Multiple steps with a slightly randomised equilateral simplex
test_multistep_simplex :: StatefulGen g m => g -> Int -> Int -> m (Either String Int)
test_multistep_simplex rnd m n = do
    simplex <- C.create_random_triangular_matrix rnd (0.99, 1.01) n
    momenta <- C.create_random_triangular_matrix rnd (0.09999, 0.10001) n
    return $ test_multistep (momenta, simplex) 1e-3 m n

-- |Generic multiple step tester. Runs m steps, followed by a
-- |test of the dimensionality, which should be n
test_multistep :: State -> Interval -> Int -> Int -> Either String Int
test_multistep s dt m n = do
    (momenta, dists) <- foldM (\a _ -> test_step a dt) s [0..m]
    dimensions <- C.count_dimensions dists
    return $ assert (dimensions == n) dimensions

main :: IO ()
main = do

    putStrLn ""
    putStrLn "Running Dynamics tests..."

    putStrLn $ "test_step_eq_simplex: " ++ (show test_step_eq_simplex)

    -- One generator shared by all tests that need random numbers
    rnd <- createSystemRandom
    step_simplex <- test_step_simplex rnd 5
    putStrLn $ "test_step_simplex: " ++ (show step_simplex)
    multistep_simplex <- test_multistep_simplex rnd 10 10
    putStrLn $ "test_multistep_simplex: " ++ (show multistep_simplex)
        
    putStrLn "...all tests passed"
    putStrLn ""
