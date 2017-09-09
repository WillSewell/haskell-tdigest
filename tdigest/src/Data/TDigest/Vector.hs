{-# LANGUAGE DataKinds             #-}
{-# LANGUAGE KindSignatures        #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE ScopedTypeVariables   #-}
module Data.TDigest.Vector where

import Control.DeepSeq        (NFData (..))
import Data.List              (sortBy)
import Data.List.NonEmpty     (NonEmpty (..), nonEmpty)
import Data.Ord               (comparing)
import Data.Proxy             (Proxy (..))
import Data.Semigroup         (Semigroup (..))
import Data.Semigroup.Reducer (Reducer (..))
import GHC.TypeLits           (KnownNat, Nat, natVal)
import Prelude ()
import Prelude.Compat

import qualified Data.Vector.Unboxed as VU

import Debug.Trace

import qualified Data.TDigest.Postprocess as PP

type Mean = Double
type Weight = Double
type Centroid = (Mean, Weight)
type Size = Int


-- | Mapping from quantile *q* to notional index *k* with compression parameter *𝛿*.
--
-- >>> ksize 42 0
-- 0.0
--
-- >>> ksize 42 1
-- 42.0
--
-- *q@ is clamped.:
--
-- >>> ksize 42 2
-- 42.0
--
ksize
    :: Double  -- ^ compression parameter, 𝛿
    -> Double  -- ^ quantile, q
    -> Double  -- ^ notional index, k
ksize comp q = comp * (asin (2 * clamp q - 1) / pi  + 0.5)

clamp :: Double -> Double
clamp x
    | x < 0.0   = 0.0
    | x > 1.0   = 1.0
    | otherwise = x

-- | Inverse of 'ksize'.
--
-- >>> ksizeInv 42 0
-- 0.0
--
-- >>> ksizeInv 42 42
-- 1.0
--
-- >>> ksizeInv 42 (ksize 42 0.3)
-- 0.3
--
ksizeInv
    :: Double  -- ^ compression parameter, 𝛿
    -> Double  -- ^ notional index, k
    -> Double  -- ^ quantile, q
ksizeInv comp k = 0.5 * (sin ((k / comp - 0.5) * pi) + 1)

merge :: Int -> Double -> [(Mean, Weight)] -> [(Mean, Weight)]
merge _   _    []     = []
merge tw' comp (y:ys) = go 0 (qLimit' 0) y ys
  where
    -- total weight
    tw = fromIntegral tw'

    qLimit' :: Double -> Double
    qLimit' q0 = ksizeInv comp (ksize comp q0 + 1)  -- k⁻¹ (k (q₀, 𝛿) + 1, 𝛿)

    go :: Double         -- q0
       -> Double         -- qLimit
       -> (Mean, Weight)   -- sigma
       -> [(Mean, Weight)]
       -> [(Mean, Weight)]
    go _q0 _qLimit sigma [] = [sigma] -- C'.append(σ)
    go  q0  qLimit sigma (x:xs)
        | q <= qLimit = go q0 qLimit (plus sigma x) xs
        | otherwise   = sigma : go q0' (qLimit' q0') x xs
-- traceShow ("q", sigma, x, q, qLimit) $
      where
        q = q0 + (snd sigma + snd x) / tw
        q0' = q0 + snd sigma / tw

    plus :: Centroid -> Centroid -> Centroid
    plus (m1,w1) (m2,w2) = ((m1 * w1 + m2 * w2) / w, w) where w = w1 + w2

data TDigest (compression :: Nat) = TDigest
    { tdigestTotalWeight :: !Int                   -- ^ sum of vector and buffer size
    , tdigestData        :: !(VU.Vector Centroid)  -- ^ actual data, *invariant:* sorted by mean
    , tdigestBufferSize  :: !Int
    , tdigestBuffer      :: [Double]               -- ^ addition buffer, elements with weight 1
    , tdigestDirection   :: !Bool                  -- ^ direction is a hack, so we merge from left and right
    }
  deriving Show

forceCompress :: forall comp. KnownNat comp => TDigest comp -> TDigest comp
forceCompress t@(TDigest s d _ b dir)
    | null b    = t
    | otherwise = TDigest s d' 0 [] (not dir)
  where
    d' = VU.fromList
       . rev
       . merge s comp             -- compress
       . rev
       . sortBy (comparing fst)   -- sort
       . (++ map (flip (,) 1) b)  -- add buffer
       . VU.toList
       $ d
    comp = fromInteger $ natVal (Proxy :: Proxy comp) * 50
    rev | dir       = id
        | otherwise = reverse

instance KnownNat comp => Semigroup (TDigest comp) where

instance KnownNat comp => Monoid (TDigest comp) where
    mempty = (TDigest 0 mempty 0 mempty True)
    mappend = (<>)

-- | Both 'cons' and 'snoc' are 'insert'
instance KnownNat comp => Reducer Double (TDigest comp) where
    cons = insert
    snoc = flip insert
    unit = singleton

-- | Insert single value into 'TDigest'.
insert
    :: KnownNat comp
    => Double  -- ^ element
    -> TDigest comp
    -> TDigest comp
insert x (TDigest s d sb b dir) = compress (TDigest (s + 1) d (sb + 1) (x : b) dir)

singleton :: Double -> TDigest comp
singleton x = TDigest 1 (VU.singleton (x, 1)) 0 [] True

compress :: forall comp. KnownNat comp => TDigest comp -> TDigest comp
compress t@(TDigest _ _ bs _ _)
    | bs > compInt = forceCompress t
    | otherwise    = t
  where
    compInt = fromInteger $ natVal (Proxy :: Proxy comp) * 50

median :: KnownNat comp => TDigest comp -> Maybe Double
median = quantile 0.5

quantile :: KnownNat comp => Double -> TDigest comp -> Maybe Double
quantile q td = case forceCompress td of
    (TDigest tw d _ _ _) -> PP.quantile' q (fromIntegral tw) <$> histogram d

histogram :: VU.Vector Centroid -> Maybe (NonEmpty PP.HistBin)
histogram = fmap PP.histogram' . nonEmpty . VU.toList

validate :: TDigest comp -> Either String (TDigest comp)
validate td@(TDigest tw d bs b dir)
    | otherwise = Right td
