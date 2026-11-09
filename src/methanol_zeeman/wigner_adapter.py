from functools import lru_cache
from sympy import Rational
from sympy.physics.wigner import wigner_3j
def initialize(): pass
@lru_cache(None)
def wig3jj(*x): return float(wigner_3j(*[Rational(int(i),2) for i in x]))
