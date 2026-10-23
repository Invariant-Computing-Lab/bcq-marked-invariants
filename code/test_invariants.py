import numpy as np
from fractions import Fraction
from bcq_invariants import symplectic_J, pairing, random_symplectic, hrr_pairing_line_bundles


def test_symplectic_pairing_exact():
    rng=np.random.default_rng(123)
    J=symplectic_J(3)
    x=np.array([1,2,-1,3,0,4],dtype=int); y=np.array([0,1,3,-2,2,1],dtype=int)
    p=pairing(x,y,J)
    for _ in range(100):
        U=random_symplectic(3,rng,stages=2)
        assert pairing(U@x,U@y,J)==p


def test_hrr_common_twist_exact():
    for a in range(-5,6):
        for b in range(-5,6):
            for t in [-3,-1,1,4]:
                assert hrr_pairing_line_bundles(a,b)==hrr_pairing_line_bundles(a+t,b+t)
