"""Kramers--Kronig counterparts and moments of thresholded Drude bands.

The subtraction is the KK transform of the removed interval [0, gap], not
a fit or an independently adjusted real response. See CONSISTENCY.md.
"""
from functools import lru_cache

import numpy as np

RESPONSE_VERSION = "threshold-kk-moment-budget-v2"


@lru_cache(maxsize=4)
def _quadrature(order):
    x, w = np.polynomial.legendre.leggauss(order)
    return (x + 1.) / 2., w / 2.


def imaginary(energy, centre, width, derivative=False):
    """Unit-amplitude imaginary Drude susceptibility, without Ep squared."""
    a, b = centre**2 - energy**2, width * energy
    den = a*a + b*b
    return 2*b**3 / den**2 if derivative else b / den


def retained_moment(centre, width, gap, derivative=False, order=64):
    """Fraction of the unit Drude energy-weighted moment above the gap."""
    centre, width = np.broadcast_arrays(centre, width)
    if (not np.isfinite(gap) or gap <= 0 or np.any(~np.isfinite(centre))
            or np.any(~np.isfinite(width)) or np.any(centre <= 0)
            or np.any(width <= 0)):
        raise ValueError("Positive finite Drude scales and gap are required")
    x, w = _quadrature(order)
    energy = float(gap) * x
    removed = (2 / np.pi * float(gap)
        * np.sum(w * energy * imaginary(energy, centre[..., None],
                                       width[..., None], derivative), axis=-1))
    return 1. - removed


def _removed_integral(energy, centre, width, gap, derivative):
    """Complex rational integral; its imaginary part is the missing KK term."""
    # Susceptibility = 1/D; derivative-Drude = (1 - gamma*d/dgamma)/D,
    # D = E0^2 - z^2 - i*gamma*z. All roots lie below the real axis.
    root = np.sqrt(4*centre**2 - width**2 + 0j)
    r2 = (-1j*width - root)/2
    r1 = -centre**2/r2  # product identity avoids cancellation when overdamped
    zero = energy == 0
    e = np.where(zero, 1., energy)
    lp = np.log(np.abs((gap - e) / e))
    lm = np.log1p(gap / e)

    def kernel(r):
        lr = np.log(1 - gap/r)
        a, b = e-r, -e-r
        value = .5*((lp-lr)/a + (lm-lr)/b)
        value = np.where(zero, lr/r, value)
        if not derivative:
            return value, None
        dlr = gap / (r*(r-gap))
        deriv = .5*(-dlr/a + (lp-lr)/a**2 - dlr/b + (lm-lr)/b**2)
        deriv = np.where(zero, dlr/r-lr/r**2, deriv)
        return value, deriv

    k1, d1 = kernel(r1)
    k2, d2 = kernel(r2)
    integral = -(k1-k2)/root
    if derivative:
        dr = -width/root
        dr2 = (-1j - dr)/2
        dr1 = -r1*dr2/r2
        dintegral = dr*(k1-k2)/root**2 - (d1*dr1-d2*dr2)/root
        integral -= width*dintegral
    return 2/np.pi * integral.imag


def removed_real(energy, centre, width, gap, derivative=False, order=64):
    """KK real contribution of the interval removed by the common onset.

    Arrays broadcast normally. The hard edge has a logarithmic real-part
    singularity; its value is represented by the nearest upper floating-point
    energy. The imaginary response is zero exactly at that edge. Values at
    zero and nearly repeated Drude poles use regularized Gauss quadrature.
    """
    e, c, g = np.broadcast_arrays(np.abs(np.asarray(energy, float)),
                                 np.asarray(centre, float),
                                 np.asarray(width, float))
    if (not np.isfinite(gap) or gap <= 0 or np.any(~np.isfinite(e))
            or np.any(~np.isfinite(c)) or np.any(~np.isfinite(g))
            or np.any(c <= 0) or np.any(g <= 0)):
        raise ValueError("Finite energies and positive Drude scales are required")
    shape = e.shape
    e, c, g = (a.ravel() for a in (e, c, g))
    e = np.where(e == gap, np.nextafter(float(gap), np.inf), e)
    result = np.empty_like(e)
    # The analytic divided differences lose precision at repeated poles and
    # at tiny E. This threshold selects a numerical representation only.
    fallback = (((e > 0) & (e < 1e-4*gap))
                | (np.abs(4*c*c-g*g) < 1e-3*(4*c*c+g*g)))
    mask = ~fallback
    if np.any(mask):
        result[mask] = _removed_integral(e[mask], c[mask], g[mask], gap, derivative)
    x, w = _quadrature(order)
    z = gap*x
    fallback_indices = np.flatnonzero(fallback)
    for start in range(0, fallback_indices.size, 2048):
        ix = fallback_indices[start:start+2048]
        ee, cc, gg = e[ix, None], c[ix, None], g[ix, None]
        h = z*imaginary(z, cc, gg, derivative)
        he = ee*imaginary(ee, cc, gg, derivative)
        den = (z-ee)*(z+ee)
        # Nodes strictly inside (0,gap); no endpoint evaluations are needed.
        with np.errstate(divide="ignore", invalid="ignore"):
            remainder = (h-he)/den
            boundary = np.log(np.abs((gap-ee)/(gap+ee)))/(2*ee)
        a, b = cc*cc-ee*ee, gg*ee
        d, dd = a*a+b*b, -4*ee*a+2*gg*b
        if derivative:
            dh = 8*gg**3*ee**3/d**2 - 4*gg**3*ee**4*dd/d**3
        else:
            dh = 2*gg*ee/d - gg*ee**2*dd/d**2
        limit = np.divide(dh, 2*ee, out=np.zeros_like(dh), where=ee != 0)
        remainder = np.where(den == 0, limit, remainder)
        boundary = np.where(ee == 0, -1/gap, boundary)
        result[ix] = 2/np.pi*(gap*np.sum(w*remainder, axis=1)
                              + (he*boundary)[:, 0])
    if np.any(~np.isfinite(result)):
        raise FloatingPointError("Thresholded Drude KK evaluation is nonfinite")
    return result.reshape(shape)
