"""Independent KK integrals, true moments, edge limits and query-grid checks."""
import numpy as np
import pytest
from scipy.integrate import quad

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python_scripts" / "physics_ice"))
import causal_drude as causal
import emfietzoglou_model_finite_q as model
import emfietzoglou_model_optical_limit as optical

IS_ION = False


def coefficients():
    return model.DispersionCoeffs([3.82, 2.47, 2.47, 3.01, 2.44],
        [.0272, .0295, .0311, .0111, .0633], [.098, .075, .074, .765, .425])


@pytest.mark.parametrize('width', [1e2, 1e8, 1e12])
@pytest.mark.parametrize('derivative', [False, True])
def test_overdamped_poles_remain_finite(width, derivative):
    energies = np.array([0., 7., np.nextafter(7., np.inf), 20., 600., 1e5])
    with np.errstate(all='raise'):
        result = causal.removed_real(energies, 8.65, width, 7., derivative)
    assert np.all(np.isfinite(result))


@pytest.mark.parametrize('derivative', [False, True])
@pytest.mark.parametrize('centre,width', [(8.65, 1.6), (16.3, 14.), (8., 16.),
                                         (8., 16.0001), (1000., 300.)])
@pytest.mark.parametrize('energy', [0., 1., 6., 7.00001, 20., 100.])
def test_removed_kk_against_independent_principal_value(derivative, centre, width, energy):
    def numerator(w):
        return 2/np.pi*w*causal.imaginary(w, centre, width, derivative)/(w+energy) if w+energy else 0.
    if 0 < energy < 7:
        expected = quad(numerator, 0., 7., weight='cauchy', wvar=energy,
                        epsabs=1e-11, epsrel=1e-10)[0]
    else:
        expected = quad(lambda w: numerator(w)/(w-energy), 0., 7.,
                        epsabs=1e-11, epsrel=1e-10)[0]
    actual = causal.removed_real(energy, centre, width, 7., derivative)
    assert float(actual) == pytest.approx(expected, rel=2e-8, abs=1e-11)


@pytest.mark.parametrize('phase', ['amorphous', 'hexagonal'])
@pytest.mark.parametrize('q', [0., 1., 8.])
@pytest.mark.parametrize('energy', [0., 8., 20., 100.])
def test_full_valence_kk(phase, q, energy):
    s, c = model.epsilon_optical(phase), coefficients()
    strength = model.oxygen_K_hydrogenic_gos_continuum_strength(q) if IS_ION else None
    fj, fi, ei, gj, gi = model._valence_parameters(s, np.array([q]), c, strength)
    def im(w):
        value = sum(fj[0,j]*causal.imaginary(w, o.E0, gj[0,j], True)
                    for j,o in enumerate(s.excitations))
        value += sum(fi[0,j]*causal.imaginary(w, ei[0,j], gi[0,j])
                     for j in range(len(s.ionizations)))
        return s.Ep**2*value
    upper = 1e5
    def smooth(w):
        return 2/np.pi*w*im(w)/(w+energy)
    if energy > s.Bmin:
        low = quad(smooth, s.Bmin, upper, weight='cauchy', wvar=energy,
                   epsabs=1e-8, epsrel=1e-9, limit=300)[0]
    else:
        low = quad(lambda w: smooth(w)/(w-energy), s.Bmin, upper,
                   epsabs=1e-8, epsrel=1e-9, limit=300)[0]
    tail = quad(lambda x: smooth(np.exp(x))/(np.exp(x)-energy)*np.exp(x),
                np.log(upper), np.log(1e12), epsabs=1e-10)[0]
    actual = model.epsilon1_valence_Eq([energy], q, s, c,
        inner_shell_strength_electrons=strength)['total'][0,0]
    assert actual == pytest.approx(1+low+tail, rel=2e-7, abs=2e-8)


@pytest.mark.parametrize('phase', ['amorphous', 'hexagonal'])
def test_integrated_elf_and_conductivity_budgets(phase):
    s, c = model.epsilon_optical(phase), coefficients()
    # Resolve the logarithmic hard-edge limit independently of the broad grid.
    edges = np.array([s.Bmin, *[o.Bth for o in s.ionizations], s.kshell.Bth])
    w = np.unique(np.r_[np.geomspace(s.Bmin, 1e9, 90001), edges,
        np.nextafter(edges, np.inf), s.Bmin+np.geomspace(1e-12, .1, 1000)])
    for q in [0., .5, 1., 2.5, 8., 30.]:
        strength = model.oxygen_K_hydrogenic_gos_continuum_strength(q) if IS_ION else None
        real = model.epsilon1_valence_Eq(w, q, s, c,
            inner_shell_strength_electrons=strength)['total'][0]
        im = model.epsilon2_valence_Eq(w, q, s, c,
            inner_shell_strength_electrons=strength)['total'][0]
        core = (model.oxygen_K_hydrogenic_gos_elf(w, q, Ep_eV=s.Ep) if IS_ION
                else model.epsilon2_Kshell_E0(w, s))
        scale = 2/(np.pi*s.Ep**2)
        conductivity = scale*np.trapezoid(w*(im+core), w)
        inverse = scale*np.trapezoid(w*(im/(real*real+im*im)+core), w)
        assert conductivity == pytest.approx(1., abs=3e-5)
        assert inverse == pytest.approx(1., abs=3e-5)


@pytest.mark.parametrize('phase', ['amorphous', 'hexagonal'])
def test_threshold_limits_sampling_independence_and_optical_parity(phase):
    s, c = model.epsilon_optical(phase), coefficients()
    w = np.array([50., 7., 25., 9., 12., 18., 1000., 50., 0.])
    actual = model.epsilon2_valence_Eq(w, [0., 1.], s, c)
    for i, value in enumerate(w):
        one = model.epsilon2_valence_Eq([value], [0., 1.], s, c)
        for key in ['excitations', 'ionizations']:
            np.testing.assert_allclose(np.array(actual[key])[...,i],
                                       np.array(one[key])[...,0], rtol=1e-14, atol=1e-15)
    assert np.all(actual['total'][:, w <= s.Bmin] == 0)
    assert np.all(np.isfinite(model.epsilon1_valence_Eq(w, [0., 1.], s, c)['total']))
    np.testing.assert_array_equal(optical.epsilon1_valence_E0(w,s)['total'],
                                  model.epsilon1_valence_Eq(w, 0., s, c)['total'][0])
    np.testing.assert_array_equal(optical.epsilon2_valence_E0(w,s)['total'],
                                  model.epsilon2_valence_Eq(w, 0., s, c)['total'][0])
