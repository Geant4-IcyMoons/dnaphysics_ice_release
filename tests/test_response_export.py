"""Electron generator provenance and a small numerical DAT round trip."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'python_scripts'/'physics_ice'))
import generate_ice_cross_sections as gen


def test_old_cache_rejected_and_current_cache_roundtrips(tmp_path):
    assert not gen._npz_matches_params(dict(NE=10, Nq=12), 10, 12)
    path = tmp_path/'cross_sections.npz'
    gen.save_cross_section_corrections_npz([1e6], [{}], path, NE=10, Nq=12)
    with np.load(path) as data:
        assert gen._npz_matches_params(data, 10, 12)
        assert not gen._npz_matches_params(data, 11, 12)
        damaged = dict(data, dielectric_response_signature='old')
        assert not gen._npz_matches_params(damaged, 10, 12)


def test_legacy_dcs_rejected(tmp_path):
    with pytest.raises(ValueError, match='provenance'):
        gen._write_dcs_tables_from_data({}, tmp_path/'exc.dat', tmp_path/'ion.dat')


def test_dcs_tcs_export_uses_corrected_response_and_numeric_format(tmp_path, monkeypatch):
    model = gen.model
    s = model.epsilon_optical('amorphous')
    c = model.DispersionCoeffs([3.82,2.47,2.47,3.01,2.44],
        [.0272,.0295,.0311,.0111,.0633], [.098,.075,.074,.765,.425])
    w = np.array([7., 8., 10., 13., 17., 32., 40.])
    d = np.array([[gen._integrate_channel_single_E(x, 100., j, 'ionization', s, c, 32)
                   for j in range(4)] for x in w])
    assert np.all(np.isfinite(d)) and np.all(d >= 0)
    data = dict(dielectric_response_signature=gen._response_signature(),
        T_line=np.full(w.shape,100.), E_line=w, exc_vals=np.zeros((w.size,5)),
        ion_vals=np.column_stack([d/gen.EMFI_DCS_SCALE_M2, np.zeros(w.size)]))
    monkeypatch.setattr(gen, '_export_to_custom_geant4', lambda paths: None)
    gen.write_emfietzoglou_dcs_tables(s,c,out_dir=tmp_path,dcs_data=data,
        parallel_channels=False,ice_label='test')
    paths=list(tmp_path.glob('*.dat'))
    assert len(paths)==4
    assert all(np.all(np.isfinite(np.loadtxt(p))) for p in paths)
    assert len(list(tmp_path.glob('dielectric_response_*.csv')))==1
    assert not list(tmp_path.glob('*.tmp'))
