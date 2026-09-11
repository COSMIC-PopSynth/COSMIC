"""Focused regression tests for modern-Python compatibility paths."""

import operator
import sys

import pytest

from cosmic import utils
from cosmic.filter import parse_column_filters
from cosmic.sample.sampler import sampler as sampler_registry
from cosmic.sample.stroopwafel.presets import merging_dco


def test_parse_column_filters_flattens_nested_definitions():
    preset = ("mass_2", operator.lt, 5.0)

    result = parse_column_filters(
        [["mass_1 >= 8"], [preset]],
        "ecc < 0.9",
    )

    assert result == [
        ("mass_1", operator.ge, 8.0),
        preset,
        ("ecc", operator.lt, 0.9),
    ]


def test_parse_inifile_preserves_constant_types(tmp_path):
    inifile = tmp_path / "constants.ini"
    inifile.write_text(
        """\
[sse]
stellar_engine = 'sse'

[bse]
integer = 7
floating = 2.5
enabled = True
missing = None
values = [1, 'two', False, None]
expression = 2 + 3 * 4

[rand_seed]
seed = 42

[filters]
select_final_state = True

[convergence]
pop_select = formation

[sampling]
sampling_method = independent
""",
        encoding="utf-8",
    )

    bse, sse, seed, filters, convergence, sampling = utils.parse_inifile(inifile)

    assert bse == {
        "integer": 7,
        "floating": 2.5,
        "enabled": True,
        "missing": None,
        "values": [1, "two", False, None],
        "expression": 14,
    }
    assert sse == {"stellar_engine": "sse"}
    assert seed == 42
    assert filters == {"select_final_state": True}
    assert convergence == {"pop_select": "formation"}
    assert sampling == {"sampling_method": "independent"}


def test_register_sampler_accepts_docless_method(monkeypatch):
    class DummyData:
        def sampler(self):
            pass

    def sample():
        return "sampled"

    monkeypatch.setattr(sampler_registry, "_SAMPLERS", {})

    sampler_registry.register_sampler("docless", DummyData, sample)

    registered = sampler_registry.get_sampler("docless", DummyData)
    assert registered is sample
    assert registered() == "sampled"
    assert DummyData.sampler.__doc__ is None


def test_merging_dco_names_the_optional_dependency_group(monkeypatch):
    # A None entry in sys.modules makes the import inside merging_dco fail the
    # same way it does for a user who never installed the optional extra.
    monkeypatch.setitem(sys.modules, "legwork", None)

    with pytest.raises(ImportError, match=r"pip install cosmic-popsynth\[merging-dco\]"):
        merging_dco(kstar_1=[14], kstar_2=[14])
