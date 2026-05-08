"""Focused tests for real-analysis prerequisite tissue selection."""

from argparse import Namespace

import pytest

from scripts.check_real_analysis_prereqs import POSTER_TISSUES, _selected_tissues


def test_selected_tissues_defaults_to_configured_tissues():
    cfg = {"gtex": {"tissues": ["Whole_Blood", "Spleen"]}}
    args = Namespace(tissues=None, poster_tissues=False)

    assert _selected_tissues(cfg, args) == ["Whole_Blood", "Spleen"]


def test_selected_tissues_uses_cli_override():
    cfg = {"gtex": {"tissues": ["Whole_Blood", "Spleen"]}}
    args = Namespace(tissues="Lung, Colon_Sigmoid", poster_tissues=False)

    assert _selected_tissues(cfg, args) == ["Lung", "Colon_Sigmoid"]


def test_selected_tissues_can_use_poster_subset():
    cfg = {"gtex": {"tissues": ["Whole_Blood", "Spleen"]}}
    args = Namespace(tissues=None, poster_tissues=True)

    assert _selected_tissues(cfg, args) == POSTER_TISSUES


def test_selected_tissues_requires_configured_tissues_when_not_overridden():
    args = Namespace(tissues=None, poster_tissues=False)

    with pytest.raises(ValueError, match="gtex.tissues"):
        _selected_tissues({}, args)