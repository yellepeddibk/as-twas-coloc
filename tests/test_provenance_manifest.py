"""
test_provenance_manifest.py
---------------------------
Focused tests for provenance manifest generation.
"""

import src.provenance as provenance


def test_build_manifest_records_external_tool_versions(monkeypatch):
    """Manifest records discovered external tool versions instead of hardcoded nulls."""
    monkeypatch.setattr(provenance, "_try_version", lambda package: f"{package}-version")
    monkeypatch.setattr(provenance, "_metaxcan_version", lambda cfg: "git:9b9f76b")
    monkeypatch.setattr(provenance, "_rscript_version", lambda: "Rscript (R) version 4.5.3")
    monkeypatch.setattr(provenance, "_r_package_version", lambda package: "6.0.0")

    manifest = provenance.build_manifest(
        {
            "coloc": {},
            "correction": {},
            "gwas": {},
            "spredixcan": {},
        }
    )

    assert manifest["tool_versions"]["MetaXcan/SPrediXcan"] == "git:9b9f76b"
    assert manifest["tool_versions"]["R"] == "Rscript (R) version 4.5.3"
    assert manifest["tool_versions"]["coloc_R_package"] == "6.0.0"