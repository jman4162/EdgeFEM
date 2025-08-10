"""Mapping CAD attributes to simulation tags."""

from __future__ import annotations

from pathlib import Path

from .. import Case


def apply_tags(case: Case, sidecar: str | Path | None = None) -> Case:
    """Populate :attr:`Case.meta` with tags.

    The implementation here is intentionally tiny: if ``sidecar`` is
    provided, the function records its path.  Real parsing of STEP
    attributes or YAML is left for future work.
    """

    if sidecar:
        case.meta.setdefault("sidecar", Path(sidecar))
    return case


__all__ = ["apply_tags"]
