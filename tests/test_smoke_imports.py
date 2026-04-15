"""Smoke tests for package-level imports and public exports."""

import HIREXSR_py as pkg


def test_package_imports_with_version() -> None:
    assert hasattr(pkg, "__version__")
    assert isinstance(pkg.__version__, str)
    assert pkg.__version__


def test_public_entry_points_exposed() -> None:
    # Validate key package-level exports from __init__.py remain available.
    assert callable(pkg.hirexsr_get_profile_py)
    assert callable(pkg.hirexsr_get_lint_profile_py)
    assert callable(pkg.hirexsr_load_result_py)
    assert callable(pkg.quality_check_lint)
    assert callable(pkg.quality_check_profile)
    assert callable(pkg.zeff_neo)
