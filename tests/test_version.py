import pytest

from feeds.version import package_version

H = "d8e9f777" + "0" * 56


def test_version_takes_the_release_else_the_input_hash_and_counts_rebuilds() -> None:
    assert (
        package_version(upstream_release="2026-11-03", input_sha256=H, revision=1) == "2026-11-03"
    )
    assert package_version(upstream_release=None, input_sha256=H, revision=1) == "d8e9f777"
    assert package_version(upstream_release=None, input_sha256=H, revision=2) == "d8e9f777-r2"
    with pytest.raises(ValueError):
        package_version(upstream_release=None, input_sha256=H, revision=0)


@pytest.mark.parametrize("bad", ["D8E9" + "0" * 60, "abc", "g" * 64])
def test_version_refuses_an_input_hash_that_is_not_lowercase_hex(bad: str) -> None:
    with pytest.raises(ValueError):
        package_version(upstream_release=None, input_sha256=bad)


@pytest.mark.parametrize("release", ["..", ".", "a/b", "a b"])
def test_version_refuses_a_release_that_is_not_a_path_segment(release: str) -> None:
    with pytest.raises(ValueError):
        package_version(upstream_release=release, input_sha256=H)
