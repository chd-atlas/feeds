import pytest

from feeds.numbers import shortest


@pytest.mark.parametrize(
    ("published", "expected"),
    [
        ("0.61099999999999999", "0.611"),
        ("2.4500000000000002", "2.45"),
        ("0.29799999999999999", "0.298"),
        ("1E-3", "0.001"),
        # Exponents survive: a p-value of 1.21e-10 must not become 0.0.
        ("1.21e-10", "1.21e-10"),
        ("0", "0.0"),
    ],
)
def test_a_published_float_is_written_at_its_shortest_exact_spelling(
    published: str, expected: str
) -> None:
    assert shortest(published) == expected
    assert float(shortest(published)) == float(published)


def test_text_that_is_not_a_number_is_refused() -> None:
    with pytest.raises(ValueError):
        shortest("n/a")
