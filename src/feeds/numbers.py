"""Number spelling shared by the recipes."""

from __future__ import annotations


def shortest(text: str) -> str:
    """Round-trip a published float to its shortest exact spelling.

    A workbook stores `0.61099999999999999`, the same IEEE-754 double as `0.611`,
    and 17 significant digits reads as false precision. `repr(float(x))` is the
    shortest string that round-trips to the identical double, so this is a
    rendering choice, not a rounding one. Raises `ValueError` on non-numeric text.
    """
    return repr(float(text))
