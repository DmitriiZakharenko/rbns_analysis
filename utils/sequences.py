"""Sequence helpers: T->U, validation, reverse-complement."""

import re

RNA_BASES = set("ACGU")
DNA_BASES = set("ACGT")


def to_rna(seq):
    """Convert DNA to RNA (T -> U)."""
    return seq.upper().replace("T", "U")


def is_valid_rna(seq, allowed_lengths=None):
    """Check sequence is non-empty, only A/C/G/U, and optionally length in allowed_lengths."""
    if not seq or not isinstance(seq, str):
        return False
    s = seq.upper().replace("T", "U")
    if not set(s) <= RNA_BASES:
        return False
    if allowed_lengths is not None and len(s) not in allowed_lengths:
        return False
    return True


def reverse_complement(seq, dna=True):
    """Return reverse complement. If dna=True use T, else U."""
    comp = {"A": "T", "T": "U" if not dna else "T", "U": "A", "G": "C", "C": "G", "N": "N"}
    seq = seq.upper().replace("U", "T") if dna else seq.upper()
    return "".join(comp.get(b, b) for b in reversed(seq))
