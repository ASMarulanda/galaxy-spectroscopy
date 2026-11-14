# Tools

The programs that **made the data**. They are provenance, not part of the
measurement pipeline, and nothing in `src/gspec/` imports them.

| file | what it does | state |
|---|---|---|
| `generate_synthetic.py` | builds `S01`–`S10` from seed `20261106` | verbatim; **verified**: `tests/test_generator.py` regenerates all ten FITS files bit for bit |
| `make_injections.py` | builds `R01`–`R06` (redacted) and `I01`–`I04` (with an exact increment) | verbatim; requires the raw archive files, which are not redistributed (a download workflow ships in v0.3) |
| `_compat.py` | the import shim that lets the two verbatim scripts run against `gspec` | — |

## Why "verbatim" matters

`generate_synthetic.py` and `make_injections.py` differ from the versions that
produced `data/` by two imports — the import of a sibling `common.py` became an
import of `_compat.py` — and, in `make_injections.py`, one docstring sentence. Rewriting them to the package's API would have meant that the
generator shipped here is not the generator that produced `data/cases/`, and the
reproducibility claim for the only exactly-known truth in the project would rest
on a program that never ran.

The shim supplies the names, not the algorithms. If it were wrong, the bit-for-bit
test would fail — it does not.

## Regenerating the synthetic data

```bash
.venv/bin/python tools/generate_synthetic.py     # writes ../cases and ../truth_synthetic.json
```

Note the script's output paths are relative to `tools/`, inherited unchanged; run
it and move the results, rather than editing the paths — the paths are part of
what is being preserved.
