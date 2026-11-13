# Tools

The programs that **made the data**. They are provenance, not part of the
measurement pipeline, and nothing in `src/gspec/` imports them.

| file | what it does | state |
|---|---|---|
| `generate_synthetic.py` | builds `S01`–`S10` from seed `20261106` | verbatim from the original project; **verified**: `tests/test_generator.py` regenerates all ten FITS files bit for bit |
| `make_injections.py` | builds `R01`–`R06` (redacted) and `I01`–`I04` (with an exact increment) | verbatim; requires the raw archive files, which are not redistributed (see `download_sdss.py`) |
| `download_sdss.py` | fetches the ten raw SDSS spectra | **not executed** during this repository's validation — it needs network |
| `_compat.py` | the import shim that lets the two verbatim scripts run against `gspec` | — |

## Why "verbatim" matters

`generate_synthetic.py` and `make_injections.py` differ from the originals by
exactly two lines: the import of a sibling `common.py` became an import of
`_compat.py`. Rewriting them to the package's API would have meant that the
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
