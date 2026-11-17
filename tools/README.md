# Tools

The programs that **made the data**. They are provenance, not part of the
measurement pipeline, and nothing in `src/gspec/` imports them.

| file | what it does | state |
|---|---|---|
| `generate_synthetic.py` | builds `S01`–`S10` from seed `20261106` | verbatim; **verified**: `tests/test_generator.py` regenerates all ten FITS files bit for bit |
| `make_injections.py` | builds `R01`–`R06` (redacted) and `I01`–`I04` (with an exact increment) | verbatim; requires the raw archive files, which are not redistributed (the retrieval workflow that fetches them is `download_sdss.py`) |
| `_compat.py` | the import shim that lets the two verbatim scripts run against `gspec` | — |
| `download_sdss.py` | fetches a spectrum by plate/MJD/fibre, verifies it, converts it | **v0.3**; network; see `docs/sdss-tutorial.md` |
| `fits_to_case.py` | converts an archive FITS into the case schema | **v0.3**; the converter, callable on its own |

## Why "verbatim" matters

`generate_synthetic.py` and `make_injections.py` differ from the versions that
produced `data/` by two imports — the import of a sibling `common.py` became an
import of `_compat.py` — and, in `make_injections.py`, one docstring sentence and one
code comment. Rewriting them to the package's API would have meant that the
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


## The V0.3 retrieval tools

These two are not provenance — they are the workflow a reader uses to measure a
spectrum that is not in `data/`. Both are standalone (stdlib + `numpy` + `astropy`);
neither is imported by `src/gspec/`, and neither touches the fitter.

```bash
python tools/download_sdss.py --plate 294 --mjd 51986 --fiber 36 --out data/raw
gspec measure --cases data/raw/cases --out data/raw/measurements.csv
```

Separating them is deliberate. `fits_to_case.py` can be pointed at any FITS file
and reports what is wrong with it without downloading anything, which is what makes
it testable offline; `download_sdss.py` owns the network and the provenance and
calls the converter only after a file has passed its own checks.

`data/raw/` is git-ignored: the archive files carry `SPECOBJ` and `SPZLINE`, the
survey pipeline's own redshift and line fluxes, and are never redistributed
(DATA.md §3).
