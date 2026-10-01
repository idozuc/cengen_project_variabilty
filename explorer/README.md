# Explorer

`build_assets.py` converts the current BASiCS master and supported WormCat results
into versioned, compact assets. The R clustering exporter writes the companion
`clustering/` assets. `app.py` reads only those two current contracts.

Use `CENGEN_OUTPUT_DIR` to override the default `outputs/explorer` directory.
See `schemas/` and `docs/output-schema.md` for interface details.
