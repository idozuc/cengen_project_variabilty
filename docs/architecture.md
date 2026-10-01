# Architecture

```text
SingleCellExperiment
├── clustering runner
│   ├── screen and resolution diagnostics
│   ├── accepted markers and visualization bundle
│   └── clustering explorer exporter
└── BASiCS runner
    ├── fit, diagnostics, LVG, experiment support
    ├── canonical master table
    └── GO/WormCat enrichment
        └── explorer asset builder

clustering assets + BASiCS assets → Streamlit explorer
```

R reusable functions live in `src/R`; executable orchestration lives in
`scripts`. Scheduler files do not contain scientific logic. Every run is
self-describing through resolved configuration, checksums, and session data.
