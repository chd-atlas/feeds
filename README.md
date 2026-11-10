# feeds

Versioned source packages for the [CHD atlas](https://chd-atlas.github.io/).

Each folder under `sources/` reads one upstream source into the atlas's table
format and publishes it as a GitHub release `<source_id>/<version>`: a
`manifest.json` and one gzipped TSV per table. The atlas pins a package in its
`mirrors/sources.yaml` and ingests it with `chd-atlas ingest`; the atlas is the
only thing that can refuse data.
