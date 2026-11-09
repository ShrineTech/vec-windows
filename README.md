# vec-windows: the external-data "protected window" rule (Rules s.10) as code

A community tool for the [Virtual Embryo Challenge](https://virtualembryo.ai/challenge). It is **not official**, and it is not affiliated with or endorsed by the organisers.

## What it is
You give it the per-cell stage labels (and, optionally, genotype labels) of an **external** dataset you want to train on. It tells you, for each label:

- **KEEP**: allowed;
- **KEEP_DISCLOSE**: allowed, but after E13.5, so the source and its stages must be stated explicitly;
- **DROP**: must be removed before training. Each DROP comes with the reason and the verbatim rule text.

It also writes a list of the cells to drop and a disclosure snippet for your method summary.

Rules encoded (Rules s.10, read 2026-10-19):

| Task | Protected window (cells inside are dropped) | Boundaries |
|---|---|---|
| T1 (extrapolation) | after E9.5 up to and including E13.5 | E9.5 kept, E13.5 dropped, E10.5 and E12.5 excluded absolutely, after E13.5 kept with disclosure |
| T2 heart | interpolation: strictly between E8.25 and E8.75; extrapolation: (E9.5, E13.5] | E8.25, E8.75 and E9.5 kept |
| T2 embryo | strictly between E7.25 and E8.0 | E7.25 and E8.0 kept |
| T3 (genotype) | Gata4 and beta-catenin (Ctnnb1) knockouts and other alleles of those genes | dropped at any stage (conservative) |

- Pass `--task` several times to get the union of the windows, for one model trained for several tasks.
- Somite counts (`25 somites`, `ss25`) and Theiler stages (`TS15`) are placed on the E-day scale as a **range**, as s.10 requires ("Fractional stagings [...] are placed on this scale and judged by the same boundary").
- A range that straddles a boundary is **dropped**.
- A label the tool cannot parse is **dropped** until you stage it with `--stage-map label,lo,hi`.

## Why it helps participants
- T1 expressly encourages external atlases.
- s.10 is the most intricate rule in the challenge: midpoint boundaries, an explicit extrapolation window with asymmetric ends, strict bracketing for interpolation, fractional stagings, and genotype phenocopies.
- Getting it wrong is a rules violation. This tool makes the filter reproducible and produces the disclosure text the Rules require.

## Related work
[rpinho/virtual-embryo-community](https://github.com/rpinho/virtual-embryo-community) (`ve-stagecheck`) triages **sample sheets** at specimen level against the stage windows, and sends Task 3 to manual condition review. vec-windows works one level down and covers the rest of s.10. It checks the **per-cell** stage and genotype labels of an AnnData obs table, so cell types inside one specimen are handled. It drops the held-out T3 knockouts (Gata4, Ctnnb1) automatically. It writes the drop list and the method-summary disclosure snippet. The two can be used together: ve-stagecheck to pick specimens, then vec-windows on the cells you keep.

## How to run
```bash
python vec_windows.py check examples/synthetic_obs.csv --stage-col stage --genotype-col genotype \
    --task T1 --source "MyAtlas v2 (URL, licence)" --drop-out drop.csv --disclosure disclosure.md
python vec_windows.py check atlas.h5ad --stage-col stage --task T2-heart --task T2-embryo   # needs anndata
```
- Exit code 4 means some cells must be dropped.
- `examples/synthetic_obs.csv` is made up and contains no real cells.
- For `.h5ad` input, only `obs` is read, in backed mode.

## Tests (`python -m unittest -v test_vec_windows`)
The tests use synthetic labels only. They cover:

- every boundary case written in the rule text (E9.5 versus E9.51, E13.5 versus E13.6, E8.25 and E8.75, E7.25 and E8.0);
- task unions;
- Theiler and somite ranges that straddle a boundary;
- genotype matching (Gata4, Ctnnb1, beta-catenin; Mab21l2 and WT are kept);
- unknown stages;
- the CLI drop list and disclosure snippet;
- `.h5ad` reading (skipped if `anndata` is absent).

## Known limits (read these)
- **The somite and Theiler mapping is APPROXIMATE and UNVERIFIED.**
  - Each Theiler stage is its nominal dpc (Theiler 1989, as tabulated by eMouseAtlas) widened to plus or minus 0.5 day, so only labels whose whole band clears a window are kept.
  - Somite counts are mapped through the same Theiler table.
  - Real embryos vary by strain and litter. If your source gives a better staging, use `--stage-map`.
- **Phenocopies cannot be detected from labels.** "a perturbation that phenocopies them" is your judgement; review your source.
- **Open question for T2 and T3.** The Data page says "Task 1 is the only task where training on external public single-cell data is allowed", while Rules s.10 permits external data for all tasks. The tool prints this warning for T2 and T3. Ask the organisers.
- An external stage that falls in another task's window (for example E7.5) is **kept** for T1, because s.10 is per task. If you reuse the same model for T2, pass both tasks.
- The Rules are the authority: "If a source sits near the boundary and you are unsure, ask before submitting."

## Data and AI disclosure
- It contains no challenge data. Its fixtures are synthetic labels.
- It was written with AI assistance (Claude Code, by Anthropic) and reviewed and filed by a human.
- Licence: MIT.
