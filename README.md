# Censor Check

99 Small Problems: Useful models for assumptions with expensive ambitions.  
No. 06 | Data-Rich, Insight-Poor · The Audit, No. 1 | v0.1.0-alpha

Companion article: Data-Rich, Insight-Poor · The Audit, No. 1, "Where the Inactives Went."

[Model](https://datarichinsightpoor.github.io/censor-check/) · [Code](https://github.com/dataRichinsightPoor/censor-check) · [Math](https://datarichinsightpoor.github.io/censor-check/methods.html) · [Math contract](web/MATH.md)

A greater-than sign is a measurement. Does your dataset keep it? Censor Check reads an IC50 table, applies the exact-value rule that defines a pChEMBL value, and reports what the rule removes. It then asks whether the conclusions you would draw from the pooled table survive inside the papers that produced it.

The bundled example covers 31,236 journal-literature IC50 records from ChEMBL for six targets: VEGFR2, ABL1, AKT1, BACE1, thrombin and HDAC1.

- **The weak tail is mostly censored.** In these records, 58.8% of the evidence for compounds weaker than 10 µM is a bound. The exact-value filter removes it, and the share of weak records falls from 17.8% to 8.4%. In 235 of 1,347 paper–target sets, the filter removes every compound the authors reported as weak.
- **The paper explains more of the spread than the compound.** Under a grouped censored model, the between-paper SD of mean pIC50 is 1.13 and the within-paper SD is 0.91. Dropping bounds overstates a paper's mean potency by a median of 0.30 log units, with an interquartile range of 0.15 to 0.53.
- **A pooled comparison can lack a stable sign while the within-paper one is robust.** Pooled, ABL1 records are 4.5 points more often weak than AKT1 records (paper-resampled 95% interval −3 to +15). Inside the 50 papers that measured both, ABL1 records are 19.7 points less often weak (−32 to −9). Most of the gap comes from which papers measure both targets, not from reweighting within them.
- **Counter-screens are censored more.** Counter-screen targets are bounds 29.6% of the time, against 13.3% for the target each series was optimized for. The filter removes selectivity evidence first.

## Robustness analyses

`node analysis/depth.mjs` reproduces the article's extended analyses: paper-resampled 95% intervals, threshold sensitivity (1 to 100 µM), era trends, the selection and weighting decomposition of pooled-versus-within comparisons, the association between each paper's bound position and its series potency, and the truncated-normal prediction of per-paper bias, \(\Delta(c)=\sigma_w\,\varphi(z_c)/(1-c)\). Results are written to `analysis/depth-results.json` and `analysis/per-paper-shifts.csv`.

## Open the model

- **[Try Censor Check](https://datarichinsightpoor.github.io/censor-check/):** no installation or sign-in. Load the bundled subset or your own CSV or ChEMBL export. Files stay in your browser.
- **[Read the mathematics](https://datarichinsightpoor.github.io/censor-check/methods.html):** likelihoods, estimators, weighting, boundaries and references. The same contract is in [web/MATH.md](web/MATH.md).
- **[Inspect the engine](web/model.js):** original dependency-free JavaScript.

## What it computes

- **Filter accounting:** exact, right-censored, left-censored and validity-flagged records, and what the pChEMBL-style rule keeps.
- **Weak-tail shares:** record-weighted, paper-weighted, and restricted to papers whose reported values reach the threshold. It also reports the censored share of weak evidence and the papers whose entire weak tail disappears.
- **Where bounds sit:** the share of bounds below, at and above the threshold. Bounds below the threshold are uninformative about weakness.
- **A pooled censored normal (Tobit) model:** with exact-only and substitution estimates beside it, the overstatement in log units and fold, and the predicted weak share.
- **A grouped censored normal model:** a per-paper mean with a shared within-paper SD, the between-paper SD, the variance share, and the distribution of per-paper overstatement.
- **Pooled versus within-paper comparison of two strata:** a Mantel–Haenszel-weighted within-paper difference, with sign-reversal flagging, for the weak share and the mean pIC50.
- **Counter-screen censoring:** in multi-target papers, the censored share by role.
- **A disguised-ceiling test:** ties at each paper's top exact value compared with the exact expected tie rate for interior values, and ceilings written both ways.
- **Compound–target pairs** known only through bounds.
- **Exports:** a summary CSV, results as JSON, and the figure as PNG.

## Files

- `web/index.html`, `web/app.js`, `web/style.css`: the interface.
- `web/model.js`: the engine. Pure functions, no dependencies.
- `web/MATH.md` and `web/methods.html`: the mathematical contract. `tools-build-methods.py` builds the HTML from the Markdown.
- `web/data/`: the bundled ChEMBL subset and its license and provenance note.
- `tests/model.test.js` (20 tests) and `tests/depth.test.js` (3 tests). They cover:
  - normal-CDF accuracy in both tails and the optimizers;
  - parsing of ChEMBL exports, and relation and unit normalization;
  - filter accounting;
  - recovery of known parameters by the pooled and grouped censored models;
  - weak-tail and reachability arithmetic;
  - a constructed Simpson reversal, and the single-paper identity;
  - counter-screen roles and the disguised-ceiling test;
  - compound collapse, CSV labeling and version consistency;
  - a regression test that reproduces every primary number published in the companion article from the bundled data;
  - the threshold gradient, the selection and weighting decomposition, and the truncated-normal prediction of per-paper bias.

Run the tests with Node 20 or later:

```
npm test
```

Serve `web/` with any static server to use the tool locally. Verification details are in `VERIFICATION.md`.

## Boundaries

The tool reads relations and values; it does not read dose-response curves, so it cannot tell whether an exact value was truly bracketed. The censored models assume normal potencies within a group and known censoring points. Laboratories choose top concentrations, often for solubility, aggregation or cytotoxicity reasons, so censoring can carry chemical information the model ignores. Left-censored records are counted, not modeled. The primary and counter-screen labels are operational. Estimates describe a dataset, not the potency of any compound. Software verification is not biological validation. The code is original; the public foundations are cited in `web/MATH.md`.

## License

Code: MIT. See `LICENSE`. Bundled data: derived from ChEMBL (EMBL-EBI) under CC BY-SA 3.0. See `web/data/README.md`.
