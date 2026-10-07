# Censor Check v0.1.0-alpha: verification record

## Automated

`npm test` passes 23 of 23 tests (Node test runner, Node 20). The tests cover:

- **Numerics:** normal-CDF accuracy at reference quantiles and in both tails (finite beyond z = −40), erfc reference values, and the golden-section and Nelder–Mead optimizers.
- **Parsing:** quoted, semicolon-delimited ChEMBL exports with quoted relations; column mapping; relation and unit normalization; a readable error when required columns are missing.
- **Filter accounting:** the pChEMBL-style rule and the analysis set.
- **Pooled censored normal:** recovery of the sample mean and SD, including the censored part, from 20,000 synthetic draws censored at 10 µM. Exact-only estimates overstate potency, and substitution falls between the two. With no bounds, the censored and exact-only estimates are identical.
- **Grouped censored normal:** recovery of within-paper and between-paper SD from 150 synthetic papers, inflation of the pooled SD, and non-negative shifts.
- **Weak tail:** threshold, reachability and bound-position arithmetic on a hand-computed table.
- **Comparisons:** a constructed Simpson reversal with an exact within-paper difference of −10 points, and the identity between the within-paper and pooled differences for a single shared paper.
- **Role and tie checks:** counter-screen role assignment, and the disguised-ceiling test on flagged and clean synthetic papers.
- **Compound collapse, CSV labeling, and version consistency** across the package, the engine and both pages.
- **Regression:** the bundled ChEMBL subset reproduces every number in the companion article. This covers accounting, weak shares (record-weighted, paper-weighted and reachable), bound positions, pooled and grouped estimates, ties, compound pairs, counter-screens, the ABL1–AKT1 and AKT1–VEGFR2 comparisons, and the VEGFR2 stratum.

The JavaScript engine was also checked against an independent Python implementation (pandas and SciPy) of the same definitions. Shared quantities agree to the reported precision.

## Browser

Headless Chromium (Playwright), at 1360 × 900 in dark and light themes and at 375 × 800.

- The bundled subset loads and runs. The headline, the metrics and every section render.
- No horizontal page overflow at any size, and no page errors.
- Changing strata recomputes the comparison. The reachability toggle redraws the bars.
- A negative threshold is rejected with a readable message.
- Summary CSV and results JSON downloads fire.
- A semicolon-delimited ChEMBL-style export with quoted relations and µM units uploads, maps all six columns and runs.
- The methods page renders 58 MathJax expressions.

## Corrections made during QA

- A synthetic-recovery test first compared the fitted mean with the generating parameter. The fit matched the sample mean to 0.002 log units; the gap reflected sampling. The test now checks against the sample and, more loosely, against the generating parameters.
- The per-paper shift distribution was first described as evidence because every shift was positive. A positive shift is guaranteed by the likelihood, so the interface and the contract now state that only the magnitude is informative.

## Extended analyses (added after the first release)

`tests/depth.test.js` adds three regression tests: the threshold gradient, the ABL1–AKT1 split into selection and weighting terms, and the correlation between the truncated-normal prediction and the fitted per-paper shifts (r > 0.9). The correlation is partly built in, because the formula approximates the same model. The test guards the implementation; it does not validate the model independently. Bootstrap intervals come from `analysis/depth.mjs` with a fixed seed.
