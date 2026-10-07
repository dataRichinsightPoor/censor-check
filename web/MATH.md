# Censor Check: mathematical contract

Data-Rich, Insight-Poor · 99 Small Problems · No. 06 · v0.1.0-alpha
Companion to The Audit, No. 1, "Where the Inactives Went."

The question is narrow. A potency table records some compounds as exact IC50 values and others as bounds. What does a filter that keeps only exact values remove? How large is the resulting bias? And does a comparison made by pooling many papers survive inside the papers that measured both arms? This contract states what the tool computes, which assumptions produce each number, and where the numbers stop meaning anything.

## Records and the filter

Each record \(i\) carries a relation \(r_i\), a concentration \(v_i\) converted to nM, a group \(g_i\) (the paper), a stratum \(s_i\) (usually the target), an optional compound identifier, and a validity flag. Relations are normalized to four classes:

- exact: "=";
- right-censored: ">", "≥", "≫". The curve was not crossed inside the tested range, so the true IC50 exceeds \(v_i\);
- left-censored: "<", "≤". Counted, not modeled;
- other or missing.

The exact-value filter mirrors the documented pChEMBL rule: relation "=", units convertible to nM, \(v_i>0\), and a validity comment that is empty or "manually validated." The analysis set \(\mathcal{A}\) is the union of records kept by the filter and right-censored records with a positive bound.

Potency is expressed as

\[
p_i = 9-\log_{10} v_i ,
\]

so a bound "\(\mathrm{IC50}>c\) nM" becomes "\(p_i<p_c\)", with \(p_c=9-\log_{10}c\).

## Weak tail

For a threshold \(T\) (default 10,000 nM, \(p_T=5\)), a record is weak if \(v_i\ge T\), whether \(v_i\) is an exact value or a bound. A bound below \(T\) cannot say whether a compound is weak and is counted as not weak, which makes every loss estimate conservative.

- Record-weighted weak share, all evidence: \(W_{\mathrm{all}}=|\{i\in\mathcal{A}:v_i\ge T\}|/|\mathcal{A}|\).
- After the filter: the same ratio restricted to exact records.
- Censored share of the weak evidence: weak bounds divided by all weak records.
- Paper-weighted shares: the mean, over paper–stratum sets with at least \(n_{\min}\) records, of each set's own weak share.
- Reachable shares: restricted to paper–stratum sets whose largest reported value or bound is at least \(T\). The largest reported value is a lower bound on the top tested concentration, so a set that never reports anything at or above \(T\) could not have recorded weakness at \(T\).
- Tail lost: sets with at least one weak bound and no weak exact value. After filtering, such a set appears to contain only active compounds.

## Pooled censored normal

Assume all potencies in \(\mathcal{A}\) are draws from \(N(\mu,\sigma^2)\). Exact values contribute their density and bounds contribute the probability of lying below them:

\[
\ell(\mu,\sigma)=\sum_{i\in E}\left[-\log\sigma-\tfrac12\log 2\pi-\frac{(p_i-\mu)^2}{2\sigma^2}\right]+\sum_{j\in C}\log\Phi\!\left(\frac{p_{c_j}-\mu}{\sigma}\right).
\]

This is Tobin's limited-dependent-variable likelihood with record-specific censoring points. It is maximized over \((\mu,\log\sigma)\) by Nelder–Mead from the exact-only mean and standard deviation. The tool reports three central estimates:

- exact only: \(\bar p_E\), what a filtered dataset reports;
- substitution: the mean of exact values and bounds taken at face value;
- censored: \(\hat\mu\).

It also reports the overstatement \(\bar p_E-\hat\mu\) in log units and as a fold change \(10^{\bar p_E-\hat\mu}\), and the predicted weak share \(\Phi((p_T-\hat\mu)/\hat\sigma)\).

Because each \(\log\Phi\) term is decreasing in \(\mu\), adding bounds can only move \(\hat\mu\) toward weaker potency. Substitution sits between the other two estimates. It treats each bound as the weakest value consistent with it, so it still overstates potency.

## Grouped censored normal

Medicinal-chemistry papers are not exchangeable draws from one distribution. Each reports a congeneric series designed around a scaffold, a hypothesis and an assay. The grouped model gives each paper–stratum set \(g\) its own mean and makes all sets share one within-set spread:

\[
p_i\sim N(\mu_{g_i},\sigma_w^2),\qquad \ell=\sum_g\left[\sum_{i\in E_g}\log\tfrac{1}{\sigma_w}\varphi\!\left(\tfrac{p_i-\mu_g}{\sigma_w}\right)+\sum_{j\in C_g}\log\Phi\!\left(\tfrac{p_{c_j}-\mu_g}{\sigma_w}\right)\right].
\]

Sets with no exact value are excluded, because their mean is unbounded below. The likelihood is maximized by coordinate ascent. Each \(\mu_g\) is updated by golden-section search on \([\bar p_{E_g}-6,\ \bar p_{E_g}+3]\), where the profile is concave. Then \(\log\sigma_w\) is updated by golden-section search on \([-3,2]\). Iteration stops when \(\sigma_w\) changes by less than \(10^{-6}\), or after 30 rounds.

The tool reports:

- the within-set SD \(\hat\sigma_w\);
- the between-set SD, the record-weighted standard deviation of \(\hat\mu_g\);
- the share of variance between sets, \(\sigma_b^2/(\sigma_b^2+\sigma_w^2)\);
- for sets with at least one bound and at least \(m\) exact values, the per-set overstatement \(\bar p_{E_g}-\hat\mu_g\), as a median and interquartile range.

The per-set overstatement is non-negative by construction. Its magnitude is informative; its sign is not.

When between-set variance dominates, a pooled single-distribution model is misspecified. Its \(\hat\sigma\) mixes chemistry with paper, and its \(\hat\mu\) depends on how many records each paper contributes.

## Pooled versus within-paper comparison

For strata \(A\) and \(B\), let \(w_{A}\) and \(w_{B}\) be weak shares. The tool reports three differences:

- pooled: \(w_A-w_B\) over every record;
- pooled over shared papers: the same difference, restricted to papers that measured both;
- within-paper: a Mantel–Haenszel-weighted mean of per-paper differences,

\[
\Delta_{\mathrm{within}}=\frac{\sum_k \omega_k\,(w_{A,k}-w_{B,k})}{\sum_k\omega_k},\qquad \omega_k=\frac{n_{A,k}\,n_{B,k}}{n_{A,k}+n_{B,k}} .
\]

A reversal is flagged when the within-paper difference has the opposite sign to the pooled difference and is at least one percentage point in magnitude. The same contrast is computed for the mean exact-value pIC50: the pooled difference, and the median of per-paper differences over papers with exact values in both strata.

The two differences estimate different quantities. The pooled difference describes the literature's composition: which series were published for each target, at what scale, and with which top concentrations. The within-paper difference holds the series fixed. In a selectivity panel, though, the off-target is weak by design, so the within-paper difference measures intended selectivity rather than the intrinsic difficulty of either target. Neither difference is a property of the targets alone.

## Counter-screens

In papers that measured two or more strata, the stratum with the highest median pIC50 (bounds taken at face value) is labeled primary, and the others are labeled counter-screens. The tool reports the censored share of records in each role and in single-stratum papers. The label is operational. It does not read the paper's stated intent.

## Disguised ceilings

If a paper recorded "inactive at the top concentration" as an exact value equal to that concentration, several compounds would share its highest exact value. For paper–stratum sets with at least 8 exact values, the tool reports:

- the share of sets in which the highest exact value is tied by at least 3 records;
- the expected share for a randomly chosen interior record: the mean over sets of the fraction of non-maximal exact records whose value is tied by at least 3 records;
- the number of sets in which the highest exact value also appears as a bound.

A top-tie rate well above the interior rate is the signature of disguised censoring.

## Compound level

Records are collapsed to compound–stratum pairs. The tool reports pairs known only through bounds, which vanish entirely under the filter, and pairs with both exact and censored records. Records are not deduplicated in any other calculation.

## Numerics

\(\log\Phi\) is computed from a log-space Chebyshev approximation to erfc, with fractional error below \(1.2\times10^{-7}\). It is finite beyond \(z=-40\). All calculations run in the browser. Uploaded files are not transmitted.

## Boundaries

- The tool reads relations and values. It cannot inspect curves, so it cannot tell whether an exact value was truly bracketed.
- Censoring points are assumed known and independent of potency given the paper. In practice, top concentrations are chosen by laboratories, often because of solubility, aggregation or cytotoxicity limits. That makes censoring informative about the chemistry.
- Normality within a set is an assumption. Congeneric series can be skewed or multimodal.
- Left-censored records are counted, not modeled.
- Estimates describe a dataset. They are not corrected potencies for any compound.
- Software verification is not biological validation.

## References

1. Tobin J. Estimation of relationships for limited dependent variables. Econometrica. 1958;26(1):24–36. https://doi.org/10.2307/1907382
2. Mantel N, Haenszel W. Statistical aspects of the analysis of data from retrospective studies of disease. Journal of the National Cancer Institute. 1959;22(4):719–748. https://doi.org/10.1093/jnci/22.4.719
3. Zdrazil B, et al. The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods. Nucleic Acids Research. 2024;52(D1):D1180–D1192. https://doi.org/10.1093/nar/gkad1004
4. Heinzke AL, et al. A compound-target pairs dataset: differences between drugs, clinical candidates and other bioactive compounds. Scientific Data. 2024;11:1160. https://doi.org/10.1038/s41597-024-03582-9
5. Landrum GA, Riniker S. Combining IC50 or Ki values from different sources is a source of significant noise. Journal of Chemical Information and Modeling. 2024;64(5):1560–1567. https://doi.org/10.1021/acs.jcim.4c00049
6. Svensson E, et al. Enhancing uncertainty quantification in drug discovery with censored regression labels. Artificial Intelligence in the Life Sciences. 2025;7:100128. https://doi.org/10.1016/j.ailsci.2025.100128
7. Simpson EH. The interpretation of interaction in contingency tables. Journal of the Royal Statistical Society, Series B. 1951;13(2):238–241. https://doi.org/10.1111/j.2517-6161.1951.tb00088.x
