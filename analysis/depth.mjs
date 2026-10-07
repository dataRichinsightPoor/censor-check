// Robustness and decomposition analyses for The Audit, No. 1. Run: node analysis/depth.mjs
import fs from 'node:fs';
import * as M from '../web/model.js';
const txt = fs.readFileSync(new URL('../web/data/chembl-six-targets-journal-ic50.csv', import.meta.url), 'utf8');
const p = M.parseDelimited(txt); const map = M.mapColumns(p.header);
const yi = p.header.indexOf('year');
const recs = M.toRecords(p, map).map((r, i) => ({ ...r, year: +p.rows[i][yi] }));
const set = M.analysisSet(recs);
const out = {};
const median = a => M.quantile([...a].sort((u, v) => u - v), 0.5);
const mean = a => a.reduce((s, x) => s + x, 0) / a.length;
const T = 10000, weak = r => r.v >= T;
const pI = v => 9 - Math.log10(v);
// seeded RNG
let seed = 20261007; const rnd = () => (seed = (seed * 1664525 + 1013904223) >>> 0) / 4294967296;
const byDoc = new Map(); for (const r of set) { if (!byDoc.has(r.group)) byDoc.set(r.group, []); byDoc.get(r.group).push(r); }
const docs = [...byDoc.keys()];
const resample = () => { const s = []; for (let i = 0; i < docs.length; i++) { const d = docs[Math.floor(rnd() * docs.length)]; const g = byDoc.get(d); const tag = d + '#' + i; for (const r of g) s.push({ ...r, group: tag }); } return s; };
const q = (a, x) => M.quantile([...a].sort((u, v) => u - v), x);
const ci = a => [q(a, 0.025), q(a, 0.975)].map(x => +x.toFixed(2));

// 1. Cluster (paper) bootstrap
const G = M.groupedEstimates(set);
const stat = s => {
  const w = M.weakTail(s, T); const c1 = M.compareStrata(s, 'ABL1', 'AKT1'); const c2 = M.compareStrata(s, 'AKT1', 'VEGFR2'); const c3 = M.compareStrata(s, 'ABL1', 'VEGFR2'); const cs = M.counterScreens(s);
  return { weakCens: w.weakCensoredPct, weakAll: w.weakAllPct, weakKept: w.weakKeptPct, abl_akt_pooled: c1.pooledWeakDiff, abl_akt_within: c1.withinWeakDiff, akt_vegf_pooled: c2.pooledWeakDiff, akt_vegf_within: c2.withinWeakDiff, abl_vegf_pooled: c3.pooledWeakDiff, abl_vegf_within: c3.withinWeakDiff, counterMinusPrimary: cs.counterCensoredPct - cs.primaryCensoredPct };
};
const B = 400, boots = [];
for (let b = 0; b < B; b++) boots.push(stat(resample()));
out.point = stat(set);
out.bootstrapCI = Object.fromEntries(Object.keys(boots[0]).map(k => [k, ci(boots.map(x => x[k]))]));
const sh = G.shifts; const med = []; for (let b = 0; b < 2000; b++) { const s = []; for (let i = 0; i < sh.length; i++) s.push(sh[Math.floor(rnd() * sh.length)]); med.push(median(s)); }
out.shiftMedian = { point: +G.shiftMedian.toFixed(3), ci: ci(med), n: sh.length };
out.reversalBootstrapSameSignPct = 100 * boots.filter(x => x.abl_akt_within < 0 && x.abl_akt_pooled > 0).length / B;

// 2. Threshold sensitivity
out.threshold = [1000, 3000, 10000, 30000, 100000].map(t => { const w = M.weakTail(set, t); return { uM: t / 1000, weakAll: +w.weakAllPct.toFixed(1), weakKept: +w.weakKeptPct.toFixed(1), weakCens: +w.weakCensoredPct.toFixed(1), reachable: +w.reachableGroupsPct.toFixed(1) }; });

// 3. Era
const eras = [[1976, 1999], [2000, 2009], [2010, 2017], [2018, 2025]];
out.era = eras.map(([a, b]) => { const s = set.filter(r => r.year >= a && r.year <= b); const w = M.weakTail(s, T); const bp = M.boundPositions(s, T);
  return { era: `${a}-${b}`, n: s.length, boundPct: +(100 * s.filter(M.isBound).length / s.length).toFixed(1), weakAll: +w.weakAllPct.toFixed(1), weakKept: +w.weakKeptPct.toFixed(1), weakCens: +w.weakCensoredPct.toFixed(1), reachable: +w.reachableGroupsPct.toFixed(1), boundsAbovePct: +bp.abovePct.toFixed(1) }; });

// 4. Decomposition of pooled vs within: all -> shared papers pooled -> within (MH)
const dec = (A, Bn) => { const c = M.compareStrata(set, A, Bn); const a = set.filter(r => r.stratum === A), b = set.filter(r => r.stratum === Bn);
  const shA = new Set(a.map(r => r.group)), shB = new Set(b.map(r => r.group));
  const wsh = (x, other, inShared) => { const y = x.filter(r => other.has(r.group) === inShared); return { n: y.length, weak: +(100 * y.filter(weak).length / y.length).toFixed(1) }; };
  return { pair: `${A}-${Bn}`, pooled: +c.pooledWeakDiff.toFixed(1), sharedPooled: +c.sharedPooledWeakDiff.toFixed(1), within: +c.withinWeakDiff.toFixed(1), selection: +(c.sharedPooledWeakDiff - c.pooledWeakDiff).toFixed(1), weighting: +(c.withinWeakDiff - c.sharedPooledWeakDiff).toFixed(1),
    A_shared: wsh(a, shB, true), A_alone: wsh(a, shB, false), B_shared: wsh(b, shA, true), B_alone: wsh(b, shA, false) }; };
out.decomposition = [dec('ABL1', 'AKT1'), dec('AKT1', 'VEGFR2'), dec('ABL1', 'VEGFR2')];

// 5. Informative censoring: does the ceiling track the series' potency?
const gm = new Map(); for (const r of set) { const k = r.group + '|' + r.stratum; if (!gm.has(k)) gm.set(k, []); gm.get(k).push(r); }
const xs = [], ys = [], rows = [];
for (const g of gm.values()) { const b = g.filter(M.isBound), e = g.filter(M.keptByFilter); if (!b.length || e.length < 5) continue;
  const ceil = median(b.map(r => r.v)); const medp = median(e.map(r => pI(r.v))); xs.push(Math.log10(ceil)); ys.push(medp); rows.push({ ceil, medp, frac: b.length / g.length }); }
const rank = a => { const idx = a.map((v, i) => [v, i]).sort((u, v) => u[0] - v[0]); const r = Array(a.length); let i = 0; while (i < idx.length) { let j = i; while (j + 1 < idx.length && idx[j + 1][0] === idx[i][0]) j++; for (let k = i; k <= j; k++) r[idx[k][1]] = (i + j) / 2; i = j + 1; } return r; };
const pear = (a, b) => { const ma = mean(a), mb = mean(b); let s = 0, sa = 0, sb = 0; for (let i = 0; i < a.length; i++) { s += (a[i] - ma) * (b[i] - mb); sa += (a[i] - ma) ** 2; sb += (b[i] - mb) ** 2; } return s / Math.sqrt(sa * sb); };
out.ceilingVsPotency = { sets: xs.length, spearman: +pear(rank(xs), rank(ys)).toFixed(3),
  medianPotencyByCeiling: [['<=10uM', r => r.ceil <= 10000], ['>10uM', r => r.ceil > 10000]].map(([l, f]) => { const z = rows.filter(f); return { ceiling: l, sets: z.length, medianSeriesPIC50: +median(z.map(r => r.medp)).toFixed(2), medianBoundFrac: +median(z.map(r => r.frac)).toFixed(2) }; }) };
// gap between the paper's best exact potency region and its bound: how far below the bound do exact values sit?
const gaps = []; for (const g of gm.values()) { const b = g.filter(M.isBound), e = g.filter(M.keptByFilter); if (!b.length || e.length < 5) continue; gaps.push(pI(median(b.map(r => r.v))) - M.quantile(e.map(r => pI(r.v)).sort((u, v) => u - v), 0.1)); }
out.boundBelowWeakestExact = { sets: gaps.length, medianLogUnitsBelow10thPctExact: +(-median(gaps)).toFixed(2) };

// 6. Sensitivity of the grouped shift to the minimum exact count
out.shiftByMinExact = [3, 5, 10, 20].map(m => { const g = M.groupedEstimates(set, { minExactForShift: m }); return { minExact: m, n: g.shiftGroups, median: +g.shiftMedian.toFixed(2), q1: +g.shiftQ1.toFixed(2), q3: +g.shiftQ3.toFixed(2) }; });
// 7. Bias in a weighted per-target view: pooled censored mean vs exact by target
out.targetBias = [...new Set(set.map(r => r.stratum))].map(s => { const pe = M.pooledEstimates(set.filter(r => r.stratum === s)); const ge = M.groupedEstimates(set.filter(r => r.stratum === s)); return { target: s, pooledBias: +pe.biasLog.toFixed(2), groupedShiftMedian: +ge.shiftMedian.toFixed(2), nShift: ge.shiftGroups, groupShare: +ge.groupShare.toFixed(2) }; });

// 8. Truncation prediction of the per-paper shift: delta = sigma_w * phi(z_c) / (1 - c), z_c = Phi^-1(c)
{
  const sw = G.withinSD, rowsT = [];
  for (const g of gm.values()) { const e = g.filter(M.keptByFilter).length, c = g.filter(M.isBound).length; if (c > 0 && e >= 5) rowsT.push({ c: c / (c + e), nExact: e }); }
  const inv = x => { let lo = -8, hi = 8; for (let i = 0; i < 100; i++) { const m = (lo + hi) / 2; if (M.Phi(m) < x) lo = m; else hi = m; } return (lo + hi) / 2; };
  const phi = z => Math.exp(-z * z / 2) / Math.sqrt(2 * Math.PI);
  rowsT.forEach((r, i) => { r.pred = sw * phi(inv(r.c)) / (1 - r.c); r.obs = G.shifts[i]; });
  const bins = [[5, 9], [10, 19], [20, Infinity]].map(([lo, hi]) => { const z = rowsT.filter(r => r.nExact >= lo && r.nExact <= hi); return { exact: `${lo}-${hi === Infinity ? '' : hi}`, sets: z.length, medianCensoredFrac: +median(z.map(r => r.c)).toFixed(3), observed: +median(z.map(r => r.obs)).toFixed(2), predicted: +median(z.map(r => r.pred)).toFixed(2) }; });
  out.truncation = { sigmaW: +sw.toFixed(3), sets: rowsT.length, medianCensoredFrac: +median(rowsT.map(r => r.c)).toFixed(3), predictedMedian: +median(rowsT.map(r => r.pred)).toFixed(3), observedMedian: +median(rowsT.map(r => r.obs)).toFixed(3), pearson: +pear(rowsT.map(r => r.pred), rowsT.map(r => r.obs)).toFixed(3), bySeriesSize: bins };
  fs.writeFileSync(new URL('./per-paper-shifts.csv', import.meta.url), 'censored_fraction,exact_n,observed_shift,predicted_shift\n' + rowsT.map(r => [r.c.toFixed(4), r.nExact, r.obs.toFixed(4), r.pred.toFixed(4)].join(',')).join('\n') + '\n');
}
console.log(JSON.stringify(out, null, 1));
fs.writeFileSync(new URL('./depth-results.json', import.meta.url), JSON.stringify(out, null, 1));
