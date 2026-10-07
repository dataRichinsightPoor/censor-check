// Censor Check: engine. Pure functions, no dependencies.
// Data-Rich, Insight-Poor · 99 Small Problems · No. 06
// MIT License. See LICENSE.

export const VERSION = '0.1.0-alpha';

/* ---------------------------------------------------------------- numerics */

const LOG_HALF = Math.log(0.5);
const LOG_SQRT_2PI = 0.5 * Math.log(2 * Math.PI);

// log(erfc(x)) for x >= 0 with fractional error < 1.2e-7 (Chebyshev fit, Numerical Recipes erfcc),
// evaluated in log space so it never underflows.
function logErfcPos(x) {
  const t = 1 / (1 + 0.5 * x);
  const poly = -x * x - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 +
    t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 +
    t * (-0.82215223 + t * 0.17087277))))))));
  return Math.log(t) + poly;
}

export function erfc(x) {
  return x >= 0 ? Math.exp(logErfcPos(x)) : 2 - Math.exp(logErfcPos(-x));
}

// log of the standard normal CDF, stable far into both tails.
export function logPhi(z) {
  const x = -z / Math.SQRT2;
  if (x >= 0) return LOG_HALF + logErfcPos(x);
  return Math.log1p(-0.5 * Math.exp(logErfcPos(-x)));
}

export function Phi(z) { return Math.exp(logPhi(z)); }
// Inverse normal CDF by bisection (adequate for 0 < q < 1 at the precision used here).
export function invPhi(q) { let lo = -9, hi = 9; for (let i = 0; i < 80; i++) { const m = (lo + hi) / 2; if (Phi(m) < q) lo = m; else hi = m; } return (lo + hi) / 2; }
// Truncated-normal approximation to the exact-only overstatement for a series whose weakest fraction c is censored.
export function truncationShift(c, sigmaW) { if (!(c > 0 && c < 1)) return 0; const z = invPhi(c); return sigmaW * Math.exp(-z * z / 2) / Math.sqrt(2 * Math.PI) / (1 - c); }

export function logNormPdf(x, mu, sigma) {
  const z = (x - mu) / sigma;
  return -LOG_SQRT_2PI - Math.log(sigma) - 0.5 * z * z;
}

// Golden-section maximization of a unimodal function on [a, b].
export function goldenMax(f, a, b, tol = 1e-7, maxIter = 200) {
  const g = (Math.sqrt(5) - 1) / 2;
  let c = b - g * (b - a), d = a + g * (b - a), fc = f(c), fd = f(d);
  for (let i = 0; i < maxIter && Math.abs(b - a) > tol; i++) {
    if (fc > fd) { b = d; d = c; fd = fc; c = b - g * (b - a); fc = f(c); }
    else { a = c; c = d; fc = fd; d = a + g * (b - a); fd = f(d); }
  }
  return (a + b) / 2;
}

// Nelder–Mead maximization in n dimensions.
export function nelderMeadMax(f, x0, step = 0.3, tol = 1e-10, maxIter = 4000) {
  const n = x0.length;
  let simplex = [x0.slice()];
  for (let i = 0; i < n; i++) { const p = x0.slice(); p[i] += step; simplex.push(p); }
  let vals = simplex.map(p => -f(p));
  for (let it = 0; it < maxIter; it++) {
    const order = vals.map((v, i) => i).sort((i, j) => vals[i] - vals[j]);
    simplex = order.map(i => simplex[i]); vals = order.map(i => vals[i]);
    if (Math.abs(vals[n] - vals[0]) < tol) break;
    const cen = new Array(n).fill(0);
    for (let i = 0; i < n; i++) for (let k = 0; k < n; k++) cen[k] += simplex[i][k] / n;
    const pt = (s) => cen.map((c, k) => c + s * (simplex[n][k] - c));
    const xr = pt(-1), fr = -f(xr);
    if (fr < vals[0]) {
      const xe = pt(-2), fe = -f(xe);
      if (fe < fr) { simplex[n] = xe; vals[n] = fe; } else { simplex[n] = xr; vals[n] = fr; }
    } else if (fr < vals[n - 1]) { simplex[n] = xr; vals[n] = fr; }
    else {
      const xc = pt(0.5), fcv = -f(xc);
      if (fcv < vals[n]) { simplex[n] = xc; vals[n] = fcv; }
      else {
        for (let i = 1; i <= n; i++) {
          simplex[i] = simplex[i].map((v, k) => simplex[0][k] + 0.5 * (v - simplex[0][k]));
          vals[i] = -f(simplex[i]);
        }
      }
    }
  }
  return simplex[0];
}

const mean = (a) => a.reduce((s, v) => s + v, 0) / a.length;
const sd = (a) => { const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / a.length); };
export function quantile(a, q) {
  const s = a.slice().sort((x, y) => x - y); if (!s.length) return NaN;
  const h = (s.length - 1) * q, lo = Math.floor(h), hi = Math.ceil(h);
  return s[lo] + (h - lo) * (s[hi] - s[lo]);
}
const median = (a) => quantile(a, 0.5);

/* ------------------------------------------------------------- parsing */

export function detectDelimiter(headerLine) {
  const counts = { ',': 0, ';': 0, '\t': 0 };
  let q = false;
  for (const ch of headerLine) { if (ch === '"') q = !q; else if (!q && ch in counts) counts[ch]++; }
  return Object.entries(counts).sort((a, b) => b[1] - a[1])[0][0];
}

export function parseDelimited(text, delimiter) {
  text = text.replace(/^\uFEFF/, '');
  const firstLine = text.slice(0, text.indexOf('\n') === -1 ? text.length : text.indexOf('\n'));
  const d = delimiter || detectDelimiter(firstLine);
  const rows = []; let row = [], field = '', q = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (q) {
      if (ch === '"') { if (text[i + 1] === '"') { field += '"'; i++; } else q = false; }
      else field += ch;
    } else if (ch === '"') q = true;
    else if (ch === d) { row.push(field); field = ''; }
    else if (ch === '\n' || ch === '\r') {
      if (ch === '\r' && text[i + 1] === '\n') i++;
      row.push(field); field = '';
      if (row.length > 1 || row[0] !== '') rows.push(row);
      row = [];
    } else field += ch;
  }
  if (field !== '' || row.length) { row.push(field); rows.push(row); }
  const header = rows.shift().map(h => h.trim());
  return { header, rows, delimiter: d };
}

// Column roles and the header names that map to them (ChEMBL web export names included).
export const COLUMN_ALIASES = {
  relation: ['relation', 'standard_relation', 'standard relation'],
  value: ['value_nm', 'value', 'standard_value', 'standard value'],
  units: ['units', 'standard_units', 'standard units'],
  group: ['document', 'group', 'paper', 'document_chembl_id', 'document chembl id'],
  compound: ['molecule', 'compound', 'molecule_chembl_id', 'molecule chembl id'],
  stratum: ['target', 'stratum', 'target_pref_name', 'target name', 'target chembl id', 'target_chembl_id'],
  validity: ['validity', 'data_validity_comment', 'data validity comment'],
};

export function mapColumns(header) {
  const low = header.map(h => h.toLowerCase().replace(/^['"]|['"]$/g, '').trim());
  const map = {};
  for (const [role, names] of Object.entries(COLUMN_ALIASES)) {
    const idx = names.map(n => low.indexOf(n)).find(i => i >= 0);
    map[role] = idx === undefined ? -1 : idx;
  }
  return map;
}

const UNIT_TO_NM = { nm: 1, 'um': 1e3, 'µm': 1e3, 'μm': 1e3, 'mm': 1e6, 'pm': 1e-3, 'm': 1e9 };

export function normalizeRelation(r) {
  const s = String(r ?? '').replace(/['"\s]/g, '');
  if (s === '=' ) return 'exact';
  if (s === '>' || s === '>=' || s === '≥' || s === '>>') return 'right';
  if (s === '<' || s === '<=' || s === '≤' || s === '<<') return 'left';
  if (s === '') return 'missing';
  return 'other';
}

// Records: { relation, rel, v (nM), group, compound, stratum, valid, unitOK }
export function toRecords(parsed, map) {
  if (map.relation < 0 || map.value < 0) throw new Error('The file needs a relation column and a value column.');
  const out = [];
  for (const row of parsed.rows) {
    const rawUnits = map.units >= 0 ? String(row[map.units] ?? '').trim().toLowerCase().replace(/['"]/g, '') : 'nm';
    const factor = UNIT_TO_NM[rawUnits];
    const raw = parseFloat(String(row[map.value] ?? '').replace(/['"]/g, ''));
    const validity = map.validity >= 0 ? String(row[map.validity] ?? '').replace(/['"]/g, '').trim() : '';
    out.push({
      rel: normalizeRelation(row[map.relation]),
      v: factor !== undefined && Number.isFinite(raw) ? raw * factor : NaN,
      unitOK: factor !== undefined,
      group: map.group >= 0 ? String(row[map.group] ?? '').trim() || '(none)' : '(all)',
      compound: map.compound >= 0 ? String(row[map.compound] ?? '').trim() : '',
      stratum: map.stratum >= 0 ? String(row[map.stratum] ?? '').trim() || '(none)' : '(all)',
      valid: validity === '' || validity.toLowerCase() === 'manually validated',
    });
  }
  return out;
}

/* ---------------------------------------------------------- classification */

// pChEMBL-style rule: relation '=', units nM, value > 0, validity empty or manually validated.
export function keptByFilter(r) { return r.rel === 'exact' && r.unitOK && r.v > 0 && r.valid; }
export function isBound(r) { return r.rel === 'right' && r.unitOK && r.v > 0; }
export const pIC50 = (vnM) => 9 - Math.log10(vnM);

export function accounting(records) {
  const n = records.length;
  const c = { records: n, exact: 0, right: 0, left: 0, missing: 0, other: 0, badUnits: 0, invalid: 0, kept: 0 };
  for (const r of records) {
    c[r.rel]++;
    if (!r.unitOK || !(r.v > 0)) c.badUnits++;
    if (!r.valid) c.invalid++;
    if (keptByFilter(r)) c.kept++;
  }
  c.removedPct = n ? 100 * (1 - c.kept / n) : NaN;
  return c;
}

// Analysis set: exact records that pass the filter plus right-censored bounds.
export function analysisSet(records) { return records.filter(r => keptByFilter(r) || isBound(r)); }

const key = (r) => r.group + '\u0001' + r.stratum;
function groupBy(arr, f) { const m = new Map(); for (const x of arr) { const k = f(x); if (!m.has(k)) m.set(k, []); m.get(k).push(x); } return m; }

/* ------------------------------------------------------------ weak tail */

export function weakTail(set, T = 10000, minGroup = 5) {
  const isWeak = (r) => r.v >= T;
  const kept = set.filter(keptByFilter);
  const bounds = set.filter(isBound);
  const weakAll = set.filter(isWeak).length, weakKept = kept.filter(isWeak).length;
  const weakBounds = bounds.filter(isWeak).length;
  // paper-weighted
  const groups = groupBy(set, key);
  let pa = [], pk = [], reachable = 0, gCount = 0, rAll = 0, rAllW = 0, rKept = 0, rKeptW = 0, fullyLost = 0, quarterLost = 0, anyBound = 0, gMin = 0;
  for (const g of groups.values()) {
    gCount++;
    const top = Math.max(...g.map(r => r.v));
    if (top >= T) {
      reachable++;
      for (const r of g) { rAll++; if (isWeak(r)) rAllW++; if (keptByFilter(r)) { rKept++; if (isWeak(r)) rKeptW++; } }
    }
    if (g.length >= minGroup) {
      gMin++;
      pa.push(g.filter(isWeak).length / g.length);
      const gk = g.filter(keptByFilter);
      if (gk.length) pk.push(gk.filter(isWeak).length / gk.length);
      const wc = g.filter(r => isBound(r) && isWeak(r)).length, we = gk.filter(isWeak).length;
      if (wc > 0 && we === 0) fullyLost++;
      if (g.some(isBound)) anyBound++;
      if (1 - gk.length / g.length >= 0.25) quarterLost++;
    }
  }
  return {
    T, n: set.length, kept: kept.length, bounds: bounds.length,
    weakAllPct: 100 * weakAll / set.length, weakKeptPct: 100 * weakKept / kept.length,
    weakCensoredPct: weakAll ? 100 * weakBounds / weakAll : NaN,
    paperWeakAllPct: 100 * mean(pa), paperWeakKeptPct: 100 * mean(pk),
    groups: gCount, groupsMin: gMin, reachableGroupsPct: 100 * reachable / gCount,
    reachableWeakAllPct: 100 * rAllW / rAll, reachableWeakKeptPct: 100 * rKeptW / rKept,
    groupsAnyBound: anyBound, groupsQuarterLost: quarterLost, groupsTailLost: fullyLost,
  };
}

// Where bounds sit relative to the threshold.
export function boundPositions(set, T = 10000) {
  const b = set.filter(isBound); const n = b.length || 1;
  return { n: b.length, belowPct: 100 * b.filter(r => r.v < T).length / n, atPct: 100 * b.filter(r => r.v === T).length / n, abovePct: 100 * b.filter(r => r.v > T).length / n };
}

/* -------------------------------------------------- censored-normal models */

// Pooled (single-distribution) censored normal in pIC50 units. A bound '> c nM' means p < pIC50(c).
export function pooledEstimates(set, T = 10000) {
  const e = set.filter(keptByFilter).map(r => pIC50(r.v));
  const c = set.filter(isBound).map(r => pIC50(r.v));
  const ll = ([mu, ls]) => {
    const s = Math.exp(ls); let L = 0;
    for (const x of e) L += logNormPdf(x, mu, s);
    for (const x of c) L += logPhi((x - mu) / s);
    return L;
  };
  const [mu, ls] = nelderMeadMax(ll, [mean(e), Math.log(sd(e) || 1)]);
  const s = Math.exp(ls), pT = pIC50(T);
  return {
    nExact: e.length, nBound: c.length,
    dropMean: mean(e), dropSD: sd(e), substituteMean: mean(e.concat(c)),
    censoredMean: mu, censoredSD: s,
    biasLog: mean(e) - mu, biasFold: Math.pow(10, mean(e) - mu),
    censoredWeakPct: 100 * Phi((pT - mu) / s),
  };
}

// Fixed-effects censored normal: each group has its own mean, all share one within-group SD.
// Coordinate ascent: golden-section on each group mean, then on log SD.
export function groupedEstimates(set, { minExactForShift = 5, iterations = 30 } = {}) {
  const groups = [];
  for (const g of groupBy(set, key).values()) {
    const e = g.filter(keptByFilter).map(r => pIC50(r.v));
    const c = g.filter(isBound).map(r => pIC50(r.v));
    if (e.length >= 1) groups.push({ e, c, n: e.length + c.length, mu: mean(e), dropMean: mean(e) });
  }
  if (!groups.length) return null;
  let s = sd(groups.flatMap(g => g.e)) || 1;
  const gLL = (g, mu, sig) => { let L = 0; for (const x of g.e) L += logNormPdf(x, mu, sig); for (const x of g.c) L += logPhi((x - mu) / sig); return L; };
  for (let it = 0; it < iterations; it++) {
    const sPrev = s;
    for (const g of groups) g.mu = goldenMax(mu => gLL(g, mu, s), g.dropMean - 6, g.dropMean + 3, 1e-6);
    const ls = goldenMax(l => { const sig = Math.exp(l); let L = 0; for (const g of groups) L += gLL(g, g.mu, sig); return L; }, -3, 2, 1e-7);
    s = Math.exp(ls);
    if (Math.abs(s - sPrev) < 1e-6) break;
  }
  const W = groups.reduce((a, g) => a + g.n, 0);
  const muBar = groups.reduce((a, g) => a + g.n * g.mu, 0) / W;
  const between = Math.sqrt(groups.reduce((a, g) => a + g.n * (g.mu - muBar) ** 2, 0) / W);
  const shiftSets = groups.filter(g => g.c.length > 0 && g.e.length >= minExactForShift);
  const shifts = shiftSets.map(g => g.dropMean - g.mu);
  const predicted = shiftSets.map(g => truncationShift(g.c.length / g.n, s));
  const corr = (a, b) => { if (a.length < 3) return NaN; const ma = mean(a), mb = mean(b); let q = 0, x = 0, y = 0; for (let i = 0; i < a.length; i++) { q += (a[i] - ma) * (b[i] - mb); x += (a[i] - ma) ** 2; y += (b[i] - mb) ** 2; } return q / Math.sqrt(x * y); };
  const pooledSD = sd(groups.flatMap(g => g.e));
  return {
    groups: groups.length, withinSD: s, betweenSD: between, pooledExactSD: pooledSD,
    groupShare: between * between / (between * between + s * s),
    shiftGroups: shifts.length, shiftMedian: median(shifts), shiftQ1: quantile(shifts, 0.25), shiftQ3: quantile(shifts, 0.75),
    shiftPositivePct: shifts.length ? 100 * shifts.filter(x => x > 0).length / shifts.length : NaN,
    shifts, predicted,
    censoredFracMedian: shiftSets.length ? median(shiftSets.map(g => g.c.length / g.n)) : NaN,
    predictedMedian: predicted.length ? median(predicted) : NaN, predictedCorr: corr(predicted, shifts),
  };
}

/* ------------------------------------------------ disguised-ceiling test */

// For groups with at least minExact exact values: share with >= k ties at their highest exact value,
// versus the expected share of interior records whose value has >= k ties (exact expectation over a random pick).
export function ceilingTies(set, { minExact = 8, k = 3 } = {}) {
  let n = 0, top = 0, interior = 0, mixed = 0;
  for (const g of groupBy(set, key).values()) {
    const e = g.filter(keptByFilter); if (e.length < minExact) continue;
    n++;
    const counts = new Map(); for (const r of e) counts.set(r.v, (counts.get(r.v) || 0) + 1);
    const mx = Math.max(...e.map(r => r.v));
    if (counts.get(mx) >= k) top++;
    const inner = e.filter(r => r.v !== mx);
    if (inner.length) interior += inner.filter(r => counts.get(r.v) >= k).length / inner.length;
    if (g.some(r => isBound(r) && r.v === mx)) mixed++;
  }
  return { groups: n, topTiedPct: 100 * top / n, interiorTiedPct: 100 * interior / n, mixedCeilings: mixed };
}

/* ------------------------------------------------------- compound level */

export function compoundLevel(set) {
  const withId = set.filter(r => r.compound);
  const pairs = groupBy(withId, r => r.compound + '\u0001' + r.stratum);
  let onlyBound = 0, both = 0, multiGroup = 0;
  for (const g of pairs.values()) {
    const k = g.filter(keptByFilter).length, b = g.filter(isBound).length;
    if (k === 0 && b > 0) onlyBound++;
    if (k > 0 && b > 0) both++;
    if (new Set(g.map(r => r.group)).size > 1) multiGroup++;
  }
  return { pairs: pairs.size, onlyBound, onlyBoundPct: 100 * onlyBound / pairs.size, both, multiGroup, recordsPerPair: withId.length / pairs.size };
}

/* -------------------------------------------- pooled vs within-group comparison */

// Compare stratum A with stratum B on the weak share and on mean pIC50 of exact values.
// pooled: all records; shared: only groups that measured both; within: Mantel–Haenszel-weighted within-group difference.
export function compareStrata(set, A, B, { T = 10000, use = 'all' } = {}) {
  const S = use === 'filtered' ? set.filter(keptByFilter) : set;
  const a = S.filter(r => r.stratum === A), b = S.filter(r => r.stratum === B);
  const weak = (x) => x.filter(r => r.v >= T).length / x.length;
  const ga = groupBy(a, r => r.group), gb = groupBy(b, r => r.group);
  const shared = [...ga.keys()].filter(k => gb.has(k));
  let num = 0, den = 0; const sa = [], sb = [];
  for (const k of shared) {
    const x = ga.get(k), y = gb.get(k), w = x.length * y.length / (x.length + y.length);
    num += w * (weak(x) - weak(y)); den += w; sa.push(...x); sb.push(...y);
  }
  // mean pIC50 among exact values
  const pa = a.filter(keptByFilter), pb = b.filter(keptByFilter);
  const mp = (x) => mean(x.map(r => pIC50(r.v)));
  const gpa = groupBy(pa, r => r.group), gpb = groupBy(pb, r => r.group);
  const sharedP = [...gpa.keys()].filter(k => gpb.has(k));
  const diffs = sharedP.map(k => mp(gpa.get(k)) - mp(gpb.get(k)));
  const pooledP = pa.length && pb.length ? mp(pa) - mp(pb) : NaN;
  return {
    A, B, use, sharedGroups: shared.length, nA: a.length, nB: b.length, nSharedA: sa.length, nSharedB: sb.length,
    pooledWeakDiff: 100 * (weak(a) - weak(b)),
    sharedPooledWeakDiff: shared.length ? 100 * (weak(sa) - weak(sb)) : NaN,
    withinWeakDiff: den ? 100 * num / den : NaN,
    reversed: den ? Math.sign(num) !== Math.sign(weak(a) - weak(b)) && Math.abs(100 * num / den) >= 1 : false,
    pooledPDiff: pooledP, withinPMedian: diffs.length ? median(diffs) : NaN, withinPGroups: diffs.length,
    sameSignPct: diffs.length ? 100 * diffs.filter(d => Math.sign(d) === Math.sign(pooledP)).length / diffs.length : NaN,
  };
}

// In groups that measured two or more strata, call the stratum with the highest median pIC50 (bounds at face value)
// the primary, and the rest counter-screens. Report censoring in each role.
export function counterScreens(set) {
  const byGroup = groupBy(set, r => r.group);
  let pn = 0, pc = 0, cn = 0, cc = 0, groups = 0, singleN = 0, singleC = 0;
  for (const g of byGroup.values()) {
    const strata = groupBy(g, r => r.stratum);
    if (strata.size < 2) { singleN += g.length; singleC += g.filter(isBound).length; continue; }
    groups++;
    let best = null, bestMed = -Infinity;
    for (const [s, x] of strata) { const m = median(x.map(r => pIC50(r.v))); if (m > bestMed) { bestMed = m; best = s; } }
    for (const [s, x] of strata) {
      const c = x.filter(isBound).length;
      if (s === best) { pn += x.length; pc += c; } else { cn += x.length; cc += c; }
    }
  }
  return { groups, primaryCensoredPct: 100 * pc / pn, counterCensoredPct: 100 * cc / cn, singleCensoredPct: 100 * singleC / singleN, primaryN: pn, counterN: cn };
}

/* --------------------------------------------------------------- driver */

export function strataSummary(set, T = 10000) {
  const out = [];
  for (const [s, x] of groupBy(set, r => r.stratum)) {
    const w = weakTail(x, T); out.push({ stratum: s, n: x.length, weakAllPct: w.weakAllPct, weakKeptPct: w.weakKeptPct, weakCensoredPct: w.weakCensoredPct, reachableWeakAllPct: w.reachableWeakAllPct, reachableWeakKeptPct: w.reachableWeakKeptPct });
  }
  return out.sort((a, b) => b.weakAllPct - a.weakAllPct);
}

export function runAll(records, opts = {}) {
  const T = opts.threshold ?? 10000;
  const set = analysisSet(records);
  if (!set.filter(keptByFilter).length) throw new Error('No exact values passed the filter, so there is nothing to compare.');
  return {
    version: VERSION, threshold: T,
    accounting: accounting(records),
    weak: weakTail(set, T, opts.minGroup ?? 5),
    bounds: boundPositions(set, T),
    pooled: pooledEstimates(set, T),
    grouped: groupedEstimates(set, { minExactForShift: opts.minExactForShift ?? 5 }),
    ties: ceilingTies(set),
    compounds: compoundLevel(set),
    strata: strataSummary(set, T),
    counter: counterScreens(set),
  };
}

export function summaryCSV(res) {
  const rows = [['stratum', 'records', 'weak_all_pct', 'weak_after_filter_pct', 'weak_evidence_censored_pct', 'reachable_weak_all_pct', 'reachable_weak_after_filter_pct']];
  for (const s of res.strata) rows.push([s.stratum, s.n, s.weakAllPct, s.weakKeptPct, s.weakCensoredPct, s.reachableWeakAllPct, s.reachableWeakKeptPct].map(v => typeof v === 'number' ? (Number.isInteger(v) ? v : v.toFixed(2)) : v));
  return `# Censor Check ${VERSION}, threshold ${res.threshold} nM\n` + rows.map(r => r.map(c => /[",\n]/.test(String(c)) ? `"${String(c).replace(/"/g, '""')}"` : c).join(',')).join('\n') + '\n';
}
