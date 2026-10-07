import { test } from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import * as M from '../web/model.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg ?? ''} expected ${b} ± ${tol}, got ${a}`);

// Deterministic normal generator (mulberry32 + Box–Muller).
function rng(seed) { let a = seed >>> 0; return () => { a = (a + 0x6D2B79F5) >>> 0; let t = a; t = Math.imul(t ^ (t >>> 15), t | 1); t ^= t + Math.imul(t ^ (t >>> 7), t | 61); return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; }
function normals(n, mu, s, seed) { const u = rng(seed), out = []; for (let i = 0; i < n; i++) { const r = Math.sqrt(-2 * Math.log(u() || 1e-12)), th = 2 * Math.PI * u(); out.push(mu + s * r * Math.cos(th)); } return out; }
const nM = (p) => Math.pow(10, 9 - p);
// Censor at a top concentration: compounds weaker than top become '> top'.
function censoredRecords(ps, topnM, group = 'g', stratum = 's') {
  return ps.map((p, i) => { const v = nM(p); return v >= topnM ? { rel: 'right', v: topnM, unitOK: true, valid: true, group, stratum, compound: `${group}-${i}` } : { rel: 'exact', v, unitOK: true, valid: true, group, stratum, compound: `${group}-${i}` }; });
}

/* numerics */
test('logPhi matches known values and stays finite in both tails', () => {
  close(M.Phi(0), 0.5, 1e-7);
  close(M.Phi(1.959963985), 0.975, 1e-6);
  close(M.Phi(-1.644853627), 0.05, 1e-6);
  close(M.logPhi(-10), -53.23128515051247, 1e-4, 'far left tail');
  close(M.logPhi(-30), -454.3212439955, 1e-3, 'very far left tail');
  close(M.logPhi(8), -6.22096e-16, 1e-15, 'right tail');
  assert.ok(Number.isFinite(M.logPhi(-40)));
});

test('erfc is accurate across its range', () => {
  close(M.erfc(0), 1, 1e-7); close(M.erfc(1), 0.157299207, 1e-7); close(M.erfc(-1), 1.842700793, 1e-7); close(M.erfc(3), 2.20904970e-5, 3e-12);
});

test('goldenMax and nelderMeadMax find known maxima', () => {
  close(M.goldenMax(x => -((x - 1.3) ** 2), -5, 5), 1.3, 1e-6);
  const [a, b] = M.nelderMeadMax(([x, y]) => -((x - 2) ** 2 + 3 * (y + 1) ** 2), [0, 0]);
  close(a, 2, 1e-4); close(b, -1, 1e-4);
});

/* parsing */
test('parses quoted, semicolon-delimited ChEMBL exports and maps columns', () => {
  const txt = '"Molecule ChEMBL ID";"Standard Relation";"Standard Value";"Standard Units";"Document ChEMBL ID";"Target Name";"Data Validity Comment"\n"CHEMBL1";"\'=\'";"12.5";"nM";"CHEMBL9";"Kinase A";""\n"CHEMBL2";"\'>\'";"10";"uM";"CHEMBL9";"Kinase A";"Outside typical range"\n';
  const p = M.parseDelimited(txt); assert.equal(p.delimiter, ';');
  const map = M.mapColumns(p.header); for (const k of ['relation', 'value', 'units', 'group', 'compound', 'stratum', 'validity']) assert.ok(map[k] >= 0, k);
  const r = M.toRecords(p, map);
  assert.equal(r[0].rel, 'exact'); assert.equal(r[1].rel, 'right'); assert.equal(r[1].v, 10000); assert.equal(r[1].valid, false);
  assert.equal(r[0].group, 'CHEMBL9'); assert.equal(r[0].stratum, 'Kinase A');
});

test('relation normalization and unit conversion', () => {
  assert.equal(M.normalizeRelation('>='), 'right'); assert.equal(M.normalizeRelation('≫'.replace('≫', '>>')), 'right');
  assert.equal(M.normalizeRelation('<='), 'left'); assert.equal(M.normalizeRelation('~'), 'other'); assert.equal(M.normalizeRelation(''), 'missing');
  const p = M.parseDelimited('relation,value,units\n=,1,uM\n=,1,pM\n=,1,furlongs\n');
  const r = M.toRecords(p, M.mapColumns(p.header));
  assert.equal(r[0].v, 1000); assert.equal(r[1].v, 0.001); assert.equal(r[2].unitOK, false);
});

test('missing required columns raise a readable error', () => {
  const p = M.parseDelimited('compound,potency\nx,1\n');
  assert.throws(() => M.toRecords(p, M.mapColumns(p.header)), /relation column/);
});

/* filter accounting */
test('filter keeps only valid exact nM values and counts what it removes', () => {
  const recs = [
    { rel: 'exact', v: 10, unitOK: true, valid: true }, { rel: 'exact', v: 10, unitOK: true, valid: false },
    { rel: 'right', v: 1e4, unitOK: true, valid: true }, { rel: 'left', v: 1, unitOK: true, valid: true }, { rel: 'exact', v: 0, unitOK: true, valid: true },
  ].map(r => ({ group: 'g', stratum: 's', compound: '', ...r }));
  const a = M.accounting(recs);
  assert.equal(a.kept, 1); assert.equal(a.right, 1); assert.equal(a.left, 1); assert.equal(a.invalid, 1); close(a.removedPct, 80, 1e-9);
  assert.equal(M.analysisSet(recs).length, 2);
});

/* censored likelihood */
test('pooled censored normal recovers known parameters that exact-only estimates miss', () => {
  const ps = normals(20000, 5.8, 1.2, 7), recs = censoredRecords(ps, 1e4); // top 10 µM = pIC50 5
  const est = M.pooledEstimates(recs, 1e4);
  const m = ps.reduce((a, b) => a + b) / ps.length, s = Math.sqrt(ps.reduce((a, b) => a + (b - m) ** 2, 0) / ps.length);
  close(est.censoredMean, m, 0.01, 'mean recovers the sample mean, including the censored part'); close(est.censoredSD, s, 0.01, 'sd');
  close(est.censoredMean, 5.8, 0.05, 'and the generating mean'); close(est.censoredSD, 1.2, 0.03, 'and the generating sd');
  assert.ok(est.dropMean > 6.1, 'dropping bounds overstates potency');
  assert.ok(est.substituteMean > est.censoredMean && est.substituteMean < est.dropMean, 'substitution sits between');
  close(est.censoredWeakPct, 100 * ps.filter(p => p <= 5).length / ps.length, 0.5, 'predicted weak share matches the sample');
});

test('with no bounds the censored estimate equals the exact-only estimate', () => {
  const recs = normals(500, 7, 1, 3).map(p => ({ rel: 'exact', v: nM(p), unitOK: true, valid: true, group: 'g', stratum: 's', compound: '' }));
  const est = M.pooledEstimates(recs);
  close(est.censoredMean, est.dropMean, 1e-4); close(est.biasLog, 0, 1e-4);
});

test('grouped model separates within-paper and between-paper spread', () => {
  const means = normals(150, 6.5, 1.0, 11); let recs = [];
  means.forEach((m, i) => { recs = recs.concat(censoredRecords(normals(30, m, 0.6, 100 + i), 1e4, 'paper' + i)); });
  const g = M.groupedEstimates(recs);
  close(g.withinSD, 0.6, 0.04, 'within'); close(g.betweenSD, 1.0, 0.15, 'between');
  assert.ok(g.pooledExactSD > g.withinSD, 'pooled SD inflated by between-paper spread');
  assert.ok(g.shifts.every(s => s >= -1e-6), 'bounds can only lower a group mean');
  assert.ok(g.groupShare > 0.6 && g.groupShare < 0.85);
});

/* weak tail and reachability */
test('weak tail counts bounds at or above threshold and reachability excludes papers that could not see weakness', () => {
  const mk = (rel, v, group) => ({ rel, v, unitOK: true, valid: true, group, stratum: 's', compound: '' });
  const recs = [
    mk('exact', 50, 'A'), mk('exact', 80, 'A'), mk('right', 1e4, 'A'), mk('right', 1e4, 'A'), mk('exact', 20000, 'A'),
    mk('exact', 5, 'B'), mk('exact', 9, 'B'), mk('right', 1000, 'B'), mk('exact', 30, 'B'), mk('exact', 40, 'B'),
  ];
  const w = M.weakTail(recs, 1e4, 5);
  close(w.weakAllPct, 30, 1e-9); close(w.weakKeptPct, 100 / 7, 1e-9); close(w.weakCensoredPct, 200 / 3, 1e-9);
  close(w.reachableGroupsPct, 50, 1e-9); close(w.reachableWeakAllPct, 60, 1e-9);
  assert.equal(w.groupsTailLost, 0); assert.equal(w.groupsAnyBound, 2);
  const b = M.boundPositions(recs, 1e4); close(b.belowPct, 100 / 3, 1e-9); close(b.atPct, 200 / 3, 1e-9);
});

/* Simpson-type reversal */
test('pooled comparison reverses within papers when composition differs', () => {
  const mk = (stratum, group, weak) => ({ rel: weak ? 'right' : 'exact', v: weak ? 1e4 : 100, unitOK: true, valid: true, group, stratum, compound: '' });
  const recs = [];
  // Paper P1: mostly A, A is weaker than B (60% vs 80%). Paper P2: mostly B, both potent; A weaker (10% vs 20%)?
  // Construct: within each paper, A has FEWER weak records than B, but A is concentrated in the paper where everything is weak.
  for (let i = 0; i < 90; i++) recs.push(mk('A', 'P1', i < 72)); // A in P1: 80% weak
  for (let i = 0; i < 10; i++) recs.push(mk('B', 'P1', i < 9));  // B in P1: 90% weak
  for (let i = 0; i < 10; i++) recs.push(mk('A', 'P2', i < 1));  // A in P2: 10% weak
  for (let i = 0; i < 90; i++) recs.push(mk('B', 'P2', i < 18)); // B in P2: 20% weak
  const c = M.compareStrata(recs, 'A', 'B');
  assert.ok(c.pooledWeakDiff > 40, 'pooled: A looks much weaker');
  assert.ok(c.withinWeakDiff < 0, 'within papers: A is less often weak');
  close(c.withinWeakDiff, -10, 1e-9); assert.equal(c.reversed, true); assert.equal(c.sharedGroups, 2);
});

test('within-paper difference equals the plain difference for a single shared paper', () => {
  const mk = (stratum, weak) => ({ rel: 'exact', v: weak ? 2e4 : 10, unitOK: true, valid: true, group: 'P', stratum, compound: '' });
  const recs = [mk('A', 1), mk('A', 0), mk('A', 0), mk('A', 0), mk('B', 1), mk('B', 1)];
  const c = M.compareStrata(recs, 'A', 'B'); close(c.withinWeakDiff, c.pooledWeakDiff, 1e-9); close(c.withinWeakDiff, -75, 1e-9);
});

test('counter-screen check assigns primary by median potency', () => {
  const mk = (stratum, rel, v) => ({ rel, v, unitOK: true, valid: true, group: 'P', stratum, compound: '' });
  const recs = [mk('On', 'exact', 5), mk('On', 'exact', 8), mk('Off', 'right', 1e4), mk('Off', 'exact', 3000), { ...mk('X', 'exact', 50), group: 'Q' }];
  const c = M.counterScreens(recs);
  assert.equal(c.groups, 1); close(c.primaryCensoredPct, 0, 1e-9); close(c.counterCensoredPct, 50, 1e-9); close(c.singleCensoredPct, 0, 1e-9);
});

test('ceiling-tie check flags disguised ceilings and passes clean data', () => {
  const mk = (rel, v, g) => ({ rel, v, unitOK: true, valid: true, group: g, stratum: 's', compound: '' });
  const disguised = [5, 9, 14, 30, 60, 120, 1e4, 1e4, 1e4, 1e4].map(v => mk('exact', v, 'D'));
  const clean = [5, 9, 14, 30, 60, 120, 400, 900].map(v => mk('exact', v, 'C')).concat([mk('right', 1e4, 'C')]);
  const t = M.ceilingTies(disguised.concat(clean));
  assert.equal(t.groups, 2); close(t.topTiedPct, 50, 1e-9); close(t.interiorTiedPct, 0, 1e-9); assert.equal(t.mixedCeilings, 0);
});

test('compound level counts pairs known only through bounds', () => {
  const mk = (compound, rel, group) => ({ rel, v: rel === 'right' ? 1e4 : 10, unitOK: true, valid: true, group, stratum: 's', compound });
  const c = M.compoundLevel([mk('a', 'exact', 'P'), mk('a', 'right', 'Q'), mk('b', 'right', 'P'), mk('c', 'exact', 'P')]);
  assert.equal(c.pairs, 3); assert.equal(c.onlyBound, 1); assert.equal(c.both, 1); assert.equal(c.multiGroup, 1);
});

test('runAll rejects a table with no exact values', () => {
  assert.throws(() => M.runAll([{ rel: 'right', v: 1e4, unitOK: true, valid: true, group: 'g', stratum: 's', compound: '' }]), /No exact values/);
});

test('summary CSV is labeled with version and threshold', () => {
  const recs = censoredRecords(normals(200, 6, 1, 5), 1e4);
  const csv = M.summaryCSV(M.runAll(recs));
  assert.match(csv, new RegExp(`Censor Check ${M.VERSION.replace(/\./g, '\\.')}, threshold 10000 nM`));
  assert.match(csv, /^stratum,records,weak_all_pct/m);
});

/* regression: the numbers published in The Audit, No. 1 */
test('bundled ChEMBL subset reproduces the published numbers', () => {
  const p = M.parseDelimited(fs.readFileSync(new URL('../web/data/chembl-six-targets-journal-ic50.csv', import.meta.url), 'utf8'));
  const recs = M.toRecords(p, M.mapColumns(p.header));
  const r = M.runAll(recs);
  assert.equal(r.accounting.records, 31236); assert.equal(r.accounting.kept, 26575); assert.equal(r.accounting.right, 3971);
  close(r.accounting.removedPct, 14.92, 0.01);
  close(r.weak.weakAllPct, 17.77, 0.01); close(r.weak.weakKeptPct, 8.41, 0.01); close(r.weak.weakCensoredPct, 58.80, 0.01);
  close(r.weak.paperWeakAllPct, 19.00, 0.01); close(r.weak.paperWeakKeptPct, 11.14, 0.01);
  close(r.weak.reachableGroupsPct, 44.42, 0.01); close(r.weak.reachableWeakAllPct, 31.53, 0.01); close(r.weak.reachableWeakKeptPct, 16.17, 0.01);
  assert.equal(r.weak.groupsMin, 1347); assert.equal(r.weak.groupsAnyBound, 594); assert.equal(r.weak.groupsTailLost, 235);
  close(r.bounds.belowPct, 19.64, 0.01);
  close(r.pooled.dropMean, 6.780, 0.001); close(r.pooled.substituteMean, 6.524, 0.001); close(r.pooled.censoredMean, 6.420, 0.002); close(r.pooled.censoredWeakPct, 17.53, 0.05);
  close(r.grouped.withinSD, 0.910, 0.003); close(r.grouped.betweenSD, 1.130, 0.003); close(r.grouped.shiftMedian, 0.304, 0.003); assert.equal(r.grouped.shiftGroups, 486);
  close(r.ties.topTiedPct, 0.71, 0.01); close(r.ties.interiorTiedPct, 3.99, 0.01); assert.equal(r.ties.groups, 985); assert.equal(r.ties.mixedCeilings, 19);
  assert.equal(r.compounds.pairs, 23780); assert.equal(r.compounds.onlyBound, 3525);
  close(r.counter.counterCensoredPct, 29.62, 0.01); close(r.counter.primaryCensoredPct, 13.35, 0.01); assert.equal(r.counter.groups, 197);
  const set = M.analysisSet(recs);
  const ab = M.compareStrata(set, 'ABL1', 'AKT1'); close(ab.pooledWeakDiff, 4.46, 0.01); close(ab.withinWeakDiff, -19.70, 0.01); assert.equal(ab.sharedGroups, 50); assert.equal(ab.reversed, true);
  const av = M.compareStrata(set, 'AKT1', 'VEGFR2'); close(av.pooledWeakDiff, 3.61, 0.01); close(av.withinWeakDiff, 30.13, 0.01);
  const v = r.strata.find(s => s.stratum === 'VEGFR2'); close(v.weakAllPct, 15.17, 0.01); close(v.weakKeptPct, 4.92, 0.01); close(v.weakCensoredPct, 72.46, 0.01);
});

test('version is consistent across package and pages', () => {
  const pkg = JSON.parse(fs.readFileSync(new URL('../package.json', import.meta.url), 'utf8'));
  assert.equal(pkg.version, M.VERSION);
  for (const f of ['../web/index.html', '../web/methods.html']) assert.ok(fs.readFileSync(new URL(f, import.meta.url), 'utf8').includes(M.VERSION), f);
});
