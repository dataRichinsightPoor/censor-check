import test from 'node:test'; import assert from 'node:assert/strict'; import fs from 'node:fs';
import * as M from '../web/model.js';
const txt = fs.readFileSync(new URL('../web/data/chembl-six-targets-journal-ic50.csv', import.meta.url), 'utf8');
const p = M.parseDelimited(txt); const set = M.analysisSet(M.toRecords(p, M.mapColumns(p.header)));
test('threshold gradient: censored share of weak evidence at 1, 10 and 100 µM', () => {
  const f = t => M.weakTail(set, t).weakCensoredPct;
  assert.ok(Math.abs(f(1000) - 33.8) < 0.1); assert.ok(Math.abs(f(10000) - 58.8) < 0.1); assert.ok(Math.abs(f(100000) - 82.6) < 0.1);
});
test('ABL1–AKT1 decomposition: selection dominates weighting', () => {
  const c = M.compareStrata(set, 'ABL1', 'AKT1');
  assert.ok(Math.abs((c.sharedPooledWeakDiff - c.pooledWeakDiff) - (-26.5)) < 0.1);
  assert.ok(Math.abs((c.withinWeakDiff - c.sharedPooledWeakDiff) - 2.3) < 0.1);
});
test('truncated-normal formula tracks per-paper shifts', () => {
  const G = M.groupedEstimates(set); const sw = G.withinSD;
  const key = r => r.group + '\u0001' + r.stratum; const gm = new Map();
  for (const r of set) { const k = key(r); if (!gm.has(k)) gm.set(k, []); gm.get(k).push(r); }
  const inv = x => { let lo = -8, hi = 8; for (let i = 0; i < 100; i++) { const m = (lo + hi) / 2; if (M.Phi(m) < x) lo = m; else hi = m; } return (lo + hi) / 2; };
  const phi = z => Math.exp(-z * z / 2) / Math.sqrt(2 * Math.PI);
  const pred = []; for (const g of gm.values()) { const e = g.filter(M.keptByFilter).length, c = g.filter(M.isBound).length; if (c > 0 && e >= 5) { const f = c / (c + e); pred.push(sw * phi(inv(f)) / (1 - f)); } }
  assert.equal(pred.length, G.shifts.length);
  const mean = a => a.reduce((s, x) => s + x, 0) / a.length; const a = pred, b = G.shifts, ma = mean(a), mb = mean(b);
  let q = 0, x = 0, y = 0; for (let i = 0; i < a.length; i++) { q += (a[i] - ma) * (b[i] - mb); x += (a[i] - ma) ** 2; y += (b[i] - mb) ** 2; }
  assert.ok(q / Math.sqrt(x * y) > 0.9);
});
