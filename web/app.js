import * as M from './model.js';

const $ = (id) => document.getElementById(id);
let records = null, result = null, set = null, sourceName = '';

const css = (v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
const f1 = (x) => Number.isFinite(x) ? x.toFixed(1) : 'n/a';
const f2 = (x) => Number.isFinite(x) ? x.toFixed(2) : 'n/a';
const int = (x) => Number.isFinite(x) ? Math.round(x).toLocaleString('en-US') : 'n/a';

function status(msg, err = false) { $('status').textContent = msg; $('status').className = err ? 'error' : ''; }

/* --------------------------------------------------------------- loading */

async function load() {
  try {
    if ($('source').value === 'bundled') {
      status('Loading bundled ChEMBL subset…');
      const txt = await (await fetch('data/chembl-six-targets-journal-ic50.csv')).text();
      ingest(txt, 'Bundled ChEMBL subset');
    } else {
      const f = $('file').files[0];
      if (!f) { status('Choose a file first.', true); return; }
      ingest(await f.text(), f.name);
    }
  } catch (e) { status(e.message, true); }
}

function ingest(text, name) {
  const parsed = M.parseDelimited(text);
  const map = M.mapColumns(parsed.header);
  records = M.toRecords(parsed, map);
  sourceName = name;
  const found = Object.entries(map).filter(([, i]) => i >= 0).map(([r, i]) => `${r} → “${parsed.header[i]}”`);
  $('columns').textContent = `${int(records.length)} rows. Matched: ${found.join('; ')}.`;
  run();
}

/* ------------------------------------------------------------------ run */

function run(ev) {
  if (ev) ev.preventDefault();
  if (!records) { load(); return; }
  const T = parseFloat($('threshold').value), mg = parseInt($('minGroup').value, 10), ms = parseInt($('minShift').value, 10);
  if (!(T > 0)) { status('The threshold must be a positive concentration in nM.', true); return; }
  if (!(mg >= 1) || !(ms >= 1)) { status('Group minimums must be whole numbers of at least 1.', true); return; }
  try {
    status('Running…');
    result = M.runAll(records, { threshold: T, minGroup: mg, minExactForShift: ms });
    set = M.analysisSet(records);
    render();
    document.body.classList.remove('stale');
    ['csv', 'json', 'png'].forEach(id => $(id).disabled = false);
    status(`Done · ${sourceName} · threshold ${int(T)} nM.`);
  } catch (e) { status(e.message, true); }
}

/* --------------------------------------------------------------- render */

function render() {
  const r = result, w = r.weak, p = r.pooled, g = r.grouped, a = r.accounting;
  $('headline').textContent = `${f1(w.weakCensoredPct)}% of the weak evidence is a bound. The exact-value filter removes it.`;
  $('badge').textContent = w.weakCensoredPct >= 33 ? 'TAIL MOSTLY CENSORED' : w.weakCensoredPct >= 10 ? 'TAIL PARTLY CENSORED' : 'TAIL MOSTLY EXACT';
  $('metrics').innerHTML = [
    ['Weak, all evidence', f1(w.weakAllPct) + '%', `records at or above ${int(r.threshold)} nM`],
    ['Weak, after filter', f1(w.weakKeptPct) + '%', 'exact values only'],
    ['Papers losing their whole weak tail', int(w.groupsTailLost), `of ${int(w.groupsMin)} with ≥ ${$('minGroup').value} records`],
    ['Median per-paper overstatement', g ? f2(g.shiftMedian) : 'n/a', g ? `log units · ${f1(Math.pow(10, g.shiftMedian))}-fold` : ''],
  ].map(([s, v, d]) => `<div class="metric"><small>${s}</small><strong>${v}</strong><p>${d}</p></div>`).join('');

  $('accounting').innerHTML = [
    ['Records read', int(a.records), 'Every row with a relation.'],
    ['Exact ("=")', int(a.exact), 'Claims a located midpoint.'],
    ['Right-censored (">", "≥", "≫")', int(a.right), 'Curve not crossed in the tested range: potency is below the bound.'],
    ['Left-censored ("<", "≤")', int(a.left), 'More potent than the lowest concentration. Not modeled here.'],
    ['Validity-flagged', int(a.invalid), 'Comment other than empty or manually validated.'],
    ['Kept by the exact-value filter', int(a.kept), `${f1(a.removedPct)}% of rows removed.`],
    ['Compound–stratum pairs known only through bounds', int(r.compounds.onlyBound), `${f1(r.compounds.onlyBoundPct)}% of ${int(r.compounds.pairs)} pairs vanish entirely.`],
  ].map(([q, v, m]) => `<tr><td>${q}</td><td class="num">${v}</td><td>${m}</td></tr>`).join('');

  drawBars();
  $('agg-cards').innerHTML = g ? [
    [g.groupShare > 0.5 ? 'flag' : 'ok', 'Paper explains the spread', `Between-paper SD ${f2(g.betweenSD)} versus within-paper SD ${f2(g.withinSD)}: ${f1(100 * g.groupShare)}% of potency variance sits between papers. A single pooled distribution is misspecified.`],
    ['flag', 'Weighting changes the answer', `Record-weighted weak share ${f1(w.weakAllPct)}% → ${f1(w.weakKeptPct)}%. Paper-weighted ${f1(w.paperWeakAllPct)}% → ${f1(w.paperWeakKeptPct)}%. Large papers set the pooled number.`],
    [w.reachableGroupsPct < 60 ? 'flag' : 'ok', 'Not every paper could see weakness', `Only ${f1(w.reachableGroupsPct)}% of paper sets report any value or bound at or above the threshold. Among those, weak share is ${f1(w.reachableWeakAllPct)}% → ${f1(w.reachableWeakKeptPct)}%.`],
  ].map(([c, h, t]) => `<article class="${c}"><h3>${h}</h3><p>${t}</p></article>`).join('') : '';
  drawShifts();
  $('shift-note').textContent = g ? `${int(g.shiftGroups)} paper sets with at least one bound and ≥ ${$('minShift').value} exact values. Median ${f2(g.shiftMedian)}, interquartile range ${f2(g.shiftQ1)} to ${f2(g.shiftQ3)}. Every shift is positive by construction: a bound can only pull a group mean down. The magnitude, not the sign, is the finding. The rightmost bar includes all larger shifts. From censored fractions alone (median ${f2(g.censoredFracMedian)}), the truncated-normal formula predicts a median of ${f2(g.predictedMedian)}; prediction and fit correlate at r = ${f2(g.predictedCorr)}.` : '';
  $('estimates').innerHTML = [
    ['Exact values only, mean', f2(p.dropMean), 'What a filtered dataset reports.'],
    ['Bounds substituted as values, mean', f2(p.substituteMean), 'Treats "> 10 µM" as "= 10 µM". Still biased toward potency.'],
    ['Pooled censored normal, mean', f2(p.censoredMean), `SD ${f2(p.censoredSD)}. Overstatement ${f2(p.biasLog)} log units (${f1(p.biasFold)}-fold).`],
    ['Pooled censored normal, predicted weak share', f1(p.censoredWeakPct) + '%', `Compare ${f1(w.weakAllPct)}% counted directly.`],
    ['Grouped model, within-paper SD', g ? f2(g.withinSD) : 'n/a', 'Shared by all papers.'],
    ['Grouped model, between-paper SD', g ? f2(g.betweenSD) : 'n/a', 'Record-weighted SD of paper means.'],
  ].map(([q, v, m]) => `<tr><td>${q}</td><td class="num">${v}</td><td>${m}</td></tr>`).join('');
  drawBoundsBar();

  const strata = r.strata.map(s => s.stratum);
  for (const id of ['cmpA', 'cmpB']) $(id).innerHTML = strata.map(s => `<option>${s}</option>`).join('');
  if (strata.length > 1) { $('cmpA').value = strata.includes('ABL1') ? 'ABL1' : strata[0]; $('cmpB').value = strata.includes('AKT1') ? 'AKT1' : strata[1]; }
  compare();
  const c = r.counter;
  $('counter-cards').innerHTML = c.groups ? [
    ['flag', 'Counter-screens are censored more', `In ${int(c.groups)} papers measuring two or more strata, ${f1(c.counterCensoredPct)}% of counter-screen records are bounds, against ${f1(c.primaryCensoredPct)}% for the stratum each series was optimized for. The filter removes selectivity evidence first.`],
    ['', 'What a within-paper contrast measures', 'Inside a selectivity panel, the off-target is weak by design. The within-paper difference measures the series' + "'" + ' intended selectivity, not the intrinsic difficulty of either target.'],
    ['', 'Single-stratum papers', `${f1(c.singleCensoredPct)}% of their records are bounds.`],
  ].map(([k, h, t]) => `<article class="${k}"><h3>${h}</h3><p>${t}</p></article>`).join('') : '<article><h3>One stratum per paper</h3><p>No paper measured more than one stratum, so no within-paper contrast or counter-screen check is possible.</p></article>';
  const t = r.ties;
  $('tie-cards').innerHTML = [
    [t.topTiedPct > t.interiorTiedPct ? 'flag' : 'ok', 'Ties at the top', `${f1(t.topTiedPct)}% of ${int(t.groups)} paper sets with ≥ 8 exact values have three or more compounds tied at their highest value.`],
    ['', 'Ties in the interior', `Expected for a random interior record: ${f1(t.interiorTiedPct)}%.`],
    [t.mixedCeilings > 0 ? 'flag' : 'ok', 'Same ceiling written both ways', `${int(t.mixedCeilings)} paper sets report their highest exact value also as a bound for other compounds.`],
  ].map(([k, h, x]) => `<article class="${k}"><h3>${h}</h3><p>${x}</p></article>`).join('');
}

function compare() {
  if (!set) return;
  const A = $('cmpA').value, B = $('cmpB').value;
  if (!A || !B || A === B) { $('cmp-note').textContent = 'Choose two different strata.'; clear($('cmp')); return; }
  const c = M.compareStrata(set, A, B, { T: result.threshold, use: $('cmpUse').value });
  drawCompare(c);
  $('cmp-note').textContent = c.sharedGroups
    ? `${A} minus ${B}, weak share in percentage points. ${int(c.sharedGroups)} papers measured both (${int(c.nSharedA)} and ${int(c.nSharedB)} records). ${c.reversed ? 'The within-paper difference has the opposite sign. A pooled sign can be fragile; resample papers (analysis/depth.mjs) before reading it as a reversal.' : 'Same sign pooled and within papers; the magnitude may still differ.'} Mean pIC50 of exact values: pooled ${f2(c.pooledPDiff)}, within-paper median ${f2(c.withinPMedian)} across ${int(c.withinPGroups)} papers.`
    : `No paper measured both ${A} and ${B}. The pooled difference (${f1(c.pooledWeakDiff)} points) cannot be checked within papers.`;
}

/* -------------------------------------------------------------- drawing */

function setup(cv) {
  const dpr = window.devicePixelRatio || 1, W = cv.clientWidth, H = cv.clientHeight;
  cv.width = W * dpr; cv.height = H * dpr; const ctx = cv.getContext('2d'); ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.fillStyle = css('--surface'); ctx.fillRect(0, 0, W, H); ctx.font = '12px DM Sans, sans-serif'; return { ctx, W, H };
}
function clear(cv) { setup(cv); }

function drawBars() {
  const { ctx, W, H } = setup($('bars')), reach = $('reach').checked;
  const S = result.strata, L = 90, R = 50, top = 10, bottom = 34, rowH = (H - top - bottom) / S.length;
  const vals = S.flatMap(s => reach ? [s.reachableWeakAllPct, s.reachableWeakKeptPct] : [s.weakAllPct, s.weakKeptPct]).filter(Number.isFinite);
  const xmax = Math.max(10, Math.ceil(Math.max(...vals) / 10) * 10), x = (v) => L + (W - L - R) * v / xmax;
  ctx.strokeStyle = css('--line'); ctx.fillStyle = css('--muted');
  for (let v = 0; v <= xmax; v += 10) { ctx.beginPath(); ctx.moveTo(x(v), top); ctx.lineTo(x(v), H - bottom); ctx.stroke(); ctx.fillText(v + '%', x(v) - 8, H - bottom + 16); }
  S.forEach((s, i) => {
    const y0 = top + i * rowH, bh = Math.min(16, rowH * 0.34);
    const [va, vk] = reach ? [s.reachableWeakAllPct, s.reachableWeakKeptPct] : [s.weakAllPct, s.weakKeptPct];
    ctx.fillStyle = css('--ink'); ctx.fillText(s.stratum.slice(0, 12), 6, y0 + rowH / 2 + 4);
    ctx.fillStyle = css('--teal'); ctx.fillRect(L, y0 + rowH / 2 - bh - 1, x(va) - L, bh);
    ctx.fillStyle = css('--amber'); ctx.fillRect(L, y0 + rowH / 2 + 1, Math.max(0, x(vk) - L), bh);
    ctx.fillStyle = css('--muted'); ctx.fillText(f1(va) + '%', x(va) + 5, y0 + rowH / 2 - 4); ctx.fillText(f1(vk) + '%', x(vk) + 5, y0 + rowH / 2 + bh - 2);
  });
  $('bars-badge').textContent = reach ? 'PAPERS THAT REACH THE THRESHOLD' : 'ALL PAPERS';
}

function drawShifts() {
  const cv = $('shifts'), { ctx, W, H } = setup(cv), g = result.grouped; if (!g || !g.shifts.length) return;
  const lo = 0, hi = Math.max(1.5, Math.ceil(quant(g.shifts, 0.98) * 4) / 4), nb = 30, bw = (hi - lo) / nb;
  const bins = new Array(nb).fill(0); for (const s of g.shifts) bins[Math.min(nb - 1, Math.max(0, Math.floor((s - lo) / bw)))]++;
  const L = 40, R = 16, top = 10, bottom = 34, ymax = Math.max(...bins), x = (v) => L + (W - L - R) * (v - lo) / (hi - lo), y = (c) => H - bottom - (H - top - bottom) * c / ymax;
  ctx.strokeStyle = css('--line'); ctx.beginPath(); ctx.moveTo(L, H - bottom); ctx.lineTo(W - R, H - bottom); ctx.stroke();
  ctx.fillStyle = css('--teal'); bins.forEach((c, i) => ctx.fillRect(x(lo + i * bw) + 1, y(c), (W - L - R) / nb - 2, H - bottom - y(c)));
  ctx.fillStyle = css('--muted'); for (let v = 0; v <= hi + 1e-9; v += 0.25) ctx.fillText(v.toFixed(2), x(v) - 10, H - bottom + 16);
  ctx.fillText('papers', 2, top + 10);
  ctx.strokeStyle = css('--amber'); ctx.setLineDash([4, 3]); ctx.beginPath(); ctx.moveTo(x(g.shiftMedian), top); ctx.lineTo(x(g.shiftMedian), H - bottom); ctx.stroke(); ctx.setLineDash([]);
  ctx.fillStyle = css('--amber'); ctx.fillText('median ' + f2(g.shiftMedian), x(g.shiftMedian) + 6, top + 12);
}
function quant(a, q) { return M.quantile(a, q); }

function drawBoundsBar() {
  const { ctx, W, H } = setup($('boundsbar')), b = result.bounds, L = 10, R = 10, y = 22, h = 30;
  let x0 = L; const tot = W - L - R;
  [[b.belowPct, '--floor', 'below threshold'], [b.atPct, '--amber', 'at threshold'], [b.abovePct, '--teal', 'above threshold']].forEach(([v, c, lab]) => {
    const w = tot * v / 100; ctx.fillStyle = css(c); ctx.fillRect(x0, y, w, h);
    ctx.fillStyle = css('--ink'); if (w > 90) ctx.fillText(`${f1(v)}% ${lab}`, x0 + 6, y + h + 18); x0 += w;
  });
}

function drawCompare(c) {
  const { ctx, W, H } = setup($('cmp'));
  const items = [['Pooled, all papers', c.pooledWeakDiff], ['Pooled, shared papers only', c.sharedPooledWeakDiff], ['Within-paper (weighted)', c.withinWeakDiff]];
  const m = Math.max(5, ...items.map(([, v]) => Math.abs(v || 0))) * 1.15, L = 190, R = 50, mid = L + (W - L - R) / 2, sc = (W - L - R) / 2 / m, rowH = (H - 30) / 3;
  ctx.strokeStyle = css('--line'); ctx.beginPath(); ctx.moveTo(mid, 6); ctx.lineTo(mid, H - 24); ctx.stroke();
  ctx.fillStyle = css('--muted'); ctx.fillText(`${c.B} weaker ←`, L, H - 6); const t = `→ ${c.A} weaker`; ctx.fillText(t, W - R - ctx.measureText(t).width, H - 6);
  items.forEach(([lab, v], i) => {
    const yy = 8 + i * rowH; ctx.fillStyle = css('--ink'); ctx.fillText(lab, 6, yy + rowH / 2);
    if (!Number.isFinite(v)) { ctx.fillStyle = css('--muted'); ctx.fillText('no shared papers', mid + 6, yy + rowH / 2); return; }
    ctx.fillStyle = i === 2 ? css('--amber') : css('--teal'); const w = v * sc; ctx.fillRect(Math.min(mid, mid + w), yy + rowH / 2 - 10, Math.abs(w), 16);
    ctx.fillStyle = css('--ink'); const s = (v > 0 ? '+' : '') + f1(v); ctx.fillText(s, v >= 0 ? mid + w + 6 : mid + w - ctx.measureText(s).width - 6, yy + rowH / 2 + 2);
  });
}

/* -------------------------------------------------------------- exports */

function download(name, blob) { const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 1000); $('export-status').textContent = `Saved ${name}.`; }

function init() {
  $('controls').addEventListener('submit', run);
  $('load').addEventListener('click', load);
  $('source').addEventListener('change', () => { $('file-label').hidden = $('source').value !== 'file'; });
  $('file').addEventListener('change', load);
  $('reach').addEventListener('change', () => result && drawBars());
  $('cmpRun').addEventListener('click', compare);
  ['cmpA', 'cmpB', 'cmpUse'].forEach(id => $(id).addEventListener('change', compare));
  ['threshold', 'minGroup', 'minShift'].forEach(id => $(id).addEventListener('input', () => { if (result) { document.body.classList.add('stale'); ['csv', 'json', 'png'].forEach(b => $(b).disabled = true); status('Inputs changed. Run the check again.'); } }));
  $('csv').addEventListener('click', () => download('censor-check-summary.csv', new Blob([M.summaryCSV(result)], { type: 'text/csv' })));
  $('json').addEventListener('click', () => { const r = structuredClone(result); if (r.grouped) delete r.grouped.shifts; download('censor-check-results.json', new Blob([JSON.stringify({ source: sourceName, ...r }, null, 2)], { type: 'application/json' })); });
  $('png').addEventListener('click', () => $('bars').toBlob(b => download('censor-check-weak-tail.png', b)));
  $('theme').addEventListener('click', () => { const d = document.documentElement, light = d.dataset.theme === 'dark'; d.dataset.theme = light ? 'light' : 'dark'; $('theme').textContent = light ? 'Dark mode' : 'Light mode'; if (result) render(); });
  window.addEventListener('resize', () => { if (result) { drawBars(); drawShifts(); drawBoundsBar(); compare(); } });
  load();
}
init();
