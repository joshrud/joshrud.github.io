// Nucleotide encoding: A=0, C=1, G=2, T=3

function rnd(n) { return Math.floor(Math.random() * n); }

function isStop(a, b, c) {
  return (a === 3 && b === 0 && (c === 0 || c === 2)) ||
         (a === 3 && b === 2 && c === 0);
}

// ── Generate ──────────────────────────────────────────────────────────────────

function generate(len) {
  const seq = new Uint8Array(len);

  // Background: human-like frequencies (~29% A, 21% C, 21% G, 29% T)
  for (let i = 0; i < len; i++) {
    const r = Math.random() * 100;
    seq[i] = r < 29 ? 0 : r < 50 ? 1 : r < 71 ? 2 : 3;
  }

  // Plant CpG islands: high GC (≈65%), CpG dinucleotide-enriched
  const numCpG = 12 + rnd(5);
  for (let n = 0; n < numCpG; n++) {
    const s = rnd(len - 4000) + 500;
    const isLen = 500 + rnd(1500);
    for (let i = s; i < s + isLen && i < len; i++) {
      const r = Math.random();
      seq[i] = r < 0.30 ? 1 : r < 0.60 ? 2 : r < 0.80 ? 0 : 3;
    }
    // Seed extra CpG dinucleotides within the island
    for (let i = s; i < s + isLen - 1 && i < len - 1; i++) {
      if (Math.random() < 0.07) { seq[i] = 1; seq[i + 1] = 2; i++; }
    }
  }

  // Plant gene-like structures: TATA box + ATG + non-stop ORF + stop codon
  const numGenes = 20 + rnd(11);
  for (let n = 0; n < numGenes; n++) {
    const s = rnd(len - 12000) + 1000;
    // TATA box: T-A-T-A-A-A
    [3, 0, 3, 0, 0, 0].forEach((v, j) => { if (s + j < len) seq[s + j] = v; });
    // ATG start codon ~30 bp downstream
    const atg = s + 30;
    if (atg + 2 < len) { seq[atg] = 0; seq[atg + 1] = 3; seq[atg + 2] = 2; }
    // Coding sequence with no in-frame stop codons
    let pos = atg + 3;
    const orfLen = 600 + rnd(1400);
    for (let i = 0; i < orfLen && pos + 2 < len; i += 3) {
      let a, b, c;
      do { a = rnd(4); b = rnd(4); c = rnd(4); } while (isStop(a, b, c));
      seq[pos] = a; seq[pos + 1] = b; seq[pos + 2] = c;
      pos += 3;
    }
    // TAA stop codon
    if (pos + 2 < len) { seq[pos] = 3; seq[pos + 1] = 0; seq[pos + 2] = 0; }
  }

  return seq;
}

// ── Parse ─────────────────────────────────────────────────────────────────────

function parse(seq) {
  const len = seq.length;

  // GC% in 500 equal bins
  const BINS = 500;
  const binSize = len / BINS;
  const gcBins = new Float32Array(BINS);
  for (let b = 0; b < BINS; b++) {
    const s = Math.floor(b * binSize), e = Math.floor((b + 1) * binSize);
    let gc = 0;
    for (let i = s; i < e; i++) if (seq[i] === 1 || seq[i] === 2) gc++;
    gcBins[b] = gc / (e - s);
  }

  // CpG islands: window=200bp, step=50bp; GC≥50%, CpG obs/exp≥0.60
  // (Gardiner-Garden & Frommer 1987 criteria)
  const cpgIslands = [];
  let inIsland = false, islandStart = 0;
  for (let i = 0; i + 200 <= len; i += 50) {
    let c = 0, g = 0, cpg = 0;
    for (let j = i; j < i + 200; j++) {
      if (seq[j] === 1) c++;
      else if (seq[j] === 2) g++;
      if (j + 1 < len && seq[j] === 1 && seq[j + 1] === 2) cpg++;
    }
    const gcPct = (c + g) / 200;
    const cpgOE = (cpg * 200) / ((c || 1) * (g || 1));
    const hit = gcPct >= 0.50 && cpgOE >= 0.60;
    if (hit && !inIsland) { inIsland = true; islandStart = i; }
    if (!hit && inIsland) {
      if (i - islandStart >= 200) cpgIslands.push({ start: islandStart, end: i });
      inIsland = false;
    }
  }
  if (inIsland && len - islandStart >= 200) cpgIslands.push({ start: islandStart, end: len });

  // ORFs: greedy forward-strand scan, ATG → first in-frame stop, min 150 bp
  const orfs = [];
  let i = 0;
  while (i < len - 150) {
    if (seq[i] === 0 && seq[i + 1] === 3 && seq[i + 2] === 2) { // ATG
      let j = i + 3;
      while (j + 2 < len && !isStop(seq[j], seq[j + 1], seq[j + 2])) j += 3;
      const orfLen = j - i;
      if (orfLen >= 150) orfs.push({ start: i, end: j + 3, long: orfLen >= 500 });
      i = j + 3;
    } else {
      i++;
    }
  }

  // TATA boxes: TATAAA motif (putative TSS ~30 bp downstream)
  const tataBoxes = [];
  for (let i = 0; i < len - 6; i++) {
    if (seq[i] === 3 && seq[i+1] === 0 && seq[i+2] === 3 &&
        seq[i+3] === 0 && seq[i+4] === 0 && seq[i+5] === 0) {
      tataBoxes.push(i);
      i += 5;
    }
  }

  let totalGC = 0;
  for (let i = 0; i < len; i++) if (seq[i] === 1 || seq[i] === 2) totalGC++;

  return { gcBins, cpgIslands, orfs, tataBoxes, gcOverall: totalGC / len };
}

// ── Render DNA ────────────────────────────────────────────────────────────────

function render(canvas, W, H, len, { gcBins, cpgIslands, orfs, tataBoxes }) {
  const ctx = canvas.getContext('2d');
  ctx.fillStyle = '#faf8f4';
  ctx.fillRect(0, 0, W, H);

  const LW = 82;       // label column width (px)
  const TW = W - LW;  // track width (px)
  const toX = pos => LW + (pos / len) * TW;

  // Track layout: [y-offset, height]
  const T = {
    gc:   [10,  58],
    cpg:  [76,  20],
    orf:  [104, 40],
    tata: [152, 20],
    bar:  [180, 18],
  };

  // Labels
  ctx.font = 'italic 10.5px Georgia, serif';
  ctx.fillStyle = '#aaa';
  ctx.textAlign = 'right';
  ctx.textBaseline = 'middle';
  [['GC %', T.gc], ['CpG islands', T.cpg], ['ORFs', T.orf], ['TATA boxes', T.tata]].forEach(
    ([lbl, [y, h]]) => ctx.fillText(lbl, LW - 6, y + h / 2)
  );

  // Track dividers
  ctx.strokeStyle = 'rgba(0,0,0,0.07)';
  ctx.lineWidth = 1;
  [T.cpg, T.orf, T.tata, T.bar].forEach(([y]) => {
    ctx.beginPath(); ctx.moveTo(LW, y - 3); ctx.lineTo(LW + TW, y - 3); ctx.stroke();
  });

  // Left axis border
  ctx.strokeStyle = 'rgba(0,0,0,0.15)';
  ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(LW, 6); ctx.lineTo(LW, T.bar[0] + T.bar[1]); ctx.stroke();

  // ── GC content (filled area + line) ───────────────────────────────────────
  const [gy, gh] = T.gc;
  const gcMin = 0.25, gcMax = 0.75;
  const gcY = gc => gy + gh - ((Math.max(gcMin, Math.min(gcMax, gc)) - gcMin) / (gcMax - gcMin)) * gh;

  ctx.beginPath();
  ctx.moveTo(LW, gy + gh);
  gcBins.forEach((gc, b) => ctx.lineTo(LW + (b / gcBins.length) * TW, gcY(gc)));
  ctx.lineTo(LW + TW, gy + gh);
  ctx.closePath();
  ctx.fillStyle = 'rgba(55,115,85,0.14)';
  ctx.fill();

  ctx.beginPath();
  gcBins.forEach((gc, b) => {
    const x = LW + (b / gcBins.length) * TW;
    b === 0 ? ctx.moveTo(x, gcY(gc)) : ctx.lineTo(x, gcY(gc));
  });
  ctx.strokeStyle = '#3a7a5a';
  ctx.lineWidth = 1;
  ctx.stroke();

  // 50% GC dashed reference line
  ctx.save();
  ctx.setLineDash([3, 3]);
  ctx.strokeStyle = 'rgba(0,0,0,0.11)';
  ctx.lineWidth = 0.5;
  const ref50y = gcY(0.5);
  ctx.beginPath(); ctx.moveTo(LW, ref50y); ctx.lineTo(LW + TW, ref50y); ctx.stroke();
  ctx.restore();

  // ── CpG islands ───────────────────────────────────────────────────────────
  const [cy, ch] = T.cpg;
  ctx.fillStyle = '#4a9a80';
  cpgIslands.forEach(({ start, end }) => {
    const x1 = toX(start), w = Math.max(toX(end) - x1, 1.5);
    ctx.fillRect(x1, cy + 2, w, ch - 4);
  });

  // ── ORFs: long ones (≥500 bp) in muted red, short in slate ───────────────
  const [oy, oh] = T.orf;
  orfs.forEach(({ start, end, long }) => {
    const x1 = toX(start), w = Math.max(toX(end) - x1, 1);
    if (long) {
      ctx.fillStyle = 'rgba(160,80,70,0.80)';
      ctx.fillRect(x1, oy + 2, w, oh - 4);
    } else {
      ctx.fillStyle = 'rgba(85,115,148,0.38)';
      ctx.fillRect(x1, oy + 10, w, oh - 20);
    }
  });

  // ── TATA boxes (tick marks) ───────────────────────────────────────────────
  const [ty, th] = T.tata;
  ctx.fillStyle = 'rgba(138,106,58,0.65)';
  tataBoxes.forEach(pos => ctx.fillRect(toX(pos), ty + 2, 1.5, th - 4));

  // ── Scale bar ─────────────────────────────────────────────────────────────
  const [by] = T.bar;
  const barMid = by + T.bar[1] / 2;
  ctx.strokeStyle = '#ccc';
  ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(LW, barMid); ctx.lineTo(LW + TW, barMid); ctx.stroke();

  ctx.font = '9.5px Georgia, serif';
  ctx.fillStyle = '#bbb';
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  for (let pos = 100000; pos < len; pos += 100000) {
    const x = toX(pos);
    ctx.strokeStyle = '#ccc';
    ctx.beginPath(); ctx.moveTo(x, barMid - 3); ctx.lineTo(x, barMid + 3); ctx.stroke();
    ctx.fillText(`${pos / 1000} kb`, x, barMid + 4);
  }
}

// ── Translate ─────────────────────────────────────────────────────────────────

// 64-entry table indexed by first*16 + second*4 + third  (A=0, C=1, G=2, T=3)
const CODON_TABLE = 'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF';

function translate(seq, start, end) {
  const aas = [];
  for (let i = start; i + 2 < end; i += 3) {
    const aa = CODON_TABLE[seq[i] * 16 + seq[i + 1] * 4 + seq[i + 2]];
    if (aa === '*') break;
    aas.push(aa);
  }
  return aas; // first element is always M (ATG start)
}

// ── Render protein ────────────────────────────────────────────────────────────

// Physicochemical classes
const AA_CLASS = {
  A:0, V:0, I:0, L:0, M:0, F:0, W:0, P:0,  // hydrophobic / nonpolar
  S:1, T:1, C:1, Y:1, N:1, Q:1,             // polar uncharged
  K:2, R:2, H:2,                            // basic (+)
  D:3, E:3,                                 // acidic (−)
  G:4,                                      // glycine (special)
};
const AA_COLORS  = ['#c09060', '#6a96b4', '#5070a8', '#b05050', '#7a9a70'];
const AA_LABELS  = ['Hydrophobic', 'Polar', 'Basic (+)', 'Acidic (−)', 'Gly'];

function renderProtein(canvas, W, H, aas) {
  const ctx = canvas.getContext('2d');
  ctx.fillStyle = '#faf8f4';
  ctx.fillRect(0, 0, W, H);

  const LW = 82;
  const TW = W - LW;
  const n = aas.length;
  const barY = 8, barH = 28;

  // Track label
  ctx.font = 'italic 10.5px Georgia, serif';
  ctx.fillStyle = '#aaa';
  ctx.textAlign = 'right';
  ctx.textBaseline = 'middle';
  ctx.fillText('protein', LW - 6, barY + barH / 2);

  // Left border
  ctx.strokeStyle = 'rgba(0,0,0,0.15)';
  ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(LW, 4); ctx.lineTo(LW, barY + barH + 4); ctx.stroke();

  // One rectangle per amino acid, colored by class
  aas.forEach((aa, i) => {
    const x1 = LW + (i / n) * TW;
    const x2 = LW + ((i + 1) / n) * TW;
    ctx.fillStyle = AA_COLORS[AA_CLASS[aa] ?? 0];
    ctx.fillRect(x1, barY, Math.max(x2 - x1, 0.5), barH);
  });
}

// ── FASTA helpers ─────────────────────────────────────────────────────────────

function seqToFasta(header, seq) {
  const NUC = 'ACGT';
  const lines = [`>${header}`];
  for (let i = 0; i < seq.length; i += 60) {
    let chunk = '';
    for (let j = i; j < Math.min(i + 60, seq.length); j++) chunk += NUC[seq[j]];
    lines.push(chunk);
  }
  return lines.join('\n');
}

function aasToFasta(header, aas) {
  const str = aas.join('');
  const lines = [`>${header}`];
  for (let i = 0; i < str.length; i += 60) lines.push(str.slice(i, i + 60));
  return lines.join('\n');
}

async function copyFasta(fasta, btn, label) {
  btn.textContent = 'Copying…';
  btn.disabled = true;
  await navigator.clipboard.writeText(fasta);
  btn.textContent = 'Copied!';
  setTimeout(() => { btn.textContent = label; btn.disabled = false; }, 2000);
}

// ── Bootstrap ─────────────────────────────────────────────────────────────────

setTimeout(() => {
  const LEN = 1_000_000;
  const t0 = performance.now();
  const seq = generate(LEN);
  const features = parse(seq);
  const ms = Math.round(performance.now() - t0);

  // ── DNA canvas ──────────────────────────────────────────────────────────────
  const dnaCanvas = document.getElementById('dna-canvas');
  const dpr = window.devicePixelRatio || 1;
  const cssW = dnaCanvas.offsetWidth;
  const dnaCssH = 208;

  dnaCanvas.width  = cssW * dpr;
  dnaCanvas.height = dnaCssH * dpr;
  dnaCanvas.style.height = dnaCssH + 'px';
  const dnaCtx = dnaCanvas.getContext('2d');
  dnaCtx.scale(dpr, dpr);
  render(dnaCanvas, cssW, dnaCssH, LEN, features);

  const { cpgIslands, orfs, tataBoxes, gcOverall } = features;
  const longOrfs = orfs.filter(o => o.long).length;

  document.getElementById('dna-title').textContent =
    `Simulated 1 Mb DNA fragment · generated in ${ms} ms · unique each visit`;
  document.getElementById('dna-stats').textContent =
    `${orfs.length.toLocaleString()} ORFs detected (${longOrfs} ≥ 500 bp) · ` +
    `${cpgIslands.length} CpG islands · ` +
    `${tataBoxes.length} putative TATA boxes · ` +
    `${(gcOverall * 100).toFixed(1)}% GC`;

  // ── Protein translation ─────────────────────────────────────────────────────
  const longestOrf = orfs.reduce((best, o) =>
    (o.end - o.start > best.end - best.start ? o : best), orfs[0]);

  const aas = longestOrf ? translate(seq, longestOrf.start, longestOrf.end) : [];

  if (aas.length > 0) {
    const protCanvas = document.getElementById('protein-canvas');
    const protCssH = 44;
    protCanvas.width  = cssW * dpr;
    protCanvas.height = protCssH * dpr;
    protCanvas.style.height = protCssH + 'px';
    const protCtx = protCanvas.getContext('2d');
    protCtx.scale(dpr, dpr);
    renderProtein(protCanvas, cssW, protCssH, aas);

    const hydrophobic = aas.filter(a => AA_CLASS[a] === 0).length;
    const charged     = aas.filter(a => AA_CLASS[a] === 2 || AA_CLASS[a] === 3).length;
    const polar       = aas.filter(a => AA_CLASS[a] === 1).length;

    document.getElementById('protein-title').textContent =
      `Longest ORF translation · ` +
      `${longestOrf.end - longestOrf.start} bp at position ` +
      `${longestOrf.start.toLocaleString()}`;

    document.getElementById('protein-stats').textContent =
      `${aas.length} amino acids · ` +
      `${((hydrophobic / aas.length) * 100).toFixed(0)}% hydrophobic · ` +
      `${((charged / aas.length) * 100).toFixed(0)}% charged · ` +
      `${((polar / aas.length) * 100).toFixed(0)}% polar`;

    // Legend color swatches (injected as HTML into a dedicated element)
    const legendEl = document.getElementById('protein-legend');
    if (legendEl) {
      legendEl.innerHTML = AA_LABELS.map((lbl, i) =>
        `<span class="swatch" style="background:${AA_COLORS[i]}"></span>${lbl}`
      ).join('');
    }
  }

  // ── Copy buttons ────────────────────────────────────────────────────────────
  document.getElementById('copy-dna').disabled = false;
  document.getElementById('copy-protein').disabled = false;

  document.getElementById('copy-dna').addEventListener('click', async function () {
    const btn = this;
    btn.textContent = 'Building…';
    btn.disabled = true;
    await new Promise(r => setTimeout(r, 0)); // yield so label updates
    const fasta = seqToFasta(
      `sim_fragment | 1,000,000 bp | generated ${new Date().toISOString().slice(0,10)}`,
      seq
    );
    await copyFasta(fasta, btn, 'Copy DNA (FASTA)');
  });

  document.getElementById('copy-protein').addEventListener('click', async function () {
    const fasta = aasToFasta(
      `sim_protein | ${aas.length} aa | orf ${longestOrf.start}..${longestOrf.end}`,
      aas
    );
    await copyFasta(fasta, this, 'Copy protein (FASTA)');
  });
}, 0);
