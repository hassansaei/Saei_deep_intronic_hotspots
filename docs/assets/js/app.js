/** Shared data loading and gene search. */

const geneCache = new Map();
let indexPromise = null;

export function dataUrl(rel) {
  return new URL(`../../data/${rel}`, import.meta.url).href;
}

export async function loadIndex() {
  if (!indexPromise) {
    indexPromise = fetch(dataUrl("index.json")).then((r) => {
      if (!r.ok) throw new Error(`Could not load index.json (${r.status})`);
      return r.json();
    });
  }
  return indexPromise;
}

export async function loadGene(symbol) {
  const key = String(symbol || "").toUpperCase();
  if (!key) throw new Error("Missing gene symbol");
  if (geneCache.has(key)) return geneCache.get(key);
  const idx = await loadIndex();
  const entry = idx.genes[key] || idx.genes[symbol];
  if (!entry) {
    const err = new Error(`Gene not found: ${symbol}`);
    err.code = "NOT_FOUND";
    throw err;
  }
  const stem = entry.f || key;
  const url = dataUrl(`genes/${stem}.json`);
  const res = await fetch(url);
  if (!res.ok) throw new Error(`Could not load ${stem}.json (${res.status})`);
  const doc = await res.json();
  doc._index = entry;
  geneCache.set(key, doc);
  return doc;
}

export function matchGenes(index, query, limit = 8) {
  const q = String(query || "").trim().toUpperCase();
  if (!q) return [];
  const keys = Object.keys(index.genes);
  const prefix = [];
  const mid = [];
  for (const g of keys) {
    const u = g.toUpperCase();
    if (u.startsWith(q)) prefix.push(g);
    else if (u.includes(q)) mid.push(g);
  }
  prefix.sort();
  mid.sort();
  return prefix.concat(mid).slice(0, limit);
}

export function formatQ(q) {
  if (q == null || Number.isNaN(Number(q))) return "—";
  const n = Number(q);
  if (n === 0) return "0";
  if (n < 0.001) return n.toExponential(2);
  return n.toPrecision(3);
}

export function formatNum(n, digits = 3) {
  if (n == null || Number.isNaN(Number(n))) return "—";
  const x = Number(n);
  if (Number.isInteger(x)) return x.toLocaleString();
  return x.toFixed(digits);
}

export function formatCoord(chrom, start, end) {
  const c = String(chrom || "").replace(/^chr/i, "");
  if (start == null || end == null) return `chr${c}`;
  return `chr${c}:${Number(start).toLocaleString()}–${Number(end).toLocaleString()}`;
}

export function chromLabel(chrom) {
  return String(chrom || "").replace(/^chr/i, "");
}

export function zenodoChromArchive(chrom, version = "v1.0") {
  const c = chromLabel(chrom) || "unplaced";
  const name = c === "unplaced" ? `genes_unplaced_${version}.zip` : `genes_chr${c}_${version}.zip`;
  return name;
}

export function computeHomeStats(index) {
  let nIntervals = 0;
  let nSplice = 0;
  let nTop = 0;
  for (const g of Object.values(index.genes)) {
    nIntervals += Number(g.h) || 0;
    nTop += Number(g.t) || 0;
    const s = g.s || {};
    nSplice += Number(s.n_spliceogenic_interior) || 0;
  }
  return {
    nGenes: index.meta.n_genes || Object.keys(index.genes).length,
    nIntervals,
    nTop,
    nSplice,
    assembly: index.meta.assembly || "GRCh38",
  };
}

export function wireHomeSearch({ form, input, msg, track }) {
  form.addEventListener("submit", (ev) => {
    ev.preventDefault();
    const q = input.value.trim();
    if (!q) return;
    goGene(q, msg, track);
  });

  let list = null;
  input.addEventListener("input", async () => {
    const idx = await loadIndex();
    const hits = matchGenes(idx, input.value, 8);
    if (!list) {
      list = document.createElement("datalist");
      list.id = "gene-suggestions";
      document.body.appendChild(list);
      input.setAttribute("list", list.id);
    }
    list.innerHTML = hits.map((g) => `<option value="${g}"></option>`).join("");
  });
}

export async function goGene(symbol, msgEl, track) {
  const idx = await loadIndex();
  const hits = matchGenes(idx, symbol, 5);
  const exact = idx.genes[symbol] || idx.genes[String(symbol).toUpperCase()];
  if (!exact) {
    if (msgEl) {
      const near = hits.length ? ` Nearest: ${hits.join(", ")}.` : "";
      msgEl.textContent = `“${symbol}” is not in this dataset.${near}`;
    }
    if (track) track(`/search-miss/${String(symbol).toUpperCase().slice(0, 20)}`);
    return;
  }
  const key = exact === idx.genes[symbol] ? symbol : String(symbol).toUpperCase();
  const canon = Object.keys(idx.genes).find((g) => idx.genes[g] === exact) || key;
  window.location.href = `gene.html?g=${encodeURIComponent(canon)}`;
}
