import { track } from "./analytics.js";
import {
  loadGene,
  loadIndex,
  formatQ,
  formatNum,
  formatCoord,
  chromLabel,
  zenodoChromArchive,
  dataUrl,
  wireHomeSearch,
  goGene,
} from "./app.js";
import { renderExonMap, renderIntronView, intronFacts } from "./plots.js";

const params = new URLSearchParams(window.location.search);
const symbol = params.get("g") || "";
let sortKey = "q";
let sortDir = "asc";
let selectedId = params.get("hs") || params.get("iv") || null;
let geneDoc = null;

function $(id) {
  return document.getElementById(id);
}

function setSelected(iv, push = true) {
  if (!iv || !geneDoc) return;
  selectedId = iv.id;
  for (const row of document.querySelectorAll("#interval-table tbody tr")) {
    row.classList.toggle("selected", row.dataset.id === iv.id);
  }
  $("intron-view-title").textContent = `${iv.id} in intron ${iv.intron}`;
  renderExonMap($("exon-map"), {
    exons: geneDoc.exons,
    intervals: geneDoc.intervals,
    selected: iv,
    strand: geneDoc.strand,
  });
  renderIntronView($("intron-view"), {
    gene: geneDoc.gene,
    interval: iv,
    variants: geneDoc.variants,
    strand: geneDoc.strand,
    chrom: geneDoc.chrom || iv.chrom,
  });
  const facts = intronFacts(iv, geneDoc.strand);
  $("fact-donor").textContent = facts.donor == null ? "—" : `+${Number(facts.donor).toLocaleString()} bp`;
  $("fact-acceptor").textContent = facts.acceptor == null ? "—" : `−${Number(facts.acceptor).toLocaleString()} bp`;
  $("fact-intron-len").textContent = facts.intronLength == null
    ? "—"
    : `${Number(facts.intronLength).toLocaleString()} bp`;
  $("fact-ratio").textContent = facts.ratio;

  if (push) {
    const url = new URL(window.location.href);
    url.searchParams.set("g", geneDoc.gene);
    url.searchParams.set("iv", iv.id);
    history.pushState({}, "", url);
    track(`/hotspot/${geneDoc.gene}/${iv.id}`, `${geneDoc.gene} ${iv.id}`);
  }
}

function sortedIntervals(intervals) {
  const arr = intervals.slice();
  const dir = sortDir === "asc" ? 1 : -1;
  arr.sort((a, b) => {
    let av = a[sortKey];
    let bv = b[sortKey];
    if (sortKey === "top_in_intron") {
      av = a.top_in_intron ? 1 : 0;
      bv = b.top_in_intron ? 1 : 0;
    }
    if (av == null && bv == null) return 0;
    if (av == null) return 1;
    if (bv == null) return -1;
    if (av < bv) return -1 * dir;
    if (av > bv) return 1 * dir;
    return 0;
  });
  return arr;
}

function renderTable() {
  const tbody = $("interval-table").querySelector("tbody");
  tbody.innerHTML = "";
  const rows = sortedIntervals(geneDoc.intervals || []);
  if (!rows.length) {
    $("interval-empty").hidden = false;
    return;
  }
  $("interval-empty").hidden = true;
  for (const iv of rows) {
    const tr = document.createElement("tr");
    tr.dataset.id = iv.id;
    if (iv.id === selectedId) tr.classList.add("selected");
    const top = iv.top_in_intron
      ? `<span class="chip chip-orange" title="Best min_BH_q in this intron">top</span>`
      : "";
    tr.innerHTML = `
      <td><button type="button" class="iv-link">${iv.id}</button> ${top}</td>
      <td>${iv.intron ?? "—"}</td>
      <td class="mono">${formatCoord(iv.chrom, iv.start, iv.end)}</td>
      <td>${iv.length ?? "—"}</td>
      <td class="mono">${formatQ(iv.q)}</td>
      <td>${iv.n_spliceogenic ?? "—"}</td>
      <td>${formatNum(iv.concentration, 3)}</td>
      <td>${formatNum(iv.max_score, 3)}</td>`;
    tr.querySelector("button").addEventListener("click", () => setSelected(iv));
    tr.addEventListener("click", (ev) => {
      if (ev.target.closest("button")) return;
      setSelected(iv);
    });
    tbody.appendChild(tr);
  }
}

async function init() {
  const compact = $("gene-search-form");
  if (compact) {
    wireHomeSearch({
      form: compact,
      input: $("gene-search"),
      msg: $("gene-msg"),
      track,
    });
  }

  if (!symbol) {
    $("gene-msg").textContent = "Add a gene symbol, for example gene.html?g=COL4A5";
    return;
  }

  try {
    geneDoc = await loadGene(symbol);
  } catch (err) {
    if (err.code === "NOT_FOUND") {
      const idx = await loadIndex();
      const { matchGenes } = await import("./app.js");
      const near = matchGenes(idx, symbol, 5);
      $("gene-msg").textContent = `“${symbol}” is not in this dataset.`
        + (near.length ? ` Nearest: ${near.join(", ")}.` : "");
      track(`/search-miss/${String(symbol).toUpperCase().slice(0, 20)}`);
      return;
    }
    $("gene-msg").textContent = err.message || "Could not load gene.";
    return;
  }

  document.title = `${geneDoc.gene} · ISHMap`;
  $("gene-symbol").textContent = geneDoc.gene;
  $("gene-transcript").textContent = geneDoc.transcript || geneDoc.summary?.transcript || "—";
  const chrom = geneDoc.chrom || chromLabel(geneDoc.summary?.chrom);
  $("gene-locus").textContent = `chr${chrom} · ${geneDoc.strand || "?"} strand`;
  $("gene-exons").textContent = geneDoc.exons?.length
    ? `${geneDoc.exons.length} exons`
    : `${geneDoc.summary?.n_introns ?? "—"} introns`;

  const nIv = geneDoc.intervals?.length || 0;
  const nTop = (geneDoc.intervals || []).filter((x) => x.top_in_intron).length;
  const bestQ = geneDoc.intervals?.[0]?.q;
  $("stat-intervals").textContent = String(nIv);
  $("stat-top").textContent = String(nTop);
  $("stat-q").textContent = formatQ(bestQ);
  $("stat-score").textContent = formatNum(geneDoc.summary?.max_atlas_score, 3);

  const stem = geneDoc._index?.f || geneDoc.gene;
  const dl = $("dl-gene-csv");
  dl.href = dataUrl(`downloads/${stem}_intervals.csv`);
  dl.download = `${stem}_intervals.csv`;
  dl.addEventListener("click", () => track(`/download/${geneDoc.gene}`));

  const ver = "v1.0";
  const archive = zenodoChromArchive(chrom, ver);
  const zen = $("dl-zenodo");
  zen.textContent = `Full chr${chromLabel(chrom)} archive on Zenodo`;
  zen.href = "#release-note";
  zen.title = "Zenodo full-release DOI not reserved yet";
  zen.addEventListener("click", (ev) => {
    ev.preventDefault();
    $("release-note").hidden = false;
    track(`/download/zenodo/${archive}`);
  });

  for (const th of document.querySelectorAll("#interval-table th[data-key]")) {
    th.addEventListener("click", () => {
      const key = th.dataset.key;
      if (sortKey === key) sortDir = sortDir === "asc" ? "desc" : "asc";
      else {
        sortKey = key;
        sortDir = key === "q" ? "asc" : "desc";
      }
      renderTable();
    });
  }

  renderTable();
  const pick = (geneDoc.intervals || []).find((x) => x.id === selectedId)
    || (geneDoc.intervals || []).find((x) => x.default)
    || (geneDoc.intervals || [])[0];
  if (pick) setSelected(pick, false);
  else {
    $("interval-empty").hidden = false;
    $("intron-view").innerHTML = `<p class="muted">No interior intervals passed the threshold for this gene.</p>`;
  }

  track(`/gene/${geneDoc.gene}`, geneDoc.gene);
  window.addEventListener("popstate", () => {
    const p = new URLSearchParams(window.location.search);
    const id = p.get("iv") || p.get("hs");
    const iv = (geneDoc.intervals || []).find((x) => x.id === id);
    if (iv) setSelected(iv, false);
  });
}

init();
