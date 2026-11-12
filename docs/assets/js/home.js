import { track } from "./analytics.js";
import {
  loadIndex,
  computeHomeStats,
  wireHomeSearch,
  goGene,
} from "./app.js";

function $(id) {
  return document.getElementById(id);
}

async function init() {
  wireHomeSearch({
    form: $("search-form"),
    input: $("gene-input"),
    msg: $("search-msg"),
    track,
  });

  for (const btn of document.querySelectorAll("#example-genes .chip")) {
    btn.addEventListener("click", () => {
      const g = btn.dataset.gene;
      $("gene-input").value = g;
      goGene(g, $("search-msg"), track);
    });
  }

  try {
    const index = await loadIndex();
    const stats = computeHomeStats(index);
    $("stat-genes").textContent = stats.nGenes.toLocaleString();
    $("stat-intervals").textContent = stats.nIntervals.toLocaleString();
    $("stat-splice").textContent = stats.nSplice.toLocaleString();
    $("stat-assembly").textContent = `${stats.assembly} · scored RefSeq`;
    $("meta-built").textContent = index.meta.built || "";
    const ver = index.meta.version || "v1.0";
    $("meta-version").textContent = String(ver).replace(/^v/i, "");
    const citeVer = $("cite-version");
    if (citeVer) citeVer.textContent = ver;
  } catch (err) {
    $("search-msg").textContent = err.message || "Could not load index.json";
  }
}

init();
