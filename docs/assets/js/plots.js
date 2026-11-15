/** Inline SVG exon map and intron interval view. */

function el(name, attrs = {}, kids = []) {
  const n = document.createElementNS("http://www.w3.org/2000/svg", name);
  for (const [k, v] of Object.entries(attrs)) {
    if (v != null) n.setAttribute(k, String(v));
  }
  for (const c of kids) n.appendChild(typeof c === "string" ? document.createTextNode(c) : c);
  return n;
}

function scoreColor(score) {
  if (score == null) return "var(--c-below)";
  if (score >= 0.8) return "var(--c-score-high)";
  if (score >= 0.5) return "var(--c-score-mid)";
  return "var(--c-score-low)";
}

export function renderExonMap(container, { exons, intervals, selected, strand }) {
  container.innerHTML = "";
  if (!exons || !exons.length) {
    container.innerHTML = `<p class="muted">No exon coordinates for this transcript. Interval view uses intron bounds when available.</p>`;
    return;
  }
  const W = 1000;
  const H = 64;
  const pad = 16;
  const nEx = exons.length;
  const y = 28;
  const minExonW = 3;

  // Exons are [start, end] on GRCh38 in genomic order; scale the full span linearly.
  const geneStart = Math.min(...exons.map(([a]) => a));
  const geneEnd = Math.max(...exons.map(([, b]) => b));
  const span = Math.max(geneEnd - geneStart, 1);
  const usable = W - 2 * pad;
  const xOf = (pos) => pad + ((pos - geneStart) / span) * usable;

  const svg = el("svg", {
    viewBox: `0 0 ${W} ${H}`,
    role: "img",
    "aria-label": `Transcript exon map with ${nEx} exons scaled to genomic length`,
  });

  const byIntron = new Map();
  for (const iv of intervals || []) {
    if (iv.intron == null) continue;
    if (!byIntron.has(iv.intron) || iv.top_in_intron) byIntron.set(iv.intron, iv);
  }

  // Intron backbone first so exons sit on top.
  for (let i = 0; i < nEx - 1; i += 1) {
    const [, eEnd] = exons[i];
    const [nStart] = exons[i + 1];
    svg.appendChild(el("line", {
      x1: xOf(eEnd), y1: y, x2: xOf(nStart), y2: y,
      stroke: "var(--c-intron)", "stroke-width": 2,
    }));
  }

  for (let i = 0; i < nEx; i += 1) {
    const [a, b] = exons[i];
    const x0 = xOf(a);
    const x1 = xOf(b);
    const w = Math.max(x1 - x0, minExonW);
    svg.appendChild(el("rect", {
      x: x0, y: y - 9, width: w, height: 18, rx: 1.5,
      fill: "var(--c-primary)",
    }));
  }

  for (let i = 0; i < nEx - 1; i += 1) {
    const [, eEnd] = exons[i];
    const [nStart] = exons[i + 1];
    const intronNum = strand === "-" ? nEx - 1 - i : i + 1;
    const iv = byIntron.get(intronNum);
    if (!iv) continue;
    const mid = iv.start != null && iv.end != null
      ? (iv.start + iv.end) / 2
      : (eEnd + nStart) / 2;
    const cx = xOf(mid);
    const sel = selected && selected.intron === intronNum;
    svg.appendChild(el("polygon", {
      points: `${cx - 5},${y + 14} ${cx + 5},${y + 14} ${cx},${y + 4}`,
      fill: sel ? "var(--c-hotspot)" : "var(--c-hotspot-soft)",
    }));
    if (sel) {
      svg.appendChild(el("text", {
        x: cx, y: H - 6, "text-anchor": "middle",
        fill: "var(--c-hotspot-ink)",
        "font-family": "var(--font-mono)",
        "font-size": "11",
      }, [`i${intronNum}`]));
    }
  }
  container.appendChild(svg);
}

export function renderIntronView(container, { gene, interval, variants, strand, chrom }) {
  container.innerHTML = "";
  if (!interval) {
    container.innerHTML = `<p class="muted">No interior interval selected.</p>`;
    return;
  }

  const W = 1000;
  const H = 250;
  const left = 56;
  const right = 24;
  const top = 28;
  const bottom = 40;
  const plotW = W - left - right;
  const plotH = H - top - bottom;

  let iStart = interval.intron_start;
  let iEnd = interval.intron_end;
  if (iStart == null || iEnd == null || iStart >= iEnd) {
    const pad = 2000;
    iStart = interval.start - pad;
    iEnd = interval.end + pad;
  }
  const span = Math.max(iEnd - iStart, 1);
  const xOf = (pos) => left + ((pos - iStart) / span) * plotW;
  const yOf = (score) => top + plotH - Math.min(Math.max(score, 0), 2) / 2 * plotH;

  const svg = el("svg", {
    viewBox: `0 0 ${W} ${H}`,
    role: "img",
    "aria-label": `${interval.id} in intron ${interval.intron}: spliceogenic variants and interval bounds`,
  });

  // axes
  svg.appendChild(el("line", {
    x1: left, y1: top + plotH, x2: left + plotW, y2: top + plotH,
    stroke: "var(--c-line-strong)",
  }));
  svg.appendChild(el("line", {
    x1: left, y1: top, x2: left, y2: top + plotH,
    stroke: "var(--c-line-strong)",
  }));
  for (const s of [0, 0.56, 1, 2]) {
    const y = yOf(s);
    svg.appendChild(el("line", {
      x1: left - 4, y1: y, x2: left, y2: y,
      stroke: "var(--c-line-strong)",
    }));
    svg.appendChild(el("text", {
      x: left - 8, y: y + 3, "text-anchor": "end",
      fill: "var(--c-ink-3)", "font-size": "11", "font-family": "var(--font-mono)",
    }, [String(s)]));
  }

  // flanking exon stubs
  const stub = 18;
  svg.appendChild(el("rect", {
    x: left - stub, y: top + plotH - 8, width: stub, height: 8, fill: "var(--c-primary)",
  }));
  svg.appendChild(el("rect", {
    x: left + plotW, y: top + plotH - 8, width: stub, height: 8, fill: "var(--c-primary)",
  }));
  svg.appendChild(el("line", {
    x1: left, y1: top + plotH - 4, x2: left + plotW, y2: top + plotH - 4,
    stroke: "var(--c-intron)", "stroke-width": 2,
  }));

  const donorLeft = strand !== "-";
  svg.appendChild(el("text", {
    x: left, y: H - 12, "text-anchor": "start",
    fill: "var(--c-ink-3)", "font-size": "11",
  }, [donorLeft ? "donor" : "acceptor"]));
  svg.appendChild(el("text", {
    x: left + plotW, y: H - 12, "text-anchor": "end",
    fill: "var(--c-ink-3)", "font-size": "11",
  }, [donorLeft ? "acceptor" : "donor"]));

  // interval shade
  const x0 = xOf(interval.start);
  const x1 = xOf(interval.end);
  svg.appendChild(el("rect", {
    x: Math.min(x0, x1),
    y: top,
    width: Math.max(Math.abs(x1 - x0), 2),
    height: plotH,
    fill: "var(--c-hotspot-fill)",
    opacity: "0.55",
  }));
  svg.appendChild(el("line", {
    x1: x0, y1: top, x2: x0, y2: top + plotH,
    stroke: "var(--c-hotspot)", "stroke-dasharray": "4 3",
  }));
  svg.appendChild(el("line", {
    x1: x1, y1: top, x2: x1, y2: top + plotH,
    stroke: "var(--c-hotspot)", "stroke-dasharray": "4 3",
  }));

  const g = el("g", {});
  for (const v of variants || []) {
    const [pos, score, ref, alt] = v;
    if (pos < iStart || pos > iEnd) continue;
    const cx = xOf(pos);
    const cy = yOf(score);
    const color = scoreColor(score);
    const circ = el("circle", {
      cx, cy, r: 3.2,
      fill: color,
      opacity: "0.9",
    });
    const donorDist = strand === "-"
      ? (iEnd - pos + 1)
      : (pos - iStart + 1);
    const accDist = strand === "-"
      ? (pos - iStart + 1)
      : (iEnd - pos + 1);
    const c = String(chrom || "").replace(/^chr/i, "chr");
    circ.appendChild(el("title", {}, [
      `${c}:${Number(pos).toLocaleString()} ${ref || "?"}>${alt || "?"} · score ${Number(score).toFixed(3)} · donor +${donorDist.toLocaleString()} / acceptor −${accDist.toLocaleString()}`,
    ]));
    g.appendChild(circ);
  }
  svg.appendChild(g);

  const label = `${interval.id} in intron ${interval.intron}`;
  svg.appendChild(el("text", {
    x: left, y: 16, fill: "var(--c-ink)", "font-size": "13", "font-weight": "600",
  }, [label]));

  container.appendChild(svg);
}

export function intronFacts(interval, strand) {
  if (!interval) {
    return { donor: "—", acceptor: "—", intronLength: null, ratio: "—" };
  }
  if (interval.intron_start == null || interval.intron_end == null) {
    return {
      donor: "—",
      acceptor: "—",
      intronLength: null,
      ratio: interval.n_spliceogenic_interior
        ? `${interval.n_spliceogenic} / ${interval.n_spliceogenic_interior}`
        : `${interval.n_spliceogenic ?? "—"}`,
    };
  }
  const centre = Math.round((interval.start + interval.end) / 2);
  const donor = strand === "-"
    ? interval.intron_end - centre + 1
    : centre - interval.intron_start + 1;
  const acceptor = strand === "-"
    ? centre - interval.intron_start + 1
    : interval.intron_end - centre + 1;
  const intronLength = interval.intron_end - interval.intron_start + 1;
  const ratio = interval.n_spliceogenic_interior
    ? `${interval.n_spliceogenic} / ${interval.n_spliceogenic_interior}`
    : `${interval.n_spliceogenic ?? "—"}`;
  return { donor, acceptor, intronLength, ratio };
}
