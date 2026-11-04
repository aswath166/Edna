(function () {
  "use strict";

  const $ = (sel, root) => (root || document).querySelector(sel);
  const $$ = (sel, root) => Array.from((root || document).querySelectorAll(sel));
  const reduceMotion = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  const PALETTE = ["#34f5c5", "#5ab8ff", "#b48cff", "#8dff7a", "#ff7ad9", "#6cf0ff", "#ffe066"];
  const NOVEL = "#ffab40";
  const SISTER = "#5ab8ff";
  const ARTIFACT = "#ff6b6b";
  const ANCHOR_COLOR = "#ffd166";
  const BASE_COLORS = ["#34f5c5", "#ffab40", "#5ab8ff", "#b48cff"];

  function esc(text) {
    const div = document.createElement("div");
    div.textContent = text == null ? "" : String(text);
    return div.innerHTML;
  }

  function hexA(hex, alpha) {
    const n = parseInt(hex.slice(1), 16);
    return `rgba(${(n >> 16) & 255}, ${(n >> 8) & 255}, ${n & 255}, ${alpha})`;
  }

  /* ---------------------------------------------------------------
     Masthead Double Helix Canvas
  --------------------------------------------------------------- */
  function startHelix() {
    const canvas = $("#helix");
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    const W = canvas.width;
    const H = canvas.height;
    const rungs = 11;

    function draw(t) {
      ctx.clearRect(0, 0, W, H);
      const phase = t / 1400;
      for (let i = 0; i < rungs; i++) {
        const y = 10 + (i * (H - 20)) / (rungs - 1);
        const a = phase + i * 0.62;
        const x1 = W / 2 + Math.sin(a) * (W / 2 - 14);
        const x2 = W / 2 - Math.sin(a) * (W / 2 - 14);
        const depth = (Math.cos(a) + 1) / 2;
        ctx.strokeStyle = `rgba(127, 165, 155, ${0.25 + depth * 0.35})`;
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(x1, y);
        ctx.lineTo(x2, y);
        ctx.stroke();

        const c1 = BASE_COLORS[i % 4];
        const c2 = BASE_COLORS[(i + 2) % 4];
        [[x1, c1, depth], [x2, c2, 1 - depth]].forEach(([x, c, d]) => {
          ctx.fillStyle = c;
          ctx.globalAlpha = 0.45 + d * 0.55;
          ctx.shadowColor = c;
          ctx.shadowBlur = 10;
          ctx.beginPath();
          ctx.arc(x, y, 3 + d * 2.2, 0, Math.PI * 2);
          ctx.fill();
        });
        ctx.globalAlpha = 1;
        ctx.shadowBlur = 0;
      }
    }

    if (reduceMotion) {
      draw(0);
      return;
    }
    function loop(t) {
      draw(t);
      requestAnimationFrame(loop);
    }
    requestAnimationFrame(loop);
  }

  /* ---------------------------------------------------------------
     Upload & Preset Execution
  --------------------------------------------------------------- */
  let currentDataset = null;
  let selectedFile = null;

  const dropZone = $("#drop-zone");
  const fileInput = $("#file-input");
  const fileName = $("#file-name");
  const analyzeBtn = $("#analyze-btn");
  const statusEl = $("#status-message");
  const form = $("#upload-form");

  function buildTicker() {
    const track = $("#ticker");
    if (!track) return;
    let html = "";
    for (let i = 0; i < 160; i++) {
      const b = "ACGT"[Math.floor(Math.random() * 4)];
      html += `<span class="b-${b.toLowerCase()}">${b}</span>`;
    }
    track.innerHTML = html + html;
  }

  function setFile(file) {
    selectedFile = file || null;
    if (selectedFile) {
      const kb = Math.max(1, Math.round(selectedFile.size / 1024));
      fileName.textContent = `${selectedFile.name} (${kb} KB)`;
      analyzeBtn.disabled = false;
      setStatus("");
    } else {
      fileName.textContent = "no file chosen";
      analyzeBtn.disabled = true;
    }
  }

  function setStatus(text, isError) {
    statusEl.textContent = text;
    statusEl.classList.toggle("error", !!isError);
  }

  dropZone.addEventListener("click", (e) => {
    if (e.target.closest("button")) return;
    fileInput.click();
  });
  fileInput.addEventListener("change", () => setFile(fileInput.files[0]));

  ["dragenter", "dragover"].forEach((ev) =>
    dropZone.addEventListener(ev, (e) => {
      e.preventDefault();
      dropZone.classList.add("dragover");
    })
  );
  ["dragleave", "drop"].forEach((ev) =>
    dropZone.addEventListener(ev, (e) => {
      e.preventDefault();
      dropZone.classList.remove("dragover");
    })
  );
  dropZone.addEventListener("drop", (e) => {
    if (e.dataTransfer && e.dataTransfer.files.length) setFile(e.dataTransfer.files[0]);
  });

  // 1-Click Preset Buttons
  $$(".btn-preset").forEach((btn) => {
    btn.addEventListener("click", () => {
      const presetId = btn.dataset.preset;
      runAnalysis({ preset_id: presetId }, btn.textContent.trim());
    });
  });

  form.addEventListener("submit", (e) => {
    e.preventDefault();
    if (!selectedFile) {
      setStatus("Please select a file first.", true);
      return;
    }
    const formData = new FormData();
    formData.append("file", selectedFile);
    runAnalysis(formData, selectedFile.name);
  });

  const STAGES = [
    "Profiling nucleotide quality & GC%...",
    "Extracting canonical k-mer manifolds...",
    "HDBSCAN density clustering...",
    "Evaluating reference centroid divergence...",
    "Quantifying unaligned dark taxa...",
    "Computing Chao1 & Shannon biodiversity...",
  ];

  async function runAnalysis(bodyData, label) {
    analyzeBtn.disabled = true;
    dropZone.classList.add("busy");
    let stage = 0;
    setStatus(STAGES[0]);
    const timer = setInterval(() => {
      stage = (stage + 1) % STAGES.length;
      setStatus(STAGES[stage]);
    }, 1100);

    try {
      let response;
      if (bodyData instanceof FormData) {
        response = await fetch("/api/analyze", { method: "POST", body: bodyData });
      } else {
        const fd = new FormData();
        Object.entries(bodyData).forEach(([k, v]) => fd.append(k, v));
        response = await fetch("/api/analyze", { method: "POST", body: fd });
      }

      const payload = await response.json();
      if (!response.ok || !payload || payload.error) {
        setStatus(payload?.error || `Server error (${response.status})`, true);
        return;
      }

      currentDataset = payload;
      renderAll(payload, label);
      setStatus(`Analysis complete: ${payload.total_sequences} reads processed. ${payload.reference_taxa_loaded} reference taxa loaded.`);
    } catch (err) {
      setStatus("Could not reach backend. Verify server is running on port 5000.", true);
    } finally {
      clearInterval(timer);
      dropZone.classList.remove("busy");
      analyzeBtn.disabled = !selectedFile;
    }
  }

  /* ---------------------------------------------------------------
     Tab Navigation
  --------------------------------------------------------------- */
  $$(".tab-btn").forEach((btn) => {
    btn.addEventListener("click", () => {
      $$(".tab-btn").forEach((b) => b.classList.remove("active"));
      $$(".tab-content").forEach((c) => c.classList.remove("active"));
      btn.classList.add("active");
      const target = $("#" + btn.dataset.tab);
      if (target) target.classList.add("active");

      // Redraw charts if switching to their tab
      if (btn.dataset.tab === "tab-qc") drawQCCharts(currentDataset?.qc_metrics);
      if (btn.dataset.tab === "tab-bio") drawWhittakerChart(currentDataset?.biodiversity_indices?.rank_abundance);
      if (btn.dataset.tab === "tab-map" && mapSonar) mapSonar.requestDraw();
    });
  });

  /* ---------------------------------------------------------------
     Master Render
  --------------------------------------------------------------- */
  let rows = [];
  let activeId = null;
  let filterMode = "all";
  let mapSonar = null;

  function renderAll(data, label) {
    rows = data.sequence_results || [];
    activeId = null;
    filterMode = "all";

    $("#results-section").classList.remove("hidden");
    $("#sample-title").textContent = label || "Environmental DNA Sample";
    $("#run-summary").textContent = `${data.total_sequences} Reads &bull; ${data.reference_taxa_loaded} Reference Anchors &bull; ${data.biodiversity_indices.richness} cOTU Lineages`;

    // 1. Overview & Sonar Map
    renderMetrics(data);
    const groups = buildGroups();
    renderComposition(groups);
    renderLegend(groups);

    if (mapSonar) mapSonar.destroy();
    mapSonar = createSonar($("#map"), $("#map-tip"), rows, data.reference_anchors, groups);

    // 2. QC Tab
    renderQCSummary(data.qc_metrics);

    // 3. Biodiversity Tab
    renderClustersTable(data.cluster_summaries);

    // 4. Benchmark Tab
    renderBenchmark(data.benchmark, data.biodiversity_indices);

    // 5. Read Inspector
    renderFilters();
    renderTable();

    // 6. Thesis Summary
    renderThesisDossier(data);

    // Scroll smoothly to results
    $("#results-section").scrollIntoView({ behavior: "smooth", block: "start" });
  }

  /* ---------------------------------------------------------------
     Biodiversity Metrics Cards
  --------------------------------------------------------------- */
  function renderMetrics(data) {
    const idx = data.biodiversity_indices;
    const novelCount = rows.filter((r) => r.category === "novel_taxon" || r.category === "sister_taxon").length;
    const tiles = [
      { label: "Observed Richness (S)", value: idx.richness, dec: 0, cap: "distinct taxonomic groups" },
      { label: "Chao1 Estimated (S_chao1)", value: idx.chao1_richness, dec: 1, cap: "asymptotic true richness", novel: true },
      { label: "Shannon-Wiener (H')", value: idx.shannon_diversity, dec: 3, cap: "entropy, higher is richer" },
      { label: "Simpson (1 - D)", value: idx.simpson_diversity, dec: 3, cap: "probability reads differ" },
      { label: "Pielou Evenness (J')", value: idx.pielou_evenness, dec: 3, cap: "community equitability [0, 1]" },
      { label: "Berger-Parker Dominance", value: idx.berger_parker_dominance, dec: 3, cap: "dominant species share" },
      { label: "Dark / Novel Taxa", value: novelCount, dec: 0, cap: "reads outside reference catalog", novel: true },
    ];

    const wrap = $("#metrics");
    wrap.innerHTML = "";
    tiles.forEach((t) => {
      const div = document.createElement("div");
      div.className = "metric" + (t.novel ? " novel" : "");
      div.innerHTML = `<div class="m-label">${esc(t.label)}</div><div class="m-value">0</div><div class="m-cap">${esc(t.cap)}</div>`;
      wrap.appendChild(div);
      animateNumber($(".m-value", div), Number(t.value), t.dec);
    });
  }

  function animateNumber(el, target, dec) {
    if (reduceMotion) {
      el.textContent = target.toFixed(dec);
      return;
    }
    const start = performance.now();
    const dur = 800;
    function frame(now) {
      const p = Math.min(1, (now - start) / dur);
      const eased = 1 - Math.pow(1 - p, 3);
      el.textContent = (target * eased).toFixed(dec);
      if (p < 1) requestAnimationFrame(frame);
    }
    requestAnimationFrame(frame);
  }

  /* ---------------------------------------------------------------
     Group & Taxonomy Composition
  --------------------------------------------------------------- */
  function groupOf(r) {
    if (r.category === "known") {
      return { key: r.closest_known_taxon, label: r.closest_known_taxon, category: "known", color: "#34f5c5" };
    }
    if (r.category === "sister_taxon") {
      return { key: "sister-" + r.cluster, label: `Sister Lineage (${r.closest_known_taxon})`, category: "sister_taxon", color: SISTER };
    }
    if (r.category === "novel_taxon") {
      return { key: "novel-" + r.cluster, label: `Candidate Dark Taxon (cOTU ${r.cluster})`, category: "novel_taxon", color: NOVEL };
    }
    return { key: "artifact", label: "Stray / Chimeric Reads", category: "artifact_outlier", color: ARTIFACT };
  }

  function buildGroups() {
    const groups = new Map();
    let colIdx = 0;
    rows.forEach((r) => {
      const g = groupOf(r);
      r._group = g.key;
      if (!groups.has(g.key)) {
        let color = g.color;
        if (g.category === "known") {
          color = PALETTE[colIdx++ % PALETTE.length];
        }
        groups.set(g.key, {
          key: g.key,
          label: g.label,
          category: g.category,
          color: color,
          count: 0,
        });
      }
      groups.get(g.key).count++;
    });
    return groups;
  }

  function renderComposition(groups) {
    const list = [...groups.values()].sort((a, b) => b.count - a.count);
    const total = rows.length || 1;
    const wrap = $("#composition");
    wrap.innerHTML = "";
    list.forEach((g) => {
      const item = document.createElement("div");
      item.className = "comp-item";
      item.style.setProperty("--c", g.color);
      const pct = ((g.count / total) * 100).toFixed(1);
      item.innerHTML =
        `<div class="comp-row"><span>${esc(g.label)}</span><span class="n">${g.count} reads (${pct}%)</span></div>` +
        `<div class="track"><div class="fill" style="width: 0%"></div></div>`;
      wrap.appendChild(item);
      requestAnimationFrame(() => {
        $(".fill", item).style.width = Math.max(3, (g.count / total) * 100) + "%";
      });
    });
  }

  function renderLegend(groups) {
    const legend = $("#map-legend");
    legend.innerHTML = `
      <span style="--c:${ANCHOR_COLOR}"><i class="dot star"></i>Reference Anchor Beacon</span>
      <span style="--c:#34f5c5"><i class="dot"></i>Known Taxon</span>
      <span style="--c:${SISTER}"><i class="dot"></i>Sister Lineage</span>
      <span style="--c:${NOVEL}"><i class="dot tri"></i>Candidate Novel Taxon</span>
      <span style="--c:${ARTIFACT}"><i class="dot"></i>Artifact / Outlier</span>
    `;
  }

  /* ---------------------------------------------------------------
     Sonar 2D Map with Reference Anchor Beacons
  --------------------------------------------------------------- */
  function createSonar(canvas, tip, allRows, refAnchors, groups) {
    const ctx = canvas.getContext("2d");
    const pts = allRows
      .filter((r) => r.x != null && r.y != null)
      .map((r, i) => ({ row: r, born: i * 25, g: groups.get(r._group) || { color: "#34f5c5" } }));

    // Reference anchor beacons
    const beacons = [];
    if (refAnchors) {
      Object.entries(refAnchors).forEach(([name, coords]) => {
        beacons.push({ name: name, x: coords[0], y: coords[1] });
      });
    }

    let w = 0, h = 0, dpr = 1;
    let hover = null;
    let raf = null;
    let alive = true;
    const t0 = performance.now();
    const animated = !reduceMotion;

    function layout() {
      const rect = canvas.getBoundingClientRect();
      w = Math.max(240, rect.width);
      h = Math.max(220, rect.height);
      dpr = window.devicePixelRatio || 1;
      canvas.width = Math.round(w * dpr);
      canvas.height = Math.round(h * dpr);
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

      const allX = pts.map((p) => p.row.x).concat(beacons.map((b) => b.x));
      const allY = pts.map((p) => p.row.y).concat(beacons.map((b) => b.y));
      const minX = Math.min.apply(null, allX.length ? allX : [0]);
      const maxX = Math.max.apply(null, allX.length ? allX : [1]);
      const minY = Math.min.apply(null, allY.length ? allY : [0]);
      const maxY = Math.max.apply(null, allY.length ? allY : [1]);

      const pad = 44;
      const spanX = maxX - minX || 1;
      const spanY = maxY - minY || 1;

      pts.forEach((p) => {
        p.px = pad + ((p.row.x - minX) / spanX) * (w - pad * 2);
        p.py = pad + ((p.row.y - minY) / spanY) * (h - pad * 2);
      });
      beacons.forEach((b) => {
        b.px = pad + ((b.x - minX) / spanX) * (w - pad * 2);
        b.py = pad + ((b.y - minY) / spanY) * (h - pad * 2);
      });
    }

    function drawStar(x, y, r, p = 5, m = 0.5) {
      ctx.beginPath();
      for (let i = 0; i < p * 2; i++) {
        const rad = (i * Math.PI) / p - Math.PI / 2;
        const dist = i % 2 === 0 ? r : r * m;
        const px = x + Math.cos(rad) * dist;
        const py = y + Math.sin(rad) * dist;
        if (i === 0) ctx.moveTo(px, py);
        else ctx.lineTo(px, py);
      }
      ctx.closePath();
    }

    function draw(t) {
      ctx.clearRect(0, 0, w, h);
      const bg = ctx.createRadialGradient(w / 2, h / 2, 10, w / 2, h / 2, Math.hypot(w, h) / 2);
      bg.addColorStop(0, "rgba(10, 48, 42, 0.45)");
      bg.addColorStop(1, "rgba(2, 10, 9, 0)");
      ctx.fillStyle = bg;
      ctx.fillRect(0, 0, w, h);

      // Sonar circles
      const cx = w / 2, cy = h / 2;
      const maxR = Math.hypot(cx, cy);
      ctx.strokeStyle = "rgba(52, 245, 197, 0.08)";
      ctx.lineWidth = 1;
      for (let r = maxR / 4; r <= maxR; r += maxR / 4) {
        ctx.beginPath();
        ctx.arc(cx, cy, r, 0, Math.PI * 2);
        ctx.stroke();
      }

      // Draw Beacons (Reference Anchors)
      beacons.forEach((b) => {
        ctx.fillStyle = ANCHOR_COLOR;
        ctx.shadowColor = ANCHOR_COLOR;
        ctx.shadowBlur = 14;
        drawStar(b.px, b.py, 7, 5, 0.45);
        ctx.fill();
        ctx.shadowBlur = 0;

        ctx.fillStyle = "rgba(255, 209, 102, 0.75)";
        ctx.font = "10px JetBrains Mono, monospace";
        ctx.fillText(b.name, b.px + 9, b.py + 3);
      });

      // Draw Sequence Reads
      pts.forEach((p) => {
        const color = p.g.color;
        const isActive = p.row.id === activeId;
        const isHover = hover === p;
        const isNovel = p.row.category === "novel_taxon";
        const base = (isNovel ? 7 : 5);

        ctx.fillStyle = color;
        ctx.shadowColor = color;
        ctx.shadowBlur = isActive || isHover ? 20 : 10;

        if (isNovel) {
          // Triangle
          ctx.beginPath();
          ctx.moveTo(p.px, p.py - base - 2);
          ctx.lineTo(p.px + base + 1, p.py + base);
          ctx.lineTo(p.px - base - 1, p.py + base);
          ctx.closePath();
          ctx.fill();
        } else if (p.row.category === "artifact_outlier") {
          // Diamond
          ctx.beginPath();
          ctx.moveTo(p.px, p.py - base);
          ctx.lineTo(p.px + base, p.py);
          ctx.lineTo(p.px, p.py + base);
          ctx.lineTo(p.px - base, p.py);
          ctx.closePath();
          ctx.fill();
        } else {
          // Circle
          ctx.beginPath();
          ctx.arc(p.px, p.py, base, 0, Math.PI * 2);
          ctx.fill();
        }
        ctx.shadowBlur = 0;

        if (isActive) {
          ctx.strokeStyle = "#ffffff";
          ctx.lineWidth = 2;
          ctx.beginPath();
          ctx.arc(p.px, p.py, base + 6, 0, Math.PI * 2);
          ctx.stroke();
        }
      });
    }

    function frame(t) {
      if (!alive) return;
      draw(t);
      raf = null;
      if (animated) raf = requestAnimationFrame(frame);
    }

    function requestDraw() {
      if (!raf && alive) raf = requestAnimationFrame(frame);
    }

    function nearest(mx, my) {
      let best = null;
      let bestD = 20;
      pts.forEach((p) => {
        const d = Math.hypot(p.px - mx, p.py - my);
        if (d < bestD) {
          bestD = d;
          best = p;
        }
      });
      return best;
    }

    function onMove(e) {
      const rect = canvas.getBoundingClientRect();
      const mx = e.clientX - rect.left;
      const my = e.clientY - rect.top;
      hover = nearest(mx, my);
      canvas.style.cursor = hover ? "pointer" : "default";

      if (hover) {
        const r = hover.row;
        tip.innerHTML = `
          <strong>${esc(r.id)}</strong><br>
          <span style="color:${hover.g.color}">● ${esc(r.category.toUpperCase().replace('_', ' '))}</span><br>
          ${r.closest_known_taxon ? `<em>Closest:</em> ${esc(r.closest_known_taxon)}<br>` : ''}
          <em>Divergence:</em> d=${r.nearest_dist != null ? r.nearest_dist.toFixed(2) : 'N/A'}<br>
          <em>cOTU Cluster:</em> ${r.cluster} &bull; <em>GC:</em> ${r.gc_pct}%
        `;
        tip.style.left = Math.min(mx + 10, rect.width - 230) + "px";
        tip.style.top = Math.min(my + 10, rect.height - 110) + "px";
        tip.style.opacity = 1;
      } else {
        tip.style.opacity = 0;
      }
      if (!animated) requestDraw();
    }

    function onClick(e) {
      const rect = canvas.getBoundingClientRect();
      const p = nearest(e.clientX - rect.left, e.clientY - rect.top);
      if (p) {
        setActive(p.row.id, true);
        // Switch to Read Inspector tab to highlight
        const readsTabBtn = $('[data-tab="tab-reads"]');
        if (readsTabBtn) readsTabBtn.click();
      }
    }

    function onResize() {
      layout();
      requestDraw();
    }

    canvas.addEventListener("mousemove", onMove);
    canvas.addEventListener("mouseleave", () => { hover = null; tip.style.opacity = 0; });
    canvas.addEventListener("click", onClick);
    window.addEventListener("resize", onResize);

    layout();
    requestDraw();

    return {
      requestDraw: requestDraw,
      destroy: () => {
        alive = false;
        if (raf) cancelAnimationFrame(raf);
        canvas.removeEventListener("mousemove", onMove);
        canvas.removeEventListener("click", onClick);
        window.removeEventListener("resize", onResize);
      },
    };
  }

  /* ---------------------------------------------------------------
     Canvas Chart Drawing Helpers (Pure Vanilla JS, 0 CDN Deps)
  --------------------------------------------------------------- */
  function drawBarChart(canvasId, bins, counts, label, barColor = "#34f5c5") {
    const canvas = document.getElementById(canvasId);
    if (!canvas || !counts || !counts.length) return;
    const ctx = canvas.getContext("2d");
    const rect = canvas.getBoundingClientRect();
    const W = Math.max(200, rect.width || 360);
    const H = Math.max(160, rect.height || 180);
    canvas.width = W;
    canvas.height = H;

    ctx.clearRect(0, 0, W, H);
    const maxVal = Math.max.apply(null, counts) || 1;
    const padX = 40, padY = 30;
    const chartW = W - padX - 16;
    const chartH = H - padY - 20;

    const n = counts.length;
    const barW = Math.max(8, (chartW / n) * 0.7);
    const gap = (chartW / n);

    // Axes
    ctx.strokeStyle = "rgba(120, 255, 214, 0.2)";
    ctx.beginPath();
    ctx.moveTo(padX, 10);
    ctx.lineTo(padX, H - padY);
    ctx.lineTo(W - 10, H - padY);
    ctx.stroke();

    // Bars
    counts.forEach((c, i) => {
      const bh = (c / maxVal) * chartH;
      const bx = padX + i * gap + (gap - barW) / 2;
      const by = H - padY - bh;

      ctx.fillStyle = barColor;
      ctx.fillRect(bx, by, barW, bh);

      // Bin label
      ctx.fillStyle = "#7fa59b";
      ctx.font = "9px JetBrains Mono, monospace";
      ctx.textAlign = "center";
      const lbl = bins[i] ? bins[i].split("-")[0] : "";
      ctx.fillText(lbl, bx + barW / 2, H - padY + 14);

      // Value label on top
      if (c > 0) {
        ctx.fillStyle = "#d8f3ec";
        ctx.fillText(c, bx + barW / 2, by - 4);
      }
    });
  }

  function drawWhittakerChart(whittakerData) {
    const canvas = document.getElementById("whittaker-chart");
    if (!canvas || !whittakerData || !whittakerData.length) return;
    const ctx = canvas.getContext("2d");
    const rect = canvas.getBoundingClientRect();
    const W = Math.max(200, rect.width || 420);
    const H = Math.max(180, rect.height || 220);
    canvas.width = W;
    canvas.height = H;

    ctx.clearRect(0, 0, W, H);
    const padL = 50, padB = 30, padR = 20, padT = 20;
    const cW = W - padL - padR;
    const cH = H - padB - padT;

    const maxPct = Math.max.apply(null, whittakerData.map((d) => d.relative_abundance_pct)) || 1;
    const n = whittakerData.length;

    // Grid lines
    ctx.strokeStyle = "rgba(120, 255, 214, 0.15)";
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(padL, padT);
    ctx.lineTo(padL, H - padB);
    ctx.lineTo(W - padR, H - padB);
    ctx.stroke();

    // Whittaker Curve Points & Line
    ctx.strokeStyle = "#34f5c5";
    ctx.lineWidth = 2;
    ctx.beginPath();

    const pts = [];
    whittakerData.forEach((d, i) => {
      const px = padL + (i / Math.max(1, n - 1)) * cW;
      const py = H - padB - (d.relative_abundance_pct / maxPct) * cH;
      pts.push([px, py]);
      if (i === 0) ctx.moveTo(px, py);
      else ctx.lineTo(px, py);
    });
    ctx.stroke();

    // Dots
    pts.forEach(([px, py], i) => {
      ctx.fillStyle = "#ffd166";
      ctx.beginPath();
      ctx.arc(px, py, 4, 0, Math.PI * 2);
      ctx.fill();

      // Axis labels
      ctx.fillStyle = "#7fa59b";
      ctx.font = "9px JetBrains Mono, monospace";
      ctx.textAlign = "center";
      ctx.fillText(`R${i + 1}`, px, H - padB + 14);
    });

    // Y Axis labels
    ctx.textAlign = "right";
    ctx.fillStyle = "#7fa59b";
    ctx.fillText(`${maxPct.toFixed(0)}%`, padL - 6, padT + 8);
    ctx.fillText(`0%`, padL - 6, H - padB);
  }

  /* ---------------------------------------------------------------
     QC Summary Cards & Charts
  --------------------------------------------------------------- */
  function renderQCSummary(qc) {
    if (!qc) return;
    const wrap = $("#qc-stat-cards");
    const stats = [
      { label: "Total Reads", val: qc.total_reads, sub: "Metabarcoding sequences" },
      { label: "Total Bases", val: `${(qc.total_bases / 1000).toFixed(1)} kbp`, sub: "Sequenced yield" },
      { label: "Mean Length", val: `${qc.mean_length} bp`, sub: `Range: ${qc.min_length}-${qc.max_length} bp` },
      { label: "Mean GC%", val: `${qc.mean_gc}%`, sub: `&plusmn; ${qc.std_gc}% std dev` },
      {
        label: "Mean Phred Score",
        val: qc.quality_summary ? `Q${qc.quality_summary.mean_phred_score}` : "N/A (FASTA)",
        sub: qc.quality_summary ? `Q30 rate: ${qc.quality_summary.q30_percentage}%` : "No quality scores",
      },
      { label: "Homopolymer Flags", val: qc.homopolymer_flagged_reads, sub: ">= 8 consecutive bases" },
    ];

    wrap.innerHTML = stats
      .map(
        (s) => `
        <div class="metric">
          <div class="m-label">${esc(s.label)}</div>
          <div class="m-value" style="font-size: 26px;">${s.val}</div>
          <div class="m-cap">${s.sub}</div>
        </div>
      `
      )
      .join("");

    drawQCCharts(qc);
  }

  function drawQCCharts(qc) {
    if (!qc) return;
    drawBarChart("gc-chart", qc.gc_distribution.bins, qc.gc_distribution.counts, "GC%", "#34f5c5");
    drawBarChart("length-chart", qc.length_distribution.bins, qc.length_distribution.counts, "Length", "#5ab8ff");

    if (qc.quality_summary) {
      $("#quality-score-panel").classList.remove("hidden");
      $("#qc-q30-label").textContent = `Q30 Pass Rate: ${qc.quality_summary.q30_percentage}%`;
      drawBarChart(
        "qscore-chart",
        qc.quality_summary.distribution.bins,
        qc.quality_summary.distribution.counts,
        "Q Score",
        "#ffd166"
      );
    } else {
      $("#quality-score-panel").classList.add("hidden");
    }
  }

  /* ---------------------------------------------------------------
     Cluster Summaries Table
  --------------------------------------------------------------- */
  function renderClustersTable(summaries) {
    const tbody = $("#clusters-table tbody");
    if (!tbody) return;
    tbody.innerHTML = "";
    (summaries || []).forEach((c) => {
      const tr = document.createElement("tr");
      let catBadge = '<span class="pill pill-known">Known</span>';
      if (c.dominant_category === "sister_taxon") catBadge = '<span class="pill pill-sister">Sister</span>';
      else if (c.dominant_category === "novel_taxon") catBadge = '<span class="pill pill-novel">Novel Clade</span>';
      else if (c.dominant_category === "artifact_outlier") catBadge = '<span class="pill pill-artifact">Outlier</span>';

      tr.innerHTML = `
        <td><strong>cOTU_${c.cluster_id}</strong></td>
        <td>${esc(c.dominant_taxon)}</td>
        <td>${catBadge}</td>
        <td>${c.size}</td>
        <td>${c.dispersion.toFixed(3)}</td>
      `;
      tbody.appendChild(tr);
    });
  }

  /* ---------------------------------------------------------------
     Benchmark: Alignment vs BioEmbed
  --------------------------------------------------------------- */
  function renderBenchmark(bench, bioIndices) {
    if (!bench) return;
    $("#blast-discard-rate").textContent = `${bench.blast_discarded_pct}%`;
    $("#blast-richness-val").textContent = `Observed Richness: ${bench.blast_richness}`;
    $("#blast-shannon-val").textContent = `Shannon Diversity: ${bench.blast_shannon}`;

    $("#bioembed-rescue-rate").textContent = `${bench.rescued_dark_reads_count} reads (${bench.rescued_dark_reads_pct}%)`;
    $("#bioembed-richness-val").textContent = `True Total Richness: ${bioIndices.richness}`;
    $("#bioembed-shannon-val").textContent = `True Shannon Diversity: ${bioIndices.shannon_diversity}`;

    const tbody = $("#benchmark-table tbody");
    tbody.innerHTML = "";
    (bench.per_read_comparison || []).forEach((r) => {
      const tr = document.createElement("tr");
      let blastPill = '<span class="pill pill-novel">Discarded (&lt;90%)</span>';
      if (r.blast_status === "assigned_species") blastPill = '<span class="pill pill-known">Assigned (Species &ge;97%)</span>';
      else if (r.blast_status === "assigned_genus") blastPill = '<span class="pill pill-sister">Assigned (Genus)</span>';

      let bioOutcome = r.rescued
        ? '<strong style="color:var(--glow)">⚡ Rescued as Novel Taxon</strong>'
        : `Classified as ${r.bioembed_category.replace('_', ' ')}`;

      tr.innerHTML = `
        <td>Read_${r.index + 1}</td>
        <td>${esc(r.blast_call)}</td>
        <td>${r.blast_identity}%</td>
        <td>${blastPill}</td>
        <td><span class="pill pill-${r.bioembed_category.startsWith('novel') ? 'novel' : 'known'}">${esc(r.bioembed_category)}</span></td>
        <td>${esc(r.bioembed_nearest || 'None')}</td>
        <td>${bioOutcome}</td>
      `;
      tbody.appendChild(tr);
    });
  }

  /* ---------------------------------------------------------------
     Beta Diversity Comparison Runner
  --------------------------------------------------------------- */
  $("#btn-run-compare")?.addEventListener("click", async () => {
    const sA = $("#select-sample-a").value;
    const sB = $("#select-sample-b").value;

    const fd = new FormData();
    fd.append("preset_a", sA);
    fd.append("preset_b", sB);

    try {
      const res = await fetch("/api/compare_samples", { method: "POST", body: fd });
      const data = await res.json();
      if (!res.ok || data.error) {
        alert(data.error || "Comparison failed");
        return;
      }

      $("#beta-results").classList.remove("hidden");
      $("#bc-val").textContent = data.beta_diversity.bray_curtis_dissimilarity.toFixed(3);
      $("#jaccard-val").textContent = data.beta_diversity.jaccard_distance.toFixed(3);

      // Overlap bar chart
      const overlapCanvas = document.getElementById("beta-overlap-chart");
      if (overlapCanvas) {
        const bins = ["Shared Taxa", `Unique to ${data.sample_a.name}`, `Unique to ${data.sample_b.name}`];
        const counts = [
          data.beta_diversity.shared_taxa_count,
          data.beta_diversity.unique_to_a_count,
          data.beta_diversity.unique_to_b_count,
        ];
        drawBarChart("beta-overlap-chart", bins, counts, "Taxa Partition", "#5ab8ff");
      }
    } catch (e) {
      alert("Could not compute beta diversity comparison: " + e.message);
    }
  });

  /* ---------------------------------------------------------------
     Read Inspector Table
  --------------------------------------------------------------- */
  function renderFilters() {
    const known = rows.filter((r) => r.category === "known").length;
    const sister = rows.filter((r) => r.category === "sister_taxon").length;
    const novel = rows.filter((r) => r.category === "novel_taxon").length;
    const artifact = rows.filter((r) => r.category === "artifact_outlier").length;

    const opts = [
      ["all", `All (${rows.length})`],
      ["known", `Known (${known})`],
      ["sister_taxon", `Sister Lineages (${sister})`],
      ["novel_taxon", `Candidate Novel (${novel})`],
      ["artifact_outlier", `Artifacts (${artifact})`],
    ];

    const wrap = $("#filters");
    wrap.innerHTML = "";
    opts.forEach(([mode, label]) => {
      const b = document.createElement("button");
      b.type = "button";
      b.className = "chip" + (mode === filterMode ? " on" : "");
      b.textContent = label;
      b.addEventListener("click", () => {
        filterMode = mode;
        renderFilters();
        renderTable();
      });
      wrap.appendChild(b);
    });
  }

  $("#search-reads")?.addEventListener("input", () => renderTable());

  function renderTable() {
    const tbody = $("#results-table tbody");
    tbody.innerHTML = "";
    const query = ($("#search-reads")?.value || "").toLowerCase().trim();

    rows.forEach((r) => {
      if (filterMode !== "all" && r.category !== filterMode) return;
      if (query && !r.id.toLowerCase().includes(query) && !(r.closest_known_taxon || "").toLowerCase().includes(query)) {
        return;
      }

      const tr = document.createElement("tr");
      tr.dataset.id = r.id;
      if (r.id === activeId) tr.className = "row-active";

      let pillClass = "pill-known";
      if (r.category === "sister_taxon") pillClass = "pill-sister";
      else if (r.category === "novel_taxon") pillClass = "pill-novel";
      else if (r.category === "artifact_outlier") pillClass = "pill-artifact";

      tr.innerHTML = `
        <td><strong>${esc(r.id)}</strong></td>
        <td>cOTU_${r.cluster}</td>
        <td><span class="pill ${pillClass}">${esc(r.category.replace('_', ' '))}</span></td>
        <td>${esc(r.closest_known_taxon || 'None')}</td>
        <td style="font-size: 11px; color: var(--ink-dim);">${esc(r.taxonomic_lineage || 'Unspecified')}</td>
        <td>${r.nearest_dist != null ? r.nearest_dist.toFixed(3) : '&mdash;'}</td>
        <td>${r.confidence != null ? (r.confidence * 100).toFixed(0) + '%' : '&mdash;'}</td>
        <td>${r.length} bp</td>
        <td>${r.gc_pct}%</td>
      `;
      tr.addEventListener("click", () => setActive(r.id, false));
      tbody.appendChild(tr);
    });
  }

  function setActive(id, fromMap) {
    activeId = activeId === id ? null : id;
    renderTable();
    if (mapSonar) mapSonar.requestDraw();
    if (fromMap && activeId) {
      const tr = document.querySelector(`#results-table tr[data-id="${CSS.escape(activeId)}"]`);
      if (tr && tr.scrollIntoView) tr.scrollIntoView({ behavior: "smooth", block: "center" });
    }
  }

  /* ---------------------------------------------------------------
     Thesis Dossier Summary
  --------------------------------------------------------------- */
  function renderThesisDossier(data) {
    const box = $("#thesis-run-stats");
    if (!box) return;
    box.innerHTML = `
      <strong>Run Parameter Verification:</strong><br>
      &bull; Total Processed Sequences: ${data.total_sequences}<br>
      &bull; Active Embedding Feature Space: Canonical 4-mer (136 dimensions, L2-normalized)<br>
      &bull; Loaded Reference Centroid Taxa: ${data.reference_taxa_loaded}<br>
      &bull; Clustered cOTUs: ${data.biodiversity_indices.richness}<br>
      &bull; Chao1 Richness Estimator: ${data.biodiversity_indices.chao1_richness}<br>
      &bull; Shannon Entropy (H'): ${data.biodiversity_indices.shannon_diversity}<br>
      &bull; Pielou Evenness (J'): ${data.biodiversity_indices.pielou_evenness}<br>
      &bull; Dark Taxa Rescued from Alignment Discard: ${data.benchmark.rescued_dark_reads_count} reads (${data.benchmark.rescued_dark_reads_pct}%)<br>
      &bull; Pipeline Timestamp: ${new Date().toLocaleString()}
    `;
  }

  /* ---------------------------------------------------------------
     CSV Export & Print
  --------------------------------------------------------------- */
  $("#btn-export-csv")?.addEventListener("click", async () => {
    if (!currentDataset || !currentDataset.sequence_results) return;
    try {
      const res = await fetch("/api/export_csv", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ sequence_results: currentDataset.sequence_results }),
      });
      const blob = await res.blob();
      const url = window.URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = "bioembed_edna_taxonomy_results.csv";
      document.body.appendChild(a);
      a.click();
      a.remove();
    } catch (e) {
      alert("CSV export failed: " + e.message);
    }
  });

  $("#btn-print-report")?.addEventListener("click", () => {
    // Switch to thesis tab first
    const thesisBtn = $('[data-tab="tab-thesis"]');
    if (thesisBtn) thesisBtn.click();
    setTimeout(() => window.print(), 300);
  });

  // Initialize
  buildTicker();
  startHelix();

  // Auto-load preset demo on first visit for immediate visual demonstration!
  setTimeout(() => {
    runAnalysis({ preset_id: "marine_coastal_edna" }, "🌊 Marine Coastal eDNA (Fish + Dark Protists)");
  }, 400);
})();