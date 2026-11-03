// Detailed CO2/C1 reaction network (no carbonylation, no NH3/urea): overview + six expanded modules
const pptxgen = require("pptxgenjs");
const { applyTheme } = require("/root/.claude/skills/synced/62cc770a-8c2c-45f8-8ae0-ca258c8631d4_456331ad-9fc1-4922-9170-b4c06151a9b5/pptx/scripts/apply_theme.js");
const OUT = process.argv[2] || "detailed.pptx";
const ONLY = process.argv[3] ? process.argv[3].split(",") : null;   // debug: render subset of slides

const THEME = {
  name: "C1 Network Detailed", headFontFace: "Microsoft YaHei", bodyFontFace: "Microsoft YaHei",
  colors: { dk1: "14232B", lt1: "FFFFFF", dk2: "3B4A54", lt2: "F3F5F6", accent1: "E07B24", accent2: "D63B3B",
            accent3: "2A78D6", accent4: "1E9E5A", accent5: "8B5A2B", accent6: "C2185B", hlink: "2A78D6", folHlink: "8B5A2B" },
};
const pres = new pptxgen();
const SW = 24, SH = 13.5;
pres.defineLayout({ name: "POSTER", width: SW, height: SH });
pres.layout = "POSTER";
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
pres.title = "CO2 加氢及 C1 催化细致反应网络（含自由基）";
const C = pres.SchemeColor;
const CAT = {
  R1: { c: C.accent1, n: "① 重整 / 合成气" }, R2: { c: C.accent2, n: "② 费托 / C–C 链增长" }, R3: { c: C.accent3, n: "③ 加氢还原" },
  R4: { c: C.accent4, n: "④ 酸催化 / 分子筛" }, R5: { c: C.accent5, n: "⑤ 烃类升级" }, R6: { c: C.accent6, n: "⑥ 自由基链反应", dash: "dash" },
  SH: { c: C.text2, n: "多类反应共享步骤" },
};
const CAT_ORDER = ["R1", "R2", "R3", "R4", "R5", "R6"];

pres.defineSlideMaster({ title: "POSTER", background: { color: C.background1 },
  objects: [{ placeholder: { options: { name: "title", type: "title", x: 0.3, y: 0.12, w: 23.4, h: 0.55, fontSize: 24, bold: true,
    color: C.text1, align: "left", valign: "middle", margin: 0 }, text: "" } }],
  slideNumber: { x: 23.1, y: 13.12, w: 0.6, h: 0.3, fontSize: 10, color: C.text2, align: "right" } });

// ------------------------------------------------------------------ text helpers
function runs(s, opt = {}) {
  const out = [], re = /_\{([^}]*)\}|_(\d+|[a-z])|\^\{([^}]*)\}/g;
  let last = 0, m;
  while ((m = re.exec(s))) {
    if (m.index > last) out.push({ text: s.slice(last, m.index), options: { ...opt } });
    if (m[3] !== undefined) out.push({ text: m[3], options: { ...opt, superscript: true } });
    else out.push({ text: m[1] ?? m[2], options: { ...opt, subscript: true } });
    last = re.lastIndex;
  }
  if (last < s.length) out.push({ text: s.slice(last), options: { ...opt } });
  return out;
}
function textW(s, pt) {   // rough rendered width in inches
  let em = 0;
  const plain = s.replace(/_\{([^}]*)\}|_(\d+|[a-z])|\^\{([^}]*)\}/g, (m, a, b, c) => "\u0001".repeat((a ?? b ?? c).length));
  for (const ch of plain) {
    const c = ch.codePointAt(0);
    if (ch === "\u0001") em += 0.46; else if (c > 0x2e80) em += 1.02; else if (/[A-Z0-9]/.test(ch)) em += 0.7;
    else if (/[a-z]/.test(ch)) em += 0.58; else if (ch === " ") em += 0.3; else em += 0.55;
  }
  return em * pt / 72;
}

// ------------------------------------------------------------------ generic drawing with a transform T = {ox, oy, s}
const TX = (T, x) => T.ox + T.s * x, TY = (T, y) => T.oy + T.s * y;
function seg(slide, T, x1, y1, x2, y2, o) {
  const X1 = TX(T, x1), Y1 = TY(T, y1), X2 = TX(T, x2), Y2 = TY(T, y2);
  slide.addShape(pres.shapes.LINE, {
    x: Math.min(X1, X2), y: Math.min(Y1, Y2), w: Math.max(Math.abs(X2 - X1), 0.001), h: Math.max(Math.abs(Y2 - Y1), 0.001),
    flipH: X2 < X1, flipV: Y2 < Y1,
    line: { color: o.color, width: Math.max(0.6, (o.width || 1.1) * T.s), dashType: o.dash || "solid",
            endArrowType: o.end ? "triangle" : undefined, beginArrowType: o.begin ? "triangle" : undefined },
    objectName: o.name,
  });
}
function txt(slide, T, x, y, w, h, text, o = {}) {
  slide.addText(Array.isArray(text) ? text : runs(text, { fontSize: (o.fs || 6) * T.s, color: o.color || C.text2, bold: !!o.bold, italic: !!o.italic }), {
    x: TX(T, x), y: TY(T, y), w: w * T.s, h: h * T.s, align: o.align || "left", valign: o.valign || "middle", margin: 0,
    isTextBox: true, objectName: o.name || "text", ...(o.fill ? { fill: { color: o.fill } } : {}),
  });
}
const NFS = 7, LFS = 5.5;   // local node / label font sizes (pt at scale 1)
function nodeSize(nd) {
  const fs = nd.fs || NFS;
  return [nd.w || Math.max(0.42, textW(nd.label, fs) + (nd.type === "rad" ? 0.3 : 0.2)), nd.h || 0.25];
}
function clip(nd, ux, uy) {
  const [w, h] = nodeSize(nd);
  if (nd.type === "rad") { const a = w / 2, b = h / 2, t = 1 / Math.sqrt((ux * ux) / (a * a) + (uy * uy) / (b * b)); return [nd.x + ux * t, nd.y + uy * t]; }
  const tx = ux ? w / 2 / Math.abs(ux) : Infinity, ty = uy ? h / 2 / Math.abs(uy) : Infinity, t = Math.min(tx, ty);
  return [nd.x + ux * t, nd.y + uy * t];
}
function drawEdge(slide, T, Nmap, ed, i, mod) {
  const A = Nmap[ed.a], B = Nmap[ed.b], cat = CAT[ed.cat];
  if (!A || !B) throw new Error(`edge ${mod}:${ed.a}->${ed.b} missing node`);
  const st = { color: cat.c, width: ed.width || 1.1, dash: ed.dash || cat.dash };
  const nm = `${mod} edge ${i + 1} ${ed.a}→${ed.b}`;
  const lab = (x, y) => {
    if (!ed.label) return;
    const fs = ed.fs || LFS, lw = ed.lw || textW(ed.label, fs) + 0.12, lh = fs / 72 * 1.5;
    txt(slide, T, x - lw / 2, y - lh / 2, lw, lh, runs(ed.label, { fontSize: fs * T.s, color: cat.c }), { align: "center", fill: C.background1, name: nm + " label" });
  };
  if (ed.pts) {
    const P = ed.pts;
    for (let k = 0; k < P.length - 1; k++)
      seg(slide, T, P[k][0], P[k][1], P[k + 1][0], P[k + 1][1], { ...st, end: k === P.length - 2, begin: ed.both && k === 0, name: `${nm} seg ${k + 1}` });
    if (ed.lpos) lab(ed.lpos[0], ed.lpos[1]);
    return;
  }
  const ax = A.x + (ed.ax || 0), ay = A.y + (ed.ay || 0), bx = B.x + (ed.bx || 0), by = B.y + (ed.by || 0);
  const dx = bx - ax, dy = by - ay, L = Math.hypot(dx, dy), ux = dx / L, uy = dy / L;
  let [x1, y1] = clip({ ...A, x: ax, y: ay }, ux, uy);
  let [x2, y2] = clip({ ...B, x: bx, y: by }, -ux, -uy);
  const g = 0.02; x1 += ux * g; y1 += uy * g; x2 -= ux * g; y2 -= uy * g;
  seg(slide, T, x1, y1, x2, y2, { ...st, end: true, begin: !!ed.both, name: nm });
  let nx = -uy, ny = ux;
  if (ny > 0 || (Math.abs(ny) < 1e-9 && nx > 0)) { nx = -nx; ny = -ny; }
  const t = ed.lt ?? 0.5, pp = ed.perp ?? 0;
  lab(x1 + (x2 - x1) * t + nx * pp + (ed.lx || 0), y1 + (y2 - y1) * t + ny * pp + (ed.ly || 0));
}
function drawNode(slide, T, nd, mod, cats) {
  const [w, h] = nodeSize(nd);
  const x = TX(T, nd.x - w / 2), y = TY(T, nd.y - h / 2), W = w * T.s, H = h * T.s;
  const fs = (nd.fs || NFS) * T.s;
  if (nd.type === "txt") {
    slide.addText(runs(nd.label, { fontSize: fs, color: nd.tc || C.text2, bold: !!nd.bold, italic: nd.italic !== false }),
      { x, y, w: W, h: H, align: nd.align || "center", valign: "middle", margin: 0, isTextBox: true, objectName: `${mod} node ${nd.id}` });
    return;
  }
  const S = {
    mol: { shape: pres.shapes.ROUNDED_RECTANGLE, fill: nd.fill || C.text2, line: nd.fill || C.text2, tc: C.background1, r: 0.05 },
    surf: { shape: pres.shapes.ROUNDED_RECTANGLE, fill: C.background1, line: C.text1, tc: C.text1, r: 0.5 },
    rad: { shape: pres.shapes.OVAL, fill: C.background1, line: C.accent6, tc: C.text1, dash: "dash" },
    ion: { shape: pres.shapes.RECTANGLE, fill: C.background1, line: C.text1, tc: C.text1 },
    hub: { shape: pres.shapes.ROUNDED_RECTANGLE, fill: C.accent1, line: C.accent1, tc: C.background1, r: 0.08 },
  }[nd.type];
  if (nd.type === "ion") {
    const d = 0.035 * T.s;
    slide.addShape(pres.shapes.RECTANGLE, { x: x - d, y: y - d, w: W + 2 * d, h: H + 2 * d, fill: { color: C.background1 },
      line: { color: C.text1, width: Math.max(0.5, 0.75 * T.s) }, objectName: `${mod} node ${nd.id} outer` });
  }
  const opt = { x, y, w: W, h: H, shape: S.shape, fill: { color: S.fill },
    line: { color: S.line, width: Math.max(0.5, (nd.type === "rad" ? 1.1 : 0.9) * T.s), dashType: S.dash || "solid" },
    align: "center", valign: "middle", margin: 0, objectName: `${mod} node ${nd.id}` };
  if (S.r !== undefined) opt.rectRadius = Math.min(S.r, 0.5) * Math.min(W, H) / Math.min(W, H) * (S.r === 0.5 ? H / 2 : S.r * T.s) / (S.r === 0.5 ? H / 2 : 1);
  if (S.r === 0.5) opt.rectRadius = H / 2; else if (S.r !== undefined) opt.rectRadius = S.r * T.s;
  slide.addText(runs(nd.label, { fontSize: fs, color: S.tc, bold: true }), opt);
  if (cats && cats.length) {
    const d = 0.075 * T.s, gap = 0.03 * T.s, tot = cats.length * d + (cats.length - 1) * gap;
    cats.forEach((k, i) => slide.addShape(pres.shapes.OVAL, { x: x + W - tot + i * (d + gap), y: y + H + (nd.type === "ion" ? 0.07 : 0.035) * T.s,
      w: d, h: d, fill: { color: CAT[k].c }, line: { color: CAT[k].c, width: 0.5 }, objectName: `${mod} node ${nd.id} dot ${k}` }));
  }
}
function drawZone(slide, T, z) {
  slide.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: TX(T, z.x), y: TY(T, z.y), w: z.w * T.s, h: z.h * T.s, rectRadius: 0.1 * T.s,
    fill: { color: z.fill || C.background2 }, line: { color: z.color, width: Math.max(0.75, 1.5 * T.s), dashType: z.dash || "solid" },
    objectName: `zone ${z.title}` });
  if (z.title) txt(slide, T, z.tx ?? z.x + 0.12, z.ty ?? z.y + 0.05, z.tw ?? z.w - 0.24, 0.28, runs(z.title, { fontSize: (z.fs || 9.5) * T.s, color: z.color, bold: true }),
    { align: z.talign || "left", name: `zone title ${z.title}` });
}
// a module = { id, title, color, W, H, nodes: [[id,type,x,y,label,opts]], edges: [[a,b,cat,label,opts]], texts: [...], zones: [...], notes: [...] }
function buildMap(M) {
  const N = {};
  for (const [id, type, x, y, label, o = {}] of M.nodes) N[id] = { id, type, x, y, label, ...o };
  const E = M.edges.map(([a, b, cat, label, o = {}]) => ({ a, b, cat, label, ...o }));
  return [N, E];
}
function drawModule(slide, M, T, o = {}) {
  const [N, E] = buildMap(M);
  if (o.frame !== false) {
    slide.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: TX(T, 0), y: TY(T, 0), w: M.W * T.s, h: M.H * T.s, rectRadius: 0.12 * T.s,
      fill: { color: M.fill || C.background2 }, line: { color: M.color, width: Math.max(0.75, 1.4 * T.s), dashType: M.dash || "solid" },
      objectName: `module ${M.id} frame` });
    txt(slide, T, 0.12, 0.05, M.W - 0.24, 0.3, runs(M.title, { fontSize: 8.5 * T.s, color: M.color, bold: true }), { name: `module ${M.id} title` });
  }
  drawBlocks(slide, T, M.blocks, M.id);
  (M.zones || []).forEach((z) => drawZone(slide, T, z));
  (M.texts || []).forEach(([x, y, w, h, t, oo = {}]) => {
    if (Array.isArray(t)) {
      const fs = (oo.fs || 5.5) * T.s, col = oo.color || C.text2;
      const r = [];
      t.forEach((line, k) => { const rr = runs(line, { fontSize: fs, color: col }); if (k < t.length - 1) rr[rr.length - 1].options.breakLine = true; r.push(...rr); });
      txt(slide, T, x, y, w, h, r, { valign: oo.valign || "top" });
    } else txt(slide, T, x, y, w, h, t, { fs: 5.5, ...oo });
  });
  E.forEach((ed, i) => drawEdge(slide, T, N, ed, i, M.id));
  Object.values(N).forEach((nd) => drawNode(slide, T, nd, M.id, M.cats ? M.cats[nd.id] : null));
  drawChips(slide, T, M.blocks, M.id);
}
function legend(slide, y, x0 = 0.3) {
  const types = [["mol", "分子 / 产物"], ["surf", "表面吸附物种 *"], ["rad", "自由基 •"], ["ion", "离子 / 碳正离子 ⁺"]];
  const T = { ox: 0, oy: 0, s: 1 };
  let x = x0;
  types.forEach(([t, lab]) => {
    drawNode(slide, T, { id: `lg-${t}`, type: t, x: x + 0.25, y: y + 0.15, label: "", w: 0.45, h: 0.2 }, "legend");
    txt(slide, T, x + 0.55, y + 0.02, 1.6, 0.26, lab, { fs: 10, color: C.text1, name: "legend text" });
    x += 1.9;
  });
  x += 0.2;
  [...CAT_ORDER, "SH"].forEach((k) => {
    seg(slide, T, x, y + 0.15, x + 0.4, y + 0.15, { color: CAT[k].c, width: 1.75, dash: CAT[k].dash, end: true, name: `legend ${k}` });
    txt(slide, T, x + 0.47, y + 0.02, 1.75, 0.26, CAT[k].n, { fs: 10, color: C.text1, name: "legend text" });
    x += 2.05;
  });
}

// ================================================================== CATALYST FORMS (background blocks)
const CATAL = {
  OX:   { fill: "DCE7F6", c: "2A5DA8", s: "还原性氧化物", n: "还原性氧化物 · 氧空位（In₂O₃ / ZnZrOₓ / GaZrOₓ）" },
  CU:   { fill: "D5EEEE", c: "1C7272", s: "Cu 基 / 界面", n: "Cu 基 / Cu–氧化物界面（Cu/ZnO）" },
  MET:  { fill: "FBE7D2", c: "A85810", s: "过渡金属", n: "过渡金属 Ni / Co / Ru / Fe（CO 解离 · CHₓ · 甲烷化 · 重整）" },
  CARB: { fill: "F6DCDC", c: "A12F2F", s: "金属碳化物", n: "金属碳化物 χ-Fe₅C₂ / Co₂C（+K / Na 助剂）" },
  DUAL: { fill: "F2DDEC", c: "86336F", s: "CO 插入双位", n: "CO 插入双位（Co⁰–Co₂C / CuCo / Cu–Fe）" },
  SUP:  { fill: "E9E6DC", c: "655D47", s: "载体 / 界面", n: "载体 CeO₂ / ZrO₂" },
  BA:   { fill: "DDEFE1", c: "1D7342", s: "B 酸分子筛", n: "Brønsted 酸分子筛（H-ZSM-5 / SAPO-34 / H-β）" },
  LA:   { fill: "E9F0CC", c: "5A7410", s: "Lewis 酸改性分子筛", n: "Lewis 酸 Ga / Zn 改性分子筛（脱氢）" },
  BIF:  { fill: "F2E6D6", c: "7E5226", s: "双功能 金属–酸", n: "双功能 金属–酸（Pt / Ni + 分子筛）" },
  EC:   { fill: "E6E0F5", c: "5444A0", s: "电催化", n: "电催化阴极（Cu / Ag / Sn / Bi）" },
  GAS:  { fill: "FBE4EC", c: "B0175A", s: "无催化剂（气相）", n: "无需催化剂：气相 / 均相（高温 · 等离子体 · 光）" },
};
const CHIP_FS = 5.2;
// block = [key, x, y, w, h, chip:[x,y] | null, customLabel?]
function drawBlocks(slide, T, blocks, mod) {
  (blocks || []).forEach(([k, x, y, w, h], i) => {
    const K = CATAL[k];
    slide.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: TX(T, x), y: TY(T, y), w: w * T.s, h: h * T.s, rectRadius: 0.08 * T.s,
      fill: { color: K.fill }, line: { color: K.fill, width: 0.5 }, objectName: `${mod} catalyst block ${i + 1} ${k}` });
  });
}
function drawChips(slide, T, blocks, mod, fs = CHIP_FS) {
  (blocks || []).forEach(([k, , , , , chip, label], i) => {
    if (!chip) return;
    const K = CATAL[k], t = label || K.n, w = textW(t, fs) * 1.1 + 0.16, h = fs / 72 * 1.6;
    slide.addText(runs(t, { fontSize: fs * T.s, color: "FFFFFF", bold: true }), { x: TX(T, chip[0]), y: TY(T, chip[1]), w: w * T.s, h: h * T.s,
      shape: pres.shapes.ROUNDED_RECTANGLE, rectRadius: h * T.s / 2, fill: { color: K.c }, line: { color: K.c, width: 0.5 },
      align: "center", valign: "middle", margin: 0, objectName: `${mod} catalyst label ${i + 1} ${k}` });
  });
}
function catalystLegend(slide, y, x0 = 0.3) {
  const T = { ox: 0, oy: 0, s: 1 };
  txt(slide, T, x0, y, 1.25, 0.26, "催化剂形式（色块）", { fs: 10, color: C.text1, bold: true, name: "catalyst legend title" });
  let x = x0 + 1.4;
  for (const [k, K] of Object.entries(CATAL)) {
    slide.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y: y + 0.03, w: 0.42, h: 0.2, rectRadius: 0.05, fill: { color: K.fill },
      line: { color: K.c, width: 1 }, objectName: `catalyst legend ${k}` });
    const w = textW(K.s, 10) + 0.12;
    txt(slide, T, x + 0.48, y, w, 0.26, K.s, { fs: 10, color: C.text1, name: "catalyst legend text" });
    x += 0.48 + w + 0.22;
  }
}

// ================================================================== MODULES (local inches at scale 1)
const MODS = [];
// ---------------------------------------------------------------- A: oxide, methanol, RWGS, vacancy cycle, electro
MODS.push({ id: "A", color: C.accent3, W: 7.0, H: 6.05,
  title: "A  氧化物上 CO_2 活化与甲醇合成（In_2O_3 / ZnZrO_x / GaZrO_x，氧空位循环）· RWGS · 电化学",
  nodes: [
    ["H2", "mol", 0.5, 0.75, "H_2(g)"], ["HH", "surf", 1.95, 0.75, "Ga/In–H^{−} + O–H^{+}"], ["VO", "surf", 3.6, 0.75, "V_{O} 氧空位"],
    ["OL", "surf", 5.05, 0.75, "O_{lat} 晶格氧"], ["H2O", "mol", 6.45, 0.75, "H_2O(g)"],
    ["CO2", "mol", 0.45, 2.0, "CO_2(g)"], ["CO2s", "surf", 1.55, 2.0, "CO_2^{δ−}*@V_{O}"], ["HCOO", "surf", 2.85, 2.0, "HCOO*（双齿）"],
    ["H2COO", "surf", 4.15, 2.0, "H_2COO*"], ["H2COOH", "surf", 5.45, 2.0, "H_2COOH*"], ["OH", "surf", 6.5, 2.0, "OH*"],
    ["H2CO", "surf", 6.35, 3.25, "H_2CO*"], ["CH3O", "surf", 5.25, 3.25, "CH_3O*"], ["CH3OHs", "surf", 4.0, 3.25, "CH_3OH*"],
    ["MeOH", "mol", 2.65, 3.25, "CH_3OH(g) → E"],
    ["COOH", "surf", 0.6, 3.35, "COOH*"], ["CO", "surf", 1.55, 4.5, "CO*"], ["COg", "mol", 0.55, 5.45, "CO(g) → B"],
    ["HCO", "surf", 3.4, 4.5, "HCO*"],
    ["E0", "txt", 3.05, 5.45, "CO_2 + e^{−}", { tc: C.accent3, bold: true }], ["ECO2", "rad", 4.25, 5.45, "CO_2•^{−}"],
    ["EHCOOH", "mol", 5.85, 5.1, "HCOO^{−} / HCOOH"], ["ECO", "mol", 5.85, 5.8, "CO + OH^{−}"],
  ],
  edges: [
    ["H2", "HH", "R3", "异裂解离（M–O 对）", { ly: -0.18 }], ["OL", "VO", "R3", "+2H → H_2O↑（还原生成 V_{O}）", { ly: -0.17 }],
    ["OH", "H2O", "R3", "+H*", { perp: 0.15 }],
    ["VO", "CO2s", "SH", "CO_2 吸附于 V_{O}", { perp: 0.16 }], ["CO2", "CO2s", "SH", ""],
    ["CO2s", "HCOO", "R3", "+H^{−}", { both: true, ly: -0.15 }], ["HCOO", "H2COO", "R3", "+H^{−}", { ly: -0.15 }],
    ["H2COO", "H2COOH", "R3", "+H^{+}", { ly: -0.15 }], ["H2COOH", "H2CO", "R3", "C–O 断裂", { perp: -0.2 }],
    ["H2COOH", "OH", "R3", "", { dash: "sysDot" }],
    ["H2CO", "CH3O", "R3", "+H^{−}", { ly: -0.15 }], ["CH3O", "CH3OHs", "R3", "+H^{+}", { ly: -0.15 }], ["CH3OHs", "MeOH", "SH", "脱附", { ly: -0.15 }],
    ["CH3OHs", "OL", "SH", "O 回填晶格（闭合循环）", { dash: "dash", perp: -0.32, lt: 0.78 }],
    ["CO2s", "COOH", "R3", "+H*", { perp: 0.15 }], ["COOH", "CO", "SH", "−OH*（缔合 RWGS）", { perp: 0.22, lt: 0.55 }],
    ["CO2s", "CO", "SH", "氧化还原 RWGS：O → V_{O}", { dash: "dash", perp: -0.62, lt: 0.82 }],
    ["CO", "COg", "SH", "脱附", { perp: 0.15 }], ["CO", "HCO", "R3", "+H*（次要）", { ly: -0.15 }], ["HCO", "H2CO", "R3", "+H*", { perp: 0.16 }],
    ["E0", "ECO2", "R3", ""], ["ECO2", "EHCOOH", "R3", "+H^{+} + e^{−}", { perp: 0.15 }], ["ECO2", "ECO", "R3", "+H^{+} + e^{−}", { perp: -0.15 }],
  ],
  texts: [],
  blocks: [["OX", 0.12, 0.42, 6.76, 3.2, [3.95, 3.43]], ["CU", 0.12, 3.65, 4.0, 1.22, [0.15, 4.67], "Cu–氧化物界面 / 氧化物（RWGS · CO 加氢）"],
           ["EC", 2.5, 4.92, 4.38, 1.06, [2.55, 4.94]]],
  notes: ["In₂O₃：Vₒ（氧空位）是 CO₂ 吸附与活化位；H₂ 在 In–O 对上异裂（In–H⁻ / O–H⁺）。", "甲酸盐路线 HCOO* → H₂COO* → H₂COOH* → H₂CO* → CH₃O* → CH₃OH，每步交替得到 H⁻ / H⁺。",
          "CH₃OH 脱附后 O 留在晶格，H₂ 再把它还原成 H₂O 并重新生成 Vₒ——氧空位循环。", "RWGS 两条路：缔合（COOH*）与氧化还原（CO₂ 把 O 填入 Vₒ）。高温时 RWGS 份额上升（表观 Ea ≈ 117 vs 103 kJ/mol）。",
          "GaZrOₓ：Ga–H 提高 HCOO* 覆盖（郭组 ACS Catal. 2025）。"],
});
// ---------------------------------------------------------------- B: metal C1 core
MODS.push({ id: "B", color: C.accent1, W: 7.0, H: 6.05,
  title: "B  金属位 C_1 核心：CO 活化 · CH_x · 甲烷化 · 重整（SMR / DRM / POM）· 积碳 / 碳化物",
  nodes: [
    ["COg", "mol", 0.5, 1.0, "CO(g)"], ["CO", "surf", 1.7, 1.0, "CO*"], ["HCO", "surf", 3.0, 1.0, "HCO*"], ["H2CO", "surf", 4.3, 1.0, "H_2CO*"],
    ["CH3O", "surf", 5.6, 1.0, "CH_3O*"], ["MeOH", "mol", 6.55, 1.75, "CH_3OH(g)"],
    ["COH", "surf", 2.45, 2.3, "COH*"], ["HCOH", "surf", 3.65, 2.3, "HCOH*"],
    ["CO2g", "mol", 0.5, 2.3, "CO_2(g)"], ["COOH", "surf", 1.1, 3.1, "COOH*"],
    ["C", "surf", 1.9, 4.0, "C*"], ["CH", "surf", 3.05, 4.0, "CH*"], ["CH2", "surf", 4.3, 4.0, "CH_2*"], ["CH3", "surf", 5.55, 4.0, "CH_3*"],
    ["CH4", "mol", 6.6, 4.0, "CH_4(g)"],
    ["coke", "mol", 0.6, 4.0, "积碳 / 石墨 C"], ["FeC", "surf", 0.75, 4.85, "χ-Fe_5C_2 碳化物"],
    ["tCH2", "txt", 4.3, 4.6, "→ D 单体", { tc: C.accent2, bold: true }], ["tCH3", "txt", 5.55, 4.6, "→ D 链引发", { tc: C.accent2, bold: true }],
    ["tCO", "txt", 1.7, 0.45, "→ D（CO 插入）", { tc: C.accent2, bold: true }],
    ["H2", "mol", 3.0, 5.5, "H_2(g)"], ["H", "surf", 3.95, 5.5, "H*"], ["H2O", "mol", 4.9, 5.5, "H_2O(g)"], ["OH", "surf", 5.8, 5.5, "OH*"],
    ["O", "surf", 6.6, 5.5, "O*"], ["O2", "mol", 6.6, 4.75, "O_2(g)"],
  ],
  edges: [
    ["COg", "CO", "SH", "", { both: true }], ["CO", "HCO", "R3", "+H*", { ly: -0.15 }], ["HCO", "H2CO", "R3", "+H*", { ly: -0.15 }],
    ["H2CO", "CH3O", "R3", "+H*", { ly: -0.15 }], ["CH3O", "MeOH", "R3", "+H*（Cu 类）", { perp: 0.16 }],
    ["CO", "COH", "R3", "+H*", { perp: 0.15 }], ["HCO", "HCOH", "R3", "+H*", { perp: 0.15 }],
    ["CO", "C", "SH", "直接解离 −O*", { both: true, lt: 0.62, perp: -0.32 }], ["COH", "C", "SH", "−OH*", { perp: 0.15 }],
    ["HCO", "CH", "SH", "−O*", { both: true, lt: 0.62, perp: -0.15 }], ["HCOH", "CH", "SH", "−OH*", { perp: 0.15, lt: 0.6 }],
    ["H2CO", "CH2", "SH", "−O*", { both: true, lt: 0.6, perp: -0.15 }], ["CH3O", "CH3", "SH", "−O*", { both: true, lt: 0.6, perp: -0.15 }],
    ["C", "CH", "R3", "+H*", { both: true, ly: -0.15 }], ["CH", "CH2", "R3", "+H*", { both: true, ly: -0.15 }],
    ["CH2", "CH3", "R3", "+H*", { both: true, ly: -0.15 }], ["CH3", "CH4", "R3", "+H* ⇌ C–H 活化", { both: true, ly: -0.2 }],
    ["CO2g", "CO", "R1", "DRM：CO_2* → CO* + O*", { perp: 0.17, lt: 0.45 }], ["CO", "COOH", "SH", "+OH*（WGS）", { perp: -0.2, lt: 0.6 }],
    ["COOH", "CO2g", "SH", "−H*", { perp: 0.14 }],
    ["C", "coke", "R1", "析碳", { ly: -0.15 }], ["C", "FeC", "R2", "渗碳 ⇌ 去碳", { both: true, perp: 0.18, lt: 0.6 }],
    ["CH2", "tCH2", "R2", ""], ["CH3", "tCH3", "R2", ""], ["CO", "tCO", "R2", ""],
    ["H2", "H", "SH", "解离"], ["H2O", "OH", "R1", "−H*（SMR）", { both: true, ly: -0.15 }], ["OH", "O", "SH", "−H*", { both: true, ly: -0.15 }],
    ["O2", "O", "R1", "POM 解离", { perp: -0.12 }],
  ],
  blocks: [["CU", 2.45, 0.55, 4.5, 1.45, [2.5, 0.57], "Cu 基 / Cu–ZnO（不解离加氢 → 甲醇）"], ["MET", 1.2, 0.55, 1.2, 1.5, null],
           ["MET", 1.42, 2.05, 5.53, 2.4, [1.47, 4.22]], ["SUP", 0.08, 2.05, 1.32, 1.45, [0.1, 3.28], "载体 CeO₂ / ZrO₂"],
           ["CARB", 0.08, 3.7, 1.32, 1.62, [0.1, 5.12], "碳化物 χ-Fe₅C₂"], ["MET", 2.45, 5.12, 4.48, 0.85, [2.5, 5.14], "金属 + 载体：H₂ 解离 · H₂O 活化 · O 移除"]],
  texts: [[0.2, 5.42, 2.2, 0.6, ["Boudouard：2CO ⇌ C* + CO_2", "SMR：CH_4 + H_2O ⇌ CO + 3H_2", "DRM：CH_4 + CO_2 ⇌ 2CO + 2H_2", "POM：CH_4 + ½O_2 → CO + 2H_2"], { valign: "top", color: C.accent1, fs: 5 }]],
  notes: ["CO 活化三条路：直接解离（CO* → C* + O*）、H 助解离（经 COH* / HCO* / HCOH*）、加氢不解离（→ CH₃O* → 甲醇，Cu 类）。",
          "C* ⇌ CH* ⇌ CH₂* ⇌ CH₃* ⇌ CH₄：正向为甲烷化，逆向为甲烷 C–H 活化（重整第一步）。", "Fe 催化剂上 C* 与 χ-Fe₅C₂ 互换（渗碳 / 去碳），CH₂* 与 CH₃* 供给链增长（模块 D）。",
          "O* / OH* 由 H* 移除成 H₂O；重整时方向相反（H₂O → OH* + H* → O*）。"],
});
// ---------------------------------------------------------------- C: radicals
MODS.push({ id: "C", color: C.accent6, dash: "dash", fill: "FFF7F9", W: 9.2, H: 6.05,
  title: "C  气相 / 均相自由基网络：引发 · 抽氢 · 加成 · β-断裂 · 氧化 · 支化 · 终止（R• 与 A / B 中的 R* 对应）",
  nodes: [
    ["CO2", "mol", 0.5, 1.35, "CO_2"], ["HOCO", "rad", 1.55, 0.85, "HOCO•"], ["HCOO", "rad", 1.55, 1.85, "HCOO•"],
    ["CO", "mol", 2.7, 0.85, "CO"], ["HCOOH", "mol", 2.7, 1.85, "HCOOH"], ["HCO", "rad", 3.85, 0.85, "HCO•"], ["HOCHOH", "rad", 3.7, 1.85, "HOC•HOH"],
    ["CH2O", "mol", 4.7, 2.85, "CH_2O"], ["CH3O", "rad", 3.4, 2.55, "CH_3O•"], ["CH2OH", "rad", 3.4, 3.25, "•CH_2OH"],
    ["MeOH", "mol", 2.1, 2.9, "CH_3OH"], ["tMe", "txt", 0.75, 2.9, "→ •CH_3 + OH•", { tc: C.accent6, bold: true, w: 0.95 }],
    ["CH4", "mol", 5.35, 0.85, "CH_4"], ["CH3", "rad", 6.55, 0.85, "•CH_3"], ["C2H6", "mol", 7.75, 0.85, "C_2H_6"], ["C2H5", "rad", 8.7, 1.6, "C_2H_5•"],
    ["C2H4", "mol", 7.75, 2.35, "C_2H_4"], ["C2H3", "rad", 8.7, 3.1, "C_2H_3•"], ["C2H2", "mol", 7.75, 3.4, "C_2H_2"],
    ["CH3OO", "rad", 6.55, 2.0, "CH_3OO•"],
    ["H2", "mol", 0.6, 4.55, "H_2"], ["H", "rad", 1.8, 4.55, "H•"], ["OH", "rad", 3.0, 4.55, "OH•"], ["H2O", "mol", 3.0, 5.5, "H_2O"],
    ["Orad", "rad", 4.2, 4.55, "O•"], ["O2", "mol", 5.4, 4.55, "O_2"], ["HO2", "rad", 5.4, 5.5, "HO_2•"], ["H2O2", "mol", 6.75, 5.5, "H_2O_2"],
  ],
  edges: [
    ["CO2", "HOCO", "R6", "+H•", { perp: 0.13 }], ["CO2", "HCOO", "R6", "+H•", { perp: -0.13 }], ["HOCO", "HCOO", "R6", "H 迁移", { both: true, perp: 0.22 }],
    ["HOCO", "CO", "R6", "−OH• / +OH•", { both: true, ly: -0.15 }], ["HCOO", "HCOOH", "R6", "+H•", { ly: -0.15 }], ["HOCO", "HCOOH", "R6", "+H•", { perp: 0.14, lt: 0.4 }],
    ["CO", "HCO", "R6", "+H• / −H•", { both: true, ly: -0.15 }], ["HCOOH", "HOCHOH", "R6", "+H•", { ly: -0.15 }],
    ["HCO", "CH2O", "R6", "+H• / 抽氢", { both: true, perp: -0.2, lt: 0.35 }], ["HOCHOH", "CH2O", "R6", "−OH•", { perp: 0.16 }],
    ["CH2O", "CH3O", "R6", "+H•", { perp: 0.14 }], ["CH2O", "CH2OH", "R6", "+H•", { perp: -0.14 }], ["CH3O", "CH2OH", "R6", "", { both: true }],
    ["CH3O", "MeOH", "R6", "+H•", { perp: 0.14 }], ["CH2OH", "MeOH", "R6", "+H• / 抽氢", { both: true, perp: -0.16 }], ["MeOH", "tMe", "R6", "C–O 均裂", { ly: -0.15 }],
    ["CH4", "CH3", "R6", "+H•/OH•/O• 抽氢；均裂", { both: true, ly: -0.16 }], ["CH3", "C2H6", "R6", "2•CH_3 偶联（OCM）", { ly: -0.16 }],
    ["C2H6", "C2H5", "R6", "抽氢 / 歧化", { both: true, perp: 0.18 }], ["C2H5", "C2H4", "R6", "β-断裂 −H• / +H•", { both: true, perp: 0.2 }],
    ["C2H4", "C2H3", "R6", "抽氢", { perp: -0.14 }], ["C2H3", "C2H2", "R6", "−H•", { perp: 0.13 }],
    ["CH3", "CH3OO", "R6", "+O_2", { perp: 0.14 }], ["CH3OO", "CH2O", "R6", "→ CH_2O + OH•", { perp: -0.15, lt: 0.45 }],
    ["H2", "H", "R6", "均裂 / 复合", { both: true, ly: -0.15 }], ["OH", "H", "R6", "+H_2 → H_2O + H•", { ly: -0.15 }], ["OH", "H2O", "R6", "+H_2 / +RH（抽氢）", { perp: 0.2 }],
    ["H", "H2O", "R6", "+OH•（终止）", { perp: -0.18 }], ["Orad", "OH", "R6", "+H_2 / +CH_4", { ly: -0.15 }],
    ["O2", "Orad", "R6", "+H• → OH• + O•（支化）", { ly: -0.16 }], ["O2", "HO2", "R6", "+H• (+M)", { perp: 0.14 }],
    ["HO2", "OH", "R6", "+H• → 2OH•", { perp: 0.16 }], ["HO2", "H2O2", "R6", "+HO_2• / RH", { ly: -0.15 }],
  ],
  texts: [[0.2, 0.4, 3.0, 0.25, "C_1 加氢自由基链（CO_2 → CH_3OH）", { color: C.accent6, bold: true, fs: 6.5 }],
          [5.0, 0.4, 3.5, 0.25, "CH_4 活化 / OCM / POM / C_2", { color: C.accent6, bold: true, fs: 6.5 }],
          [0.2, 4.0, 4.5, 0.25, "H / O 自由基池：引发 · 传递 · 支化 · 终止", { color: C.accent6, bold: true, fs: 6.5 }],
          [7.35, 5.85, 1.8, 0.2, "H_2O_2 → 2OH•（高温）", { color: C.accent6 }], [5.0, 3.75, 3.6, 0.2, "另：•CH_3 + O• → CH_2O + H•", { color: C.accent6 }]],
  blocks: [["GAS", 0.1, 0.42, 9.0, 5.56, [5.6, 4.0]]],
  notes: ["左上：CO₂ + H• 生成 HOCO•（羧基）或 HCOO•（甲酸根），分别通向 CO 与 HCOOH；HCO• → CH₂O → CH₃O• / •CH₂OH → CH₃OH。",
          "右上：CH₄ 被 H• / OH• / O• 抽氢或高温均裂成 •CH₃，•CH₃ 偶联得 C₂H₆（OCM），再脱氢到 C₂H₄、C₂H₂；•CH₃ + O₂ → CH₃OO• → CH₂O（POM）。",
          "下方：H• + O₂ → OH• + O• 是链支化步；HO₂•、H₂O₂ 是低温储存态；H• + OH• → H₂O 是终止。", "C1 显式枚举：10 个自由基、141 个反应（data/network_full_C1.csv）。"],
});
// ---------------------------------------------------------------- D: FT chain growth (explicit C1-C4)
const dX = (n) => 1.05 + (n - 1) * 1.35;
const dN = [], dE = [];
{
  const L = { 1: ["CH_3*", "CH_4", null, null], 2: ["C_2H_5*", "C_2H_6", "C_2H_4*", "C_2H_4"], 3: ["n-C_3H_7*", "C_3H_8", "C_3H_6*", "C_3H_6"],
              4: ["n-C_4H_9*", "n-C_4H_{10}", "C_4H_8*", "1-C_4H_8"], 5: ["C_nH_{2n+1}*", "C_nH_{2n+2}", "C_nH_{2n}*", "α-C_nH_{2n}"] };
  const AC = { 1: ["CH_3CO*", "C_2H_5OH"], 2: ["C_2H_5CO*", "C_3H_7OH"], 3: ["C_3H_7CO*", "C_4H_9OH"], 4: ["C_4H_9CO*", "C_5H_{11}OH"], 5: ["C_{n+1} 酰基*", "C_{n+1} 醇"] };
  for (let n = 1; n <= 5; n++) {
    const x = dX(n), [al, pa, os, og] = L[n], [ac, alc] = AC[n];
    dN.push([`A${n}`, "surf", x, 2.6, al], [`P${n}`, "mol", x, 1.6, pa], [`AC${n}`, "surf", x + 0.33, 3.5, ac], [`AL${n}`, "mol", x + 0.33, 4.4, alc]);
    if (os) dN.push([`OS${n}`, "surf", x - 0.33, 3.5, os], [`OG${n}`, "mol", x - 0.33, 4.4, og]);
    dE.push([`A${n}`, `P${n}`, "R2", n === 1 ? "+H*（终止）" : "+H*", { perp: 0.15 }]);
    dE.push([`A${n}`, `AC${n}`, "R2", n === 1 ? "+CO*（插入）" : "+CO*", { perp: -0.17, lt: 0.45 }]);
    dE.push([`AC${n}`, `AL${n}`, "R3", n === 1 ? "+3H*" : "", { perp: 0.17 }]);
    if (os) { dE.push([`A${n}`, `OS${n}`, "R2", n === 2 ? "β-H 消除" : "", { perp: 0.17, lt: 0.45 }]); dE.push([`OS${n}`, `OG${n}`, "SH", n === 2 ? "脱附 / 再吸附" : "", { both: true, perp: -0.22 }]); }
    if (n < 5) dE.push([`A${n}`, `A${n + 1}`, "R2", n === 4 ? "… +CH_2*" : "+CH_2*", { ly: -0.15 }]);
  }
  dN.push(["mono", "surf", 0.42, 3.5, "CH_2* 单体"], ["COs", "surf", 0.42, 4.4, "CO*"], ["CHs", "surf", 0.42, 5.3, "CH*"],
          ["vin", "surf", 1.45, 5.3, "CH=CH_2* 乙烯基"], ["ally", "surf", 2.95, 5.3, "CH_2CH=CH_2* 烯丙基"],
          ["iC4", "surf", 4.43, 1.95, "i-C_4H_9*"], ["iC4p", "mol", 4.43, 1.2, "i-C_4H_{10} / i-C_4H_8"], ["C4i", "mol", 4.77, 5.3, "2-C_4H_8"],
          ["wax", "mol", 6.4, 5.3, "C_{20+} 蜡 → F"], ["tB", "txt", 0.42, 2.95, "← B", { tc: C.accent1, bold: true, w: 0.4 }]);
  dE.push(["CHs", "vin", "R2", "+CH_2*", { ly: -0.15 }], ["vin", "OS2", "R3", "+H*", { perp: 0.15, lt: 0.8 }], ["vin", "ally", "R2", "+CH_2*", { ly: -0.15 }],
          ["ally", "OS3", "R3", "+H*", { perp: -0.15, lt: 0.5 }], ["A3", "iC4", "R2", "+CH_3*（支化）", { perp: 0.17 }], ["iC4", "iC4p", "R2", "", {}],
          ["OG4", "C4i", "SH", "双键异构", { perp: -0.16 }], ["A5", "wax", "R2", "链增长 → C_{20+}", { pts: [[6.9, 2.6], [7.25, 2.6], [7.25, 5.3], [7.0, 5.3]], lpos: [7.0, 4.95] }]);
}
MODS.push({ id: "D", color: C.accent2, W: 7.4, H: 6.05,
  title: "D  金属位 C–C 链增长（Fe_5C_2 / Co_2C）：CH_x 插入 · CO 插入 · 烯基机理（C_1 → C_4 显式 + 通式）",
  nodes: dN, edges: dE,
  blocks: [["CARB", 0.1, 0.42, 7.2, 5.56, [0.15, 0.44]],
           ...[1, 2, 3, 4, 5].map((n) => ["DUAL", dX(n) + 0.33 - (n === 5 ? 0.45 : 0.41), 3.3, n === 5 ? 0.9 : 0.82, 1.3, n === 1 ? [0.15, 5.64] : null])],
  texts: [[3.85, 5.7, 3.4, 0.25, "ASF：W_n = n(1−α)^{2} α^{n−1}，α = r_p / (r_p + r_t)", { color: C.accent2 }]],
  notes: ["碳化物机理：CHₓ* 逐个插入烷基链（CH₃* → C₂H₅* → C₃H₇* → …），每个链节点都有三条出路：+H* 终止成烷烃、β-H 消除成 α-烯烃、+CO* 插入成酰基 → 醇。",
          "烯基机理：CH* + CH₂* → CH=CH₂* → 烯丙基 → 丙烯，解释 Fe 上较高的烯烃选择性。", "CH₃* 插入仲碳得到甲基支链（i-C₄）；α-烯烃再吸附可重新引发链或双键异构成内烯烃。",
          "郭组 Fe 基 CO₂-FTS：K/Na 助剂抑制 +H* 终止、促进 β-H 消除，提高 C₂–C₄= 与 C₅+。"],
});
// ---------------------------------------------------------------- E: zeolite dual cycle
const eN = [
  ["MeOH", "mol", 0.55, 0.75, "CH_3OH(g)"], ["MeOH2", "ion", 1.85, 0.75, "CH_3OH_2^{+}–Z"], ["SMS", "surf", 3.25, 0.75, "Z–OCH_3（SMS）"],
  ["DME", "mol", 4.55, 0.75, "DME"], ["TMO", "ion", 5.8, 0.75, "(CH_3)_3O^{+}"], ["Acyl", "surf", 7.15, 0.75, "CH_3CO–Z 表面乙酰基"],
  ["Ket", "surf", 8.75, 0.75, "CH_2=C=O 乙烯酮"], ["COg", "mol", 7.15, 1.45, "CO（来自 A）"], ["MA", "mol", 8.75, 1.45, "乙酸甲酯"],
  ["HCHO", "mol", 3.25, 1.45, "HCHO + CH_4"], ["MMC", "ion", 4.55, 1.45, "CH_3OCH_2^{+}"], ["Et", "surf", 5.85, 1.45, "CH_3OCH_2CH_3 / Z–OC_2H_5"],
];
const OLE = ["C_2H_4", "C_3H_6", "C_4H_8", "C_5H_{10}", "C_6H_{12}", "C_{7+}^{=}"], CAR = ["C_3H_7^{+}", "C_4H_9^{+}", "C_5H_{11}^{+}", "C_6H_{13}^{+}", "C_7H_{15}^{+}"],
      PAR = ["C_2H_6", "C_3H_8", "C_4H_{10}", "C_5H_{12}", "C_6H_{14}"];
const eE = [
  ["MeOH", "MeOH2", "R4", "质子化（B 酸）", { ly: -0.15 }], ["MeOH2", "SMS", "R4", "−H_2O", { ly: -0.15 }], ["SMS", "DME", "R4", "+CH_3OH", { both: true, ly: -0.15 }],
  ["DME", "TMO", "R4", "+SMS", { ly: -0.15 }], ["SMS", "HCHO", "R4", "+CH_3OH 氢转移", { perp: 0.17 }], ["DME", "MMC", "R4", "+SMS −CH_4", { perp: -0.17 }],
  ["TMO", "MMC", "R4", "−H^{+}（叶立德）", { perp: 0.16 }], ["MMC", "Et", "R4", "+CH_3OH / DME", { ly: 0.17 }],
  ["COg", "Acyl", "R4", "+SMS（Koch）", { perp: 0.15 }], ["Acyl", "Ket", "R4", "−H^{+}", { both: true, ly: -0.15 }], ["Acyl", "MA", "R4", "+CH_3OH", { perp: -0.15 }],
  ["Et", "O0", "R4", "第一个 C–C 键 → C_2H_4", { pts: [[5.85, 1.58], [5.85, 1.98], [0.62, 1.98], [0.62, 2.27]], lpos: [3.6, 1.98] }],
  ["Ket", "O0", "R4", "乙烯酮 → C_2H_4（OEO）", { pts: [[9.67, 0.75], [10.15, 0.75], [10.15, 1.83], [0.46, 1.83], [0.46, 2.27]], lpos: [9.35, 1.83], dash: "dash" }],
];
OLE.forEach((l, i) => eN.push([`O${i}`, "mol", 0.55 + 1.3 * i, 2.4, l]));
CAR.forEach((l, i) => eN.push([`K${i}`, "ion", 1.2 + 1.3 * i, 3.2, l]));
PAR.forEach((l, i) => eN.push([`P${i}`, "mol", 0.55 + 1.3 * i, 4.0, l]));
for (let i = 0; i < 5; i++) {
  eE.push([`O${i}`, `K${i}`, "R4", "", {}]);
  eE.push([`K${i}`, `O${i + 1}`, "R4", "", {}]);
  eE.push([`O${i}`, `P${i}`, "R3", "", {}]);
}
eN.push(["Crack", "mol", 8.45, 2.4, "β-断裂碎片 C_3^{=} / C_4^{=}"]);
eE.push(["K3", "Crack", "R4", "β-断裂", { perp: 0.14, lt: 0.6 }], ["K4", "Crack", "R4", "", {}], ["O5", "Crack", "R4", "", {}],
        ["Crack", "O1", "R4", "回到烯烃循环", { pts: [[8.45, 2.27], [8.45, 2.13], [1.95, 2.13], [1.95, 2.27]], lpos: [7.2, 2.13] }],
        ["K0", "K3", "R4", "齐聚：C_3^{+} + C_3^{=} → C_6^{+}", { pts: [[1.2, 3.34], [1.2, 4.42], [5.1, 4.42], [5.1, 3.34]], lpos: [3.0, 4.42] }]);
const ARO = [["Cyc", "ion", 1.0, "环烯基^{+} C_6–C_8"], ["Bz", "mol", 2.3, "苯"], ["Tol", "mol", 3.25, "甲苯"], ["Xyl", "mol", 4.4, "二甲苯（o/m/p）"],
             ["TMB", "mol", 5.65, "三甲苯"], ["TeMB", "mol", 6.6, "四甲苯"], ["PMB", "mol", 7.55, "五甲苯"], ["HMB", "mol", 8.5, "六甲苯"], ["Hepta", "ion", 9.6, "七甲基苯鎓^{+}"]];
ARO.forEach(([id, t, x, l]) => eN.push([id, t, x, 5.05, l]));
for (let i = 1; i < ARO.length; i++) eE.push([ARO[i - 1][0], ARO[i][0], i === 1 ? "R4" : "R4", i === 1 ? "−3H" : i === 2 ? "+CH_3^{+}" : "", { ly: -0.15 }]);
eE.push(["K4", "Cyc", "R4", "环化（C_6–C_8 碳正离子）", { pts: [[6.4, 3.34], [6.4, 4.62], [1.0, 4.62], [1.0, 4.92]], lpos: [5.25, 4.62] }]);
eN.push(["CP", "ion", 9.6, 4.3, "环戊烯基^{+}"], ["Side", "ion", 8.5, 5.75, "外亚甲基 / 乙基侧链^{+}"], ["tSide", "txt", 7.05, 5.75, "→ C_2H_4 / C_3H_6", { tc: C.accent4, bold: true }],
        ["Naph", "mol", 4.4, 5.75, "甲基萘"], ["Coke", "mol", 5.75, 5.75, "多环芳烃 / 积碳"],
        ["tHT", "txt", 1.65, 5.8, "−3H：氢转移（3 烯烃 → 3 烷烃）或 Ga / Zn 脱氢 −H_2", { tc: C.accent4, w: 2.6, align: "left" }]);
eE.push(["Hepta", "CP", "R4", "缩环（paring）", { perp: -0.15 }], ["CP", "Crack", "R4", "释放 C_3^{=} / i-C_4^{=}", { perp: -0.18, lt: 0.45 }],
        ["HMB", "Side", "R4", "侧链甲基化", { perp: -0.16 }], ["Side", "tSide", "R4", ""], ["TMB", "Naph", "R4", "缩合环化", { perp: 0.16 }], ["Naph", "Coke", "R4", "", {}]);
MODS.push({ id: "E", color: C.accent4, W: 10.4, H: 6.05,
  title: "E  分子筛酸位：甲醇活化 · 第一个 C–C 键 · 烯烃循环（显式 C_2–C_7）· 芳烃循环（苯 → 七甲基苯鎓⁺）· 氢转移 / 积碳",
  nodes: eN, edges: eE,
  blocks: [["BA", 0.1, 0.42, 10.2, 5.56, [0.15, 1.25], "Brønsted 酸分子筛（H-ZSM-5 / SAPO-34）"], ["OX", 6.35, 0.45, 3.95, 1.5, [9.3, 1.58], "氧化物（OEO）"],
           ["OX", 0.12, 3.78, 7.1, 0.52, [6.45, 3.95], "氧化物加氢"], ["LA", 0.32, 4.75, 2.43, 0.58, [0.35, 5.36], "Lewis 酸 Ga / Zn（脱氢）"]],
  texts: [[7.35, 3.55, 2.9, 0.62, ["↘ 甲基化 +CH_3^{+}（SMS / CH_3OH / DME）", "↗ 去质子 −H^{+}", "↓ +H（氧化物 H 溢流）/ 氢转移 → 烷烃"], { color: C.accent4 }]],
  notes: [],
});
// ---------------------------------------------------------------- F: acid upgrading to jet
MODS.push({ id: "F", color: C.accent5, W: 5.0, H: 6.05,
  title: "F  酸位升级：双功能加氢异构 / 裂化 · 烯烃齐聚制航煤",
  nodes: [
    ["NC12", "mol", 0.7, 0.8, "n-C_{12}H_{26}（← D）"], ["C12e", "mol", 2.1, 0.8, "n-C_{12}H_{24}"], ["Sec", "ion", 3.6, 0.8, "sec-C_{12}H_{25}^{+}"],
    ["Mono", "ion", 3.6, 1.8, "单支链 C_{12}^{+}"], ["Di", "ion", 3.6, 2.8, "双支链 C_{12}^{+}"], ["Tri", "ion", 3.6, 3.8, "三支链 C_{12}^{+}"],
    ["I1", "mol", 1.25, 1.8, "单支链 i-C_{12}"], ["I2", "mol", 1.25, 2.8, "双支链 i-C_{12}"], ["I3", "mol", 1.25, 3.8, "三支链 i-C_{12}"],
    ["Frag", "ion", 3.6, 4.8, "碎片 C_5–C_7^{+} + C_5–C_7^{=}"], ["Light", "mol", 3.6, 5.6, "C_5–C_8 异构烷烃"],
    ["Jet", "mol", 1.25, 4.8, "C_8–C_{16} 异构烷烃航煤（SPK）", { fill: C.accent5 }],
    ["C34", "mol", 0.65, 5.6, "C_3^{=} / C_4^{=}（← E）"], ["Olig", "mol", 2.05, 5.6, "C_8–C_{16}^{=}"],
  ],
  edges: [
    ["NC12", "C12e", "R5", "金属 −H_2 / +H_2", { both: true, ly: -0.15 }], ["C12e", "Sec", "R5", "+H^{+} / −H^{+}", { both: true, ly: -0.15 }],
    ["Sec", "Mono", "R5", "PCP 支化", { perp: -0.25 }], ["Mono", "Di", "R5", "PCP", { perp: -0.15 }], ["Di", "Tri", "R5", "PCP", { perp: -0.15 }],
    ["Mono", "I1", "R5", "−H^{+} → +H_2（金属）", { ly: -0.15 }], ["Di", "I2", "R5", "", {}], ["Tri", "I3", "R5", "", {}],
    ["Tri", "Frag", "R5", "A 型（最快）", { perp: -0.18 }], ["Di", "Frag", "R5", "B 型", { ax: 0.35, bx: 0.35, perp: -0.16, lt: 0.3 }],
    ["Mono", "Frag", "R5", "C 型（慢）", { ax: 0.62, bx: 0.62, perp: -0.2, lt: 0.2 }],
    ["Frag", "Light", "R5", "−H^{+} +H_2", { perp: -0.18 }], ["I3", "Jet", "R5", "→ 航煤", { perp: 0.15 }],
    ["C34", "Olig", "R5", "酸催化齐聚", { ly: -0.15 }], ["Olig", "Jet", "R3", "+H_2", { perp: -0.15 }],
  ],
  texts: [[0.2, 4.35, 2.2, 0.25, "单 / 双 / 三支链 i-C_{12} 合并为航煤馏分", { color: C.accent5 }]],
  blocks: [["BIF", 0.1, 0.42, 4.8, 4.75, [0.15, 4.97]], ["MET", 0.12, 0.5, 2.78, 0.6, [0.15, 0.95], "金属位：脱氢 / 加氢"],
           ["MET", 0.5, 1.55, 1.5, 2.5, [0.52, 1.52], "金属位：加氢"], ["BA", 2.95, 0.5, 1.93, 4.55, [2.97, 4.95], "酸位：质子化 · PCP · β-断裂"],
           ["BA", 0.1, 5.25, 2.85, 0.72, [0.12, 5.27], "酸：齐聚"]],
  notes: ["双功能机理：金属位脱氢成烯烃 → 酸位质子化成仲碳正离子 → PCP（质子化环丙烷）逐级支化 → 去质子、金属位加氢得到异构烷烃。",
          "β-断裂速率 A（叔→叔）≫ B（仲↔叔）> C（仲→仲）≫ D：支化越多越容易裂解，异构与裂化存在竞争。",
          "烯烃齐聚制航煤（OTJ）：C₃= / C₄=（来自模块 E 或 FTO）在酸上二聚 / 三聚成 C₈–C₁₆ 烯烃，再加氢得到 SPK。"],
});

// ================================================================== OVERVIEW (previous poster without carbonylation)
const OV = { id: "OV", W: 20, H: 12.5, nodes: [], edges: [], zones: [], texts: [] };
{
  const n = (id, type, x, y, label, w, o = {}) => OV.nodes.push([id, type, x, y, label, { w, h: type === "rad" ? 0.42 : 0.44, fs: type === "hub" ? 13 : 12, ...o }]);
  n("fCO2", "mol", 0.95, 3.05, "CO_2", 1.05); n("fH2", "mol", 0.95, 4.35, "H_2", 1.05); n("fH2O", "mol", 0.95, 5.4, "H_2O", 1.05);
  n("fO2", "mol", 0.95, 6.45, "O_2", 1.05); n("fCH4", "mol", 0.95, 7.55, "CH_4", 1.05);
  n("eCO2m", "rad", 3.0, 1.75, "CO_2•^{−}", 1.05); n("mHCOOH", "mol", 4.6, 1.75, "HCOOH", 1.15);
  n("sCO2", "surf", 3.35, 3.05, "CO_2*", 1.05); n("sHCOO", "surf", 4.6, 3.05, "HCOO*", 1.05); n("sH2COO", "surf", 5.85, 3.05, "H_2COO*", 1.1);
  n("sH2COOH", "surf", 7.08, 3.05, "H_2COOH*", 1.2); n("sH2CO", "surf", 8.25, 3.05, "H_2CO*", 0.95); n("sCH3O", "surf", 9.25, 3.05, "CH_3O*", 0.9);
  n("sH", "surf", 3.0, 4.35, "H*", 0.8); n("sOH", "surf", 3.0, 5.4, "OH*", 0.8); n("sO", "surf", 3.0, 6.45, "O*", 0.8);
  n("sCOOH", "surf", 4.6, 4.35, "COOH*", 1.05); n("hubCO", "hub", 6.45, 4.25, "CO(g) 合成气", 1.7); n("sHCO", "surf", 8.8, 4.35, "HCO*", 0.9);
  n("sCO", "surf", 7.9, 5.5, "CO*", 0.9); n("sCOH", "surf", 6.55, 6.3, "COH*", 0.9);
  n("sC", "surf", 7.9, 7.45, "C*", 0.8); n("sCH", "surf", 6.6, 7.45, "CH*", 0.8); n("sCH2", "surf", 5.3, 7.45, "CH_2*", 0.85); n("sCH3", "surf", 4.0, 7.45, "CH_3*", 0.85);
  n("tCO2H", "txt", 2.85, 9.45, "CO_2 + H•", 0.95, { tc: C.accent6, bold: true, italic: false, fs: 11 });
  n("rHOCO", "rad", 3.95, 9.1, "HOCO•", 1.0); n("rHCOO", "rad", 3.95, 9.8, "HCOO•", 1.0); n("rHCO", "rad", 5.25, 9.1, "HCO•", 1.0);
  n("rCH2O", "mol", 6.45, 9.1, "CH_2O", 0.9); n("rCH3O", "rad", 7.65, 9.1, "CH_3O•", 1.0); n("rCH2OH", "rad", 7.65, 9.8, "•CH_2OH", 1.0);
  n("rMeOH", "mol", 9.0, 9.45, "CH_3OH", 1.0); n("rHCOOH", "txt", 5.25, 9.8, "→ HCOOH", 1.0, { tc: C.accent6, bold: true, italic: false, fs: 11 });
  n("rCH3", "rad", 3.1, 10.5, "•CH_3", 1.0); n("rC2", "mol", 4.55, 10.5, "C_2H_6/C_2H_4", 1.5); n("rCH3OO", "rad", 6.45, 10.5, "CH_3OO•", 1.05);
  n("rH", "rad", 7.6, 10.5, "H•", 0.8); n("rHO2", "rad", 8.5, 10.5, "HO_2•", 0.85); n("rOH", "rad", 9.35, 10.5, "OH•", 0.8);
  n("mMeOH", "mol", 10.75, 2.35, "CH_3OH", 1.05); n("zOCH3", "surf", 12.25, 2.35, "Z–OCH_3 表面甲氧基", 1.85);
  n("mDME", "mol", 13.95, 2.35, "DME", 0.85); n("mARO", "mol", 15.65, 2.35, "芳烃 BTX / C_9–C_{10}", 1.95);
  n("zKET", "surf", 11.3, 3.9, "CH_2=C=O 乙烯酮", 1.65); n("iPOOL", "ion", 14.0, 3.9, "烃池：多甲基苯^{+} / 环戊烯基^{+}", 2.75);
  n("mC24", "mol", 17.95, 3.9, "C_2–C_4 烯烃", 1.45); n("zOEt", "surf", 12.55, 5.35, "Z–OC_2H_5 乙氧基", 1.75);
  n("mC5O", "mol", 15.6, 5.35, "C_{5+} 烯烃", 1.25); n("mGAS", "mol", 17.95, 5.35, "C_5–C_{11} 汽油", 1.6);
  n("mEtOH", "mol", 11.0, 7.52, "C_2H_5OH", 1.1); n("sCH3CO", "surf", 11.0, 8.24, "CH_3CO*", 1.05);
  n("sCnH", "surf", 13.35, 9.15, "C_nH_{2n+1}* 烷基链", 1.85); n("mAOle", "mol", 15.6, 7.6, "α-烯烃 C_nH_{2n}（FTO）", 2.05);
  n("mNPar", "mol", 15.75, 9.15, "正构烷烃 C_nH_{2n+2}", 1.85); n("mWax", "mol", 15.75, 10.4, "C_{20+} 蜡", 1.25); n("mHA", "mol", 12.3, 10.4, "C_{2+} 醇（CO 插入）", 1.85);
  n("iCARB", "ion", 18.0, 8.3, "C_nH_{2n+1}^{+} 碳正离子", 2.15); n("mJET", "mol", 18.0, 10.2, "C_8–C_{16} 航煤 / 柴油", 2.05);
  const e = (a, b, cat, label, o = {}) => OV.edges.push([a, b, cat, label, { width: 1.6, fs: 10, ...o }]);
  e("fCO2", "sCO2", "SH", "吸附", { lt: 0.45 }); e("fH2", "sH", "SH", "解离", { lt: 0.45 }); e("fH2O", "sOH", "R1", "SMR", { both: true, lt: 0.45 });
  e("fO2", "sO", "R1", "POM", { lt: 0.45 }); e("fCH4", "sCH3", "SH", "C–H 活化 ⇌ 甲烷化", { both: true, lt: 0.4 });
  e("sO", "sOH", "SH", "+H*", { perp: 0.25 }); e("sOH", "sH", "SH", "");
  e("fCO2", "eCO2m", "R3", "+e^{−}（电还原）", { lt: 0.45, perp: 0.25 }); e("eCO2m", "mHCOOH", "R3", "+H^{+} + e^{−}", { ly: -0.27 }); e("mHCOOH", "sHCOO", "R3", "", { both: true });
  e("sCO2", "sHCOO", "R3", ""); e("sHCOO", "sH2COO", "R3", ""); e("sH2COO", "sH2COOH", "R3", ""); e("sH2COOH", "sH2CO", "R3", ""); e("sH2CO", "sCH3O", "R3", ""); e("sCH3O", "mMeOH", "R3", "");
  e("sCO2", "sCOOH", "R3", "+H*", { perp: -0.22 }); e("sCOOH", "sCO", "SH", "−OH*（RWGS）", { lt: 0.42, perp: -0.2, lw: 1.6 }); e("hubCO", "sCO", "SH", "", { both: true });
  e("sCO", "sHCO", "R3", "+H*", { perp: -0.2 }); e("sHCO", "sH2CO", "R3", "+H*", { perp: -0.2 }); e("sCO", "sCOH", "R3", "+H*", { perp: 0.2 });
  e("sCOH", "sC", "SH", "H 助 C–O 解离", { lt: 0.55, perp: -0.2, lw: 1.6 }); e("sCO", "sC", "SH", "C* + O* ⇌ CO*", { both: true, perp: -0.75 });
  e("sC", "sCH", "R3", "", { both: true }); e("sCH", "sCH2", "R3", "", { both: true }); e("sCH2", "sCH3", "R3", "", { both: true });
  e("sCH2", "sCnH", "R2", "CH_2* 插入 → 链增长（F–T / FTO）", { pts: [[5.3, 7.67], [5.3, 8.08], [10.0, 8.08], [10.0, 9.15], [12.4, 9.15]], lpos: [6.55, 8.08] });
  e("sCH3", "sCH3CO", "R2", "+CO* 插入（还原偶联 / 同系化）", { pts: [[4.0, 7.67], [4.0, 8.24], [10.45, 8.24]], lpos: [8.7, 8.24] });
  e("mMeOH", "zOCH3", "R4", "脱水", { perp: 0.21 }); e("zOCH3", "mDME", "R4", ""); e("zOCH3", "iPOOL", "R4", "甲基化", { perp: 0.22, lt: 0.45 });
  e("mDME", "iPOOL", "R4", "MTH", { perp: -0.22 }); e("sCO", "zKET", "R4", "OEO：CO* + CH_3* → 乙烯酮", { dash: "dash", lt: 0.55, perp: -0.2, lw: 2.75 });
  e("zKET", "mC24", "R4", "乙烯酮 → 乙烯", { pts: [[11.3, 4.12], [11.3, 4.65], [17.45, 4.65], [17.45, 4.12]], lpos: [12.9, 4.65], dash: "dash" });
  e("iPOOL", "mC24", "R4", "侧链消除 / 脱烷基", { perp: 0.2 }); e("iPOOL", "mARO", "R4", "氢转移 / 脱氢", { perp: 0.22 });
  e("mC24", "mC5O", "R4", "甲基化 / 齐聚 ⇄ 裂解", { both: true, perp: 0.3, lt: 0.35 }); e("mC5O", "iPOOL", "R4", "环化", { perp: -0.2, lt: 0.4 });
  e("mC5O", "mGAS", "R3", "+H（氧化物加氢）", { ly: -0.28 }); e("mEtOH", "zOEt", "R4", "脱水", { perp: 0.22 }); e("zOEt", "mC24", "R4", "");
  e("sCH3CO", "mEtOH", "R3", "+4H*（经 CH_3CHO）", { perp: -1.4 }); e("sCH3CO", "sCnH", "R2", "链增长", { perp: 0.2 });
  e("sCnH", "mAOle", "R2", "β-H 消除", { perp: -0.25 }); e("sCnH", "mNPar", "R2", "+H* 终止", { ly: -0.25 }); e("sCnH", "mWax", "R2", "链增长", { perp: -0.2, lt: 0.6 });
  e("sCnH", "mHA", "R2", "CO* 插入", { perp: 0.22 }); e("mAOle", "iCARB", "R5", "质子化", { perp: 0.2 }); e("mWax", "iCARB", "R5", "金属脱氢 → 质子化", { perp: -0.22, lt: 0.45 });
  e("iCARB", "mJET", "R5", "β-断裂 / 异构 / 加氢", { perp: 0.25 });
  e("iCARB", "mC24", "R5", "裂解 / 低聚", { pts: [[19.08, 8.3], [19.35, 8.3], [19.35, 3.9], [18.68, 3.9]], lpos: [18.85, 7.55] });
  e("fCH4", "rCH3", "R6", "高温均裂 / +O_2", { pts: [[0.95, 7.77], [0.95, 10.5], [2.6, 10.5]], lpos: [1.75, 10.25] });
  e("tCO2H", "rHOCO", "R6", ""); e("tCO2H", "rHCOO", "R6", ""); e("rHOCO", "rHCO", "R6", "−OH• → CO；+H•", { ly: 0.27 }); e("rHCOO", "rHCOOH", "R6", "");
  e("rHCO", "rCH2O", "R6", "+H•", { ly: 0.27 }); e("rCH2O", "rCH3O", "R6", "+H•", { ly: 0.27 }); e("rCH2O", "rCH2OH", "R6", "+H•", { perp: -0.22 });
  e("rCH3O", "rMeOH", "R6", "+H•", { perp: 0.2 }); e("rCH2OH", "rMeOH", "R6", ""); e("rCH3O", "rCH2OH", "R6", "", { both: true });
  e("rCH3", "rC2", "R6", "OCM", { ly: -0.25 }); e("rCH3", "rCH3OO", "R6", "+O_2", { pts: [[3.1, 10.71], [3.1, 10.88], [6.45, 10.88], [6.45, 10.71]], lpos: [5.75, 10.88] });
  e("rCH3OO", "rCH2O", "R6", "−OH•", { perp: 0.25 }); e("rH", "rHO2", "R6", "+O_2", { ly: -0.25 }); e("rHO2", "rOH", "R6", "+H•", { ly: -0.25 });
  e("rOH", "rH", "R6", "+H_2 → H_2O", { pts: [[9.35, 10.71], [9.35, 10.88], [7.6, 10.88], [7.6, 10.71]], lpos: [8.5, 10.88] });
  OV.zones.push({ x: 2.45, y: 2.32, w: 7.35, h: 5.68, color: C.accent1, title: "C_1 表面中间体核心（金属 / 金属–氧化物界面）→ 详见 A、B", tx: 4.3, tw: 5.35, talign: "right", fs: 12 },
                { x: 2.45, y: 8.42, w: 7.35, h: 2.63, color: C.accent6, title: "⑥ 无需催化剂：气相 / 均相自由基（POM、OCM、高温）→ 详见 C", dash: "dash", fill: CATAL.GAS.fill, fs: 12 },
                { x: 10.1, y: 1.15, w: 9.5, h: 5.5, color: C.accent4, title: "④ 分子筛酸位：表面甲氧基 / 烃池 / 碳正离子（OX-ZEO）→ 详见 E", fs: 12, ty: 1.22 },
                { x: 10.1, y: 6.88, w: 9.5, h: 4.17, color: C.accent2, title: "② 金属位 C–C 链增长 → 详见 D · ⑤ 酸位烃类升级 → 详见 F", fs: 12, ty: 6.92, tx: 11.75, tw: 7.7 });
  OV.blocks = [["EC", 2.4, 1.32, 3.0, 0.86, [4.25, 1.33], "电催化阴极"],
    ["OX", 2.52, 2.7, 7.21, 2.05, [4.25, 3.72], "还原性氧化物 / Cu–ZnO（甲醇 · RWGS）"],
    ["MET", 2.52, 4.78, 7.21, 3.17, [6.15, 7.72], "过渡金属 Ni / Co / Ru / Fe（CO 解离 · 甲烷化 · 重整）"],
    ["BA", 10.18, 1.55, 9.34, 5.02, [14.2, 6.33], "Brønsted 酸分子筛（H-ZSM-5 / SAPO-34）"],
    ["OX", 10.4, 3.6, 1.8, 0.6, [10.42, 4.22], "氧化物（OEO）"], ["LA", 14.6, 2.05, 2.1, 0.6, [16.72, 2.1], "Ga / Zn（脱氢）"],
    ["OX", 14.85, 5.05, 4.0, 0.6, [16.45, 5.67], "氧化物加氢"],
    ["CARB", 10.18, 7.25, 6.62, 3.73, [13.3, 10.74], "碳化物 Fe₅C₂ / Co₂C"], ["DUAL", 10.2, 7.2, 1.65, 1.35, [10.22, 8.58], "CO 插入双位"],
    ["DUAL", 11.3, 10.12, 2.0, 0.56, null], ["BIF", 16.85, 7.25, 2.67, 3.73, [16.9, 10.72], "双功能 金属–酸"]];
  OV.texts.push([0.45, 1.62, 1.9, 0.28, "电化学 / 光化学 → A", { fs: 11, color: C.accent3, bold: true }],
                [4.25, 3.42, 3.95, 0.26, "甲酸盐路线：每步 +H*（H_2COOH* → H_2CO* 为 −OH*）", { fs: 10, color: C.accent3, align: "center" }]);
  OV.cats = {};
  const NC = { fCO2: "R1 R3", fH2: "R2 R3", fH2O: "R1", fO2: "R1 R6", fCH4: "R1 R3 R6", eCO2m: "R3 R6", mHCOOH: "R3",
    sCO2: "R1 R3", sHCOO: "R3", sH2COO: "R3", sH2COOH: "R3", sH2CO: "R3", sCH3O: "R3", sH: "R1 R2 R3", sOH: "R1 R3", sO: "R1 R3",
    sCOOH: "R3", hubCO: "R1 R2 R3 R4", sHCO: "R2 R3", sCO: "R1 R2 R3 R4", sCOH: "R2 R3", sC: "R1 R2 R3", sCH: "R1 R2 R3", sCH2: "R2 R3", sCH3: "R1 R2 R3",
    rCH3: "R1 R6", mMeOH: "R3 R4", zOCH3: "R4", mDME: "R4", mARO: "R4 R5", mGAS: "R4 R5", zKET: "R4", iPOOL: "R4", mC24: "R2 R4 R5",
    zOEt: "R4", mC5O: "R4 R5", mEtOH: "R2 R3 R4", sCH3CO: "R2 R3", sCnH: "R2", mAOle: "R2 R5", mNPar: "R2", mWax: "R2 R5", mHA: "R2 R3", iCARB: "R4 R5", mJET: "R5" };
  for (const [k, v] of Object.entries(NC)) OV.cats[k] = v.split(" ");
}

// ================================================================== SLIDES
const want = (k) => !ONLY || ONLY.includes(k);
// 1. overview
if (want("1")) {
  const s = pres.addSlide({ masterName: "POSTER" });
  s.addText("CO₂ 加氢及 C1 催化反应网络 · 总览（模块 A–F 在后续页展开）", { placeholder: "title" });
  const T = { ox: 0.8, oy: -0.54, s: 1.12 };
  (OV.zones).forEach((z) => drawZone(s, T, z));
  drawBlocks(s, T, OV.blocks, "OV");
  OV.texts.forEach(([x, y, w, h, t, o]) => txt(s, T, x, y, w, h, t, o));
  const [N, E] = buildMap(OV);
  E.forEach((ed, i) => drawEdge(s, T, N, ed, i, "OV"));
  Object.values(N).forEach((nd) => drawNode(s, T, nd, "OV", OV.cats[nd.id]));
  drawChips(s, T, OV.blocks, "OV", 9);
  catalystLegend(s, 12.72);
  legend(s, 13.1);
  s.addNotes("总览：去掉羰基化后的完整网络骨架。橙框（C1 核心）在模块 A（氧化物、甲醇合成）与 B（金属、CO 活化 / 甲烷化 / 重整）展开；红色虚线框（自由基）在模块 C 展开；绿框（分子筛）在模块 E 展开；红框（链增长与升级）在模块 D 与 F 展开。");
}
// 1b. catalyst classification table
if (want("1b")) {
  const s = pres.addSlide({ masterName: "POSTER" });
  s.addText("催化剂形式分类：网络图中背景色块的含义", { placeholder: "title" });
  const R = [
    ["OX", "还原性氧化物 · 氧空位", "In₂O₃、ZnZrOₓ、GaZrOₓ、ZnCrOₓ", "Vₒ 吸附活化 CO₂；M–O 对异裂 H₂；H 溢流加氢烯烃", "甲酸盐路线 → 甲醇；氧化还原 RWGS；OEO 乙烯酮；烯烃加氢（A、E）"],
    ["CU", "Cu 基 / Cu–氧化物界面", "Cu/ZnO/Al₂O₃、Cu/ZrO₂、Pd–Cu", "Cu⁰–ZnOₓ 界面；CO / 甲酸盐不解离加氢", "甲醇合成、缔合 RWGS、CO 加氢（A、B）"],
    ["MET", "过渡金属", "Ni、Co、Ru、Rh、Fe⁰", "C–O 解离、C–H 活化、H₂ 解离", "甲烷化、SMR / DRM / POM 表面步骤、积碳（B）"],
    ["CARB", "金属碳化物（+K / Na）", "χ-Fe₅C₂、θ-Fe₃C、Co₂C", "生成并插入 CHₓ；碱金属抑制终止加氢、促进 β-H 消除", "费托 / FTO 链增长、α-烯烃、支化（D）"],
    ["DUAL", "CO 插入双位", "Co⁰–Co₂C、CuCo、Cu–Fe、Rh", "解离位供 CHₓ*，非解离位供 CO*", "烷基 + CO* → 酰基 → C₂₊ 醇（D）"],
    ["SUP", "载体 / 金属–载体界面", "CeO₂、ZrO₂、MgO、La₂O₃", "碱性位与氧空位活化 CO₂；参与 WGS", "DRM 中 CO₂ 活化、WGS（B）"],
    ["BA", "Brønsted 酸分子筛", "H-ZSM-5、SAPO-34、H-β、MOR", "质子化、表面甲氧基、碳正离子、孔道择形", "甲醇 → DME → 烯烃 / 芳烃双循环、裂解、齐聚（E、F）"],
    ["LA", "Lewis 酸 / 金属改性分子筛", "Ga-ZSM-5、Zn-ZSM-5", "脱氢放出 H₂，绕过氢转移的芳烃上限", "脱氢芳构化（E）"],
    ["BIF", "双功能 金属–酸", "Pt、Pd、Ni + 分子筛 / 无定形硅铝", "金属脱氢 / 加氢 + 酸位质子化、PCP 异构、β-断裂", "蜡与正构烷加氢异构 / 裂化 → 航煤（F）"],
    ["EC", "电催化阴极", "Cu、Ag、Au、Sn、Bi", "单电子还原 CO₂ → CO₂•⁻", "HCOOH / CO（A）"],
    ["GAS", "无需催化剂", "高温、等离子体、光解", "自由基链引发与传递", "OCM、POM、气相自由基加氢链（C）"],
  ];
  const th = ["色块", "催化剂形式", "代表材料", "活性位 / 作用", "对应反应路径（模块）"].map((t) => ({ text: t, options: { bold: true, color: C.background1, fill: { color: C.text2 } } }));
  s.addTable([th, ...R.map(([k, a, b, c2, d]) => [
    { text: CATAL[k].s, options: { fill: { color: CATAL[k].fill }, color: CATAL[k].c, bold: true } },
    { text: a, options: { bold: true } }, { text: b }, { text: c2 }, { text: d }])],
    { x: 0.3, y: 0.95, w: 23.4, colW: [2.6, 3.6, 4.6, 6.4, 6.2], rowH: 0.92, fontSize: 15, color: C.text1,
      border: { type: "solid", pt: 0.5, color: "D6DADC" }, margin: [0, 0.12, 0, 0.12], valign: "middle", objectName: "catalyst table" });
  s.addNotes("色块表示该段反应路径需要的催化剂形式；同一个物种可能位于两种色块的交界（例如 CO* 同时在过渡金属与碳化物上反应）。");
}
// 2. all modules on one slide
const POS = { A: [0.3, 0.8], B: [7.15, 0.8], C: [14.0, 0.8], D: [0.3, 6.72], E: [7.55, 6.72], F: [17.65, 6.72] };
const SP = 0.95;
if (want("2")) {
  const s = pres.addSlide({ masterName: "POSTER" });
  s.addText("细致反应网络总图：六个模块全部展开（可放大查看）", { placeholder: "title" });
  for (const M of MODS) drawModule(s, M, { ox: POS[M.id][0], oy: POS[M.id][1], s: SP });
  catalystLegend(s, 12.72);
  legend(s, 13.1);
  s.addNotes("六个模块的完整细节放在同一页，字号较小，适合放大或打印查看；每个模块在后面单独放大。");
}
// 3-8. one module per slide
MODS.forEach((M, i) => {
  if (!want(String(3 + i))) return;
  const s = pres.addSlide({ masterName: "POSTER" });
  s.addText(runs(M.title.replace(/^([A-F])  /, "模块 $1 · "), { fontSize: 24, bold: true, color: C.text1 }), { placeholder: "title" });
  const availW = 23.4, availH = 11.78;
  let sc = Math.min(availW / M.W, availH / M.H);
  const panel = M.notes && M.notes.length && availW - M.W * sc >= 3.6;
  if (!panel && M.notes && M.notes.length) sc = Math.min(sc, (availW - 0) / M.W);
  const ox = 0.3, oy = 0.85;
  drawModule(s, M, { ox, oy, s: sc }, { frame: true });
  if (panel) {
    const px = ox + M.W * sc + 0.35, pw = SW - 0.3 - px;
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: px, y: oy, w: pw, h: M.H * sc, rectRadius: 0.12, fill: { color: C.background2 },
      line: { color: C.background2 }, objectName: "notes panel" });
    s.addText([{ text: "要点", options: { bold: true, fontSize: 18, color: M.color, breakLine: true } },
               ...M.notes.map((t, k) => ({ text: t, options: { bullet: true, fontSize: 15, color: C.text1, breakLine: k < M.notes.length - 1, paraSpaceAfter: 10 } }))],
      { x: px + 0.25, y: oy + 0.25, w: pw - 0.5, h: M.H * sc - 0.5, valign: "top", margin: 0, isTextBox: true, objectName: "notes text" });
  }
  catalystLegend(s, 12.72);
  legend(s, 13.1);
  s.addNotes((M.notes && M.notes.length ? M.notes.join("\n") : "") || M.title);
});
// 9. radicals table
const TBL_T = { ox: 0, oy: 0, s: 1 };
if (want("9")) {
  const s = pres.addSlide({ masterName: "POSTER" });
  s.addText("网络中的全部自由基：生成、消耗与对应的表面物种", { placeholder: "title" });
  const RADS = [
    ["H•", "H₂ 均裂；OH• + H₂ → H₂O + H•", "加成到 CO₂ / CO / CH₂O / C₂H₄；+O₂ → HO₂• 或 OH• + O•", "H*", "通用"],
    ["OH•", "HOCO• → CO + OH•；HO₂• + H•；CH₃OO• 分解；H₂O₂ → 2OH•", "+H₂ → H₂O + H•；抽氢 RH → R• + H₂O", "OH*", "通用"],
    ["O•", "O₂ 高温解离；H• + O₂ → OH• + O•（支化）", "+H₂ / +CH₄ → OH• + H• / •CH₃", "O*", "POM / 高温"],
    ["HO₂•", "H• + O₂ (+M)", "+H• → 2OH•；+HO₂• / RH → H₂O₂", "—", "POM"],
    ["CO₂•⁻", "CO₂ + e⁻", "+H⁺ + e⁻ → HCOO⁻ / HCOOH 或 CO + OH⁻", "CO₂(δ−)*", "电化学"],
    ["HOCO•", "CO₂ + H•；CO + OH•", "→ CO + OH•；+H• → HCOOH；⇌ HCOO•", "COOH*", "RWGS"],
    ["HCOO•", "CO₂ + H•", "+H• → HCOOH；⇌ HOCO•", "HCOO*", "甲酸盐路线"],
    ["HOC•HOH", "HCOOH + H•", "→ CH₂O + OH•", "H₂COOH*", "甲酸盐路线"],
    ["HCO•", "CO + H•；CH₂O 被抽氢", "+H• → CH₂O；→ CO + H•", "HCO*", "CO 加氢"],
    ["CH₃O•", "CH₂O + H•", "+H• → CH₃OH；⇌ •CH₂OH", "CH₃O*", "甲醇"],
    ["•CH₂OH", "CH₂O + H•；CH₃OH 被抽氢", "+H• → CH₃OH；⇌ CH₃O•", "CH₂OH*", "甲醇"],
    ["•CH₃", "CH₄ 抽氢或均裂；CH₃OH → •CH₃ + OH•", "2•CH₃ → C₂H₆；+O₂ → CH₃OO•；+H• → CH₄", "CH₃*", "OCM / POM / 甲烷化"],
    ["CH₃OO•", "•CH₃ + O₂", "→ CH₂O + OH•", "—", "POM"],
    ["C₂H₅•", "C₂H₆ 被抽氢；C₂H₄ + H•", "β-断裂 → C₂H₄ + H•；歧化", "C₂H₅*", "OCM"],
    ["C₂H₃•", "C₂H₄ 被抽氢", "→ C₂H₂ + H•", "CHCH₂*", "OCM / 高温"],
  ];
  const th = ["自由基", "生成", "消耗", "对应表面物种", "场景"].map((t) => ({ text: t, options: { bold: true, color: C.background1, fill: { color: C.text2 } } }));
  s.addTable([th, ...RADS.map((r, i) => r.map((t, j) => ({ text: t, options: { fill: { color: i % 2 ? C.background1 : C.background2 }, bold: j === 0,
    color: j === 0 ? C.accent6 : C.text1 } })))], { x: 0.3, y: 0.95, w: 23.4, colW: [2.4, 7.4, 8.0, 2.8, 2.8], rowH: 0.7, fontSize: 15,
    border: { type: "none" }, margin: [0, 0.12, 0, 0.12], valign: "middle", objectName: "radical table" });
  s.addText("C1 显式枚举：10 个自由基、141 个反应（data/network_full_C1.csv）；C≤12 全空间约 3×10¹¹ 个自由基。", { x: 0.3, y: 12.7, w: 20, h: 0.35,
    fontSize: 13, color: C.text2, margin: 0, isTextBox: true });
}
// 10. counts per module
if (want("10")) {
  const s = pres.addSlide({ masterName: "POSTER" });
  s.addText("细致网络的规模：各模块的物种与反应边", { placeholder: "title" });
  const rows = [["模块", "内容", "分子", "表面物种", "自由基", "离子", "反应边"]];
  for (const M of MODS) {
    const [N, E] = buildMap(M);
    const c = { mol: 0, surf: 0, rad: 0, ion: 0 };
    Object.values(N).forEach((x) => { if (c[x.type] !== undefined) c[x.type]++; });
    rows.push([M.id, M.title.replace(/^[A-F]  /, "").replace(/_\{([^}]*)\}|_(\d+|[a-z])/g, "$1$2").replace(/\^\{([^}]*)\}/g, "$1").split("：")[0].split("（")[0],
               String(c.mol), String(c.surf), String(c.rad), String(c.ion), String(E.length)]);
  }
  const tot = rows.slice(1).reduce((a, r) => a.map((v, i) => (i >= 2 ? v + Number(r[i]) : v)), ["合计", "", 0, 0, 0, 0, 0]);
  rows.push(tot.map(String));
  s.addTable(rows.map((r, i) => r.map((t, j) => ({ text: t, options: { bold: i === 0 || j === 0 || i === rows.length - 1,
    color: i === 0 ? C.background1 : C.text1, fill: { color: i === 0 ? C.text2 : i % 2 ? C.background1 : C.background2 }, align: j >= 2 ? "right" : "left" } }))),
    { x: 0.3, y: 1.0, w: 15.5, colW: [1.0, 8.0, 1.2, 1.4, 1.2, 1.2, 1.5], rowH: 0.75, fontSize: 16, border: { type: "none" }, margin: [0, 0.15, 0, 0.15],
      valign: "middle", objectName: "module counts" });
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 16.3, y: 1.0, w: 7.4, h: 6.75, rectRadius: 0.12, fill: { color: C.background2 }, line: { color: C.background2 }, objectName: "scope card" });
  s.addText([
    { text: "这张网络是什么", options: { bold: true, fontSize: 20, color: C.text1, breakLine: true } },
    { text: "机理骨架：每条边是一个基元步骤或一类等价步骤（如每步 +CH₃⁺ 甲基化），C₁–C₄ / C₂–C₇ 显式写出，更长链用通式。", options: { bullet: true, fontSize: 15, color: C.text1, breakLine: true, paraSpaceAfter: 8 } },
    { text: "逐条枚举的全部基元反应在 data/network_full_C2.csv（C≤2）与 data/network_reduced_*.csv（按产物剪枝）。", options: { bullet: true, fontSize: 15, color: C.text1, breakLine: true, paraSpaceAfter: 8 } },
    { text: "已去掉羰基化（亚硝酸甲酯 → 乙二醇）与 C–N 偶联（氨 / 尿素）部分。", options: { bullet: true, fontSize: 15, color: C.text1 } },
  ], { x: 16.6, y: 1.25, w: 6.8, h: 6.3, valign: "top", margin: 0, isTextBox: true, objectName: "scope text" });
}

(async () => { await pres.writeFile({ fileName: OUT }); await applyTheme(OUT, THEME); console.log("wrote", OUT); })();
