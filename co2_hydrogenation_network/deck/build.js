// CO2 -> jet fuel OX-ZEO reaction network deck
const pptxgen = require("pptxgenjs");
const fs = require("fs");
const { applyTheme } = require("/root/.claude/skills/synced/62cc770a-8c2c-45f8-8ae0-ca258c8631d4_456331ad-9fc1-4922-9170-b4c06151a9b5/pptx/scripts/apply_theme.js");

const D = JSON.parse(fs.readFileSync(__dirname + "/deckdata.json", "utf8"));
const OUT = process.argv[2] || "deck.pptx";
const [ARC_W, ARC_H] = JSON.parse(fs.readFileSync(__dirname + "/arcs_size.json", "utf8"));

const THEME = {
  name: "CO2 Network",
  headFontFace: "Microsoft YaHei",
  bodyFontFace: "Microsoft YaHei",
  colors: {
    dk1: "14232B", lt1: "FFFFFF", dk2: "0F3D3E", lt2: "EEF3F2",
    accent1: "2A78D6", accent2: "EB6834", accent3: "1BAF7A", accent4: "D99300",
    accent5: "6250D6", accent6: "5F6B6D", hlink: "2A78D6", folHlink: "6250D6",
  },
};
const HEX = THEME.colors;
const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE"; // 13.33 x 7.5
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
pres.title = "CO2 加氢制航煤的反应网络";
const C = pres.SchemeColor;
const CLS = { par: C.accent1, ole: C.accent2, nap: C.accent3, aro: C.accent4 };
const FOOT = "CO₂ 加氢制航煤 · OX-ZEO 反应网络";

// ------------------------------------------------------------------ layouts
pres.defineSlideMaster({
  title: "TITLE_DARK",
  background: { color: C.text2 },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.7, y: 2.2, w: 7.6, h: 1.7, fontSize: 40, bold: true,
        color: C.background1, valign: "bottom", align: "left", margin: 0 }, text: "" } },
    { placeholder: { options: { name: "body", type: "body", x: 0.7, y: 4.15, w: 7.4, h: 1.6, fontSize: 16,
        color: C.background2, valign: "top", margin: 0 }, text: "" } },
  ],
});
pres.defineSlideMaster({
  title: "CONTENT",
  background: { color: C.background1 },
  margin: [0.5, 0.5, 0.6, 0.5],
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.55, y: 0.32, w: 12.2, h: 0.85, fontSize: 28, bold: true,
        color: C.text2, valign: "middle", align: "left", margin: 0 }, text: "" } },
    { text: { text: FOOT, options: { x: 0.55, y: 7.0, w: 6, h: 0.3, fontSize: 10, color: C.accent6, margin: 0 } } },
  ],
  slideNumber: { x: 12.35, y: 7.0, w: 0.5, h: 0.3, fontSize: 10, color: C.accent6, align: "right" },
});
pres.defineSlideMaster({
  title: "CLOSING_DARK",
  background: { color: C.text2 },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.7, y: 0.5, w: 12, h: 0.9, fontSize: 32, bold: true,
        color: C.background1, valign: "middle", align: "left", margin: 0 }, text: "" } },
  ],
});

// ------------------------------------------------------------------ helpers
// "CO_2 + H_2" -> runs with subscripts (digits after "_")
function runs(s, opt = {}) {
  const out = [];
  const re = /_(\d+)/g;
  let last = 0, m;
  while ((m = re.exec(s))) {
    if (m.index > last) out.push({ text: s.slice(last, m.index), options: { ...opt } });
    out.push({ text: m[1], options: { ...opt, subscript: true } });
    last = re.lastIndex;
  }
  if (last < s.length) out.push({ text: s.slice(last), options: { ...opt } });
  return out;
}
function clip(cx, cy, w, h, dx, dy, shape) {
  // point on node boundary in direction (dx,dy) from centre
  if (shape === "ellipse") {
    const a = w / 2, b = h / 2;
    const t = 1 / Math.sqrt((dx * dx) / (a * a) + (dy * dy) / (b * b));
    return [cx + dx * t, cy + dy * t];
  }
  const tx = dx !== 0 ? (w / 2) / Math.abs(dx) : Infinity;
  const ty = dy !== 0 ? (h / 2) / Math.abs(dy) : Infinity;
  const t = Math.min(tx, ty);
  return [cx + dx * t, cy + dy * t];
}
function line(slide, x1, y1, x2, y2, o) {
  const opts = {
    x: Math.min(x1, x2), y: Math.min(y1, y2), w: Math.max(Math.abs(x2 - x1), 0.001), h: Math.max(Math.abs(y2 - y1), 0.001),
    flipH: x2 < x1, flipV: y2 < y1,
    line: { color: o.color, width: o.width || 1.5, dashType: o.dash || "solid",
            endArrowType: o.head === false ? undefined : "triangle",
            beginArrowType: o.both ? "triangle" : undefined },
    objectName: o.name,
  };
  slide.addShape(pres.shapes.LINE, opts);
}
function edge(slide, A, B, o) {
  const dx = B.cx - A.cx, dy = B.cy - A.cy, L = Math.hypot(dx, dy), ux = dx / L, uy = dy / L;
  const gap = o.gap ?? 0.04;
  let [x1, y1] = clip(A.cx, A.cy, A.w, A.h, ux, uy, A.shape);
  let [x2, y2] = clip(B.cx, B.cy, B.w, B.h, -ux, -uy, B.shape);
  x1 += ux * gap; y1 += uy * gap; x2 -= ux * gap; y2 -= uy * gap;
  line(slide, x1, y1, x2, y2, o);
  if (o.label) {
    const t = o.lt ?? 0.5;
    let nx = -uy, ny = ux;                       // unit normal; point it up (or left when vertical)
    if (ny > 0 || (Math.abs(ny) < 1e-9 && nx > 0)) { nx = -nx; ny = -ny; }
    const pp = o.perp || 0;
    const mx = x1 + (x2 - x1) * t + (o.lx || 0) + nx * pp, my = y1 + (y2 - y1) * t + (o.ly || 0) + ny * pp;
    const lw = o.lw || 0.9, lh = o.lh || 0.28;
    slide.addText(runs(o.label, { fontSize: o.lfs || 11, color: o.lcolor || o.color, bold: !!o.lbold }), {
      x: mx - lw / 2, y: my - lh / 2, w: lw, h: lh, align: "center", valign: "middle", margin: 0,
      ...(o.lfill === null ? {} : { fill: { color: o.lfill || C.background1 } }), isTextBox: true, objectName: (o.name || "edge") + " label",
    });
  }
}
function node(slide, n) {
  const shape = n.shape === "ellipse" ? pres.shapes.OVAL : pres.shapes.ROUNDED_RECTANGLE;
  const opt = {
    x: n.cx - n.w / 2, y: n.cy - n.h / 2, w: n.w, h: n.h, shape,
    fill: n.fill ? { color: n.fill } : { color: C.background1 },
    line: { color: n.line || n.fill || C.accent6, width: n.lineW || 1.25, dashType: n.dash || "solid" },
    align: "center", valign: "middle", margin: 0.03, objectName: n.name,
  };
  if (shape === pres.shapes.ROUNDED_RECTANGLE) opt.rectRadius = n.r ?? 0.1;
  const text = Array.isArray(n.text) ? n.text : runs(n.text || "", { fontSize: n.fs || 14, bold: n.bold ?? true, color: n.tc || C.background1 });
  slide.addText(text, opt);
  return n;
}
const W = (f, max, lo = 1, hi = 9) => lo + (hi - lo) * Math.sqrt(f / max);
const fmt = (v) => (v >= 0.2 ? v.toFixed(2) : v.toFixed(3));

// ================================================================== slides
// 1. title
pres.addSection({ title: "总览" });
let s = pres.addSlide({ masterName: "TITLE_DARK", sectionTitle: "总览" });
s.addText("CO₂ 加氢制航煤的反应网络", { placeholder: "title" });
s.addText([
  { text: "氧化物/分子筛（OX-ZEO）甲醇中间体路线", options: { breakLine: true } },
  { text: "集总网络：72 个物种 · 299 个反应 · 11 个反应族", options: { breakLine: true } },
  { text: "通量按郭新闻课题组 GaZrOx/Ga-ZSM-5 数据校准" },
], { placeholder: "body" });
// motif: small network on the right
{
  const pts = [[9.4, 2.0, CLS.ole], [10.6, 1.5, CLS.ole], [11.8, 2.2, CLS.aro], [10.2, 3.1, C.accent5], [11.4, 3.5, CLS.nap],
               [12.4, 3.2, CLS.aro], [9.6, 4.4, CLS.par], [10.9, 4.8, CLS.ole], [12.1, 4.7, CLS.par], [10.3, 5.9, CLS.aro], [11.6, 6.0, CLS.ole]];
  const links = [[0, 1], [1, 2], [0, 3], [1, 3], [3, 4], [4, 2], [4, 5], [2, 5], [3, 6], [3, 7], [7, 4], [7, 8], [6, 9], [7, 9], [9, 10], [8, 10], [5, 8]];
  links.forEach(([a, b], i) => line(s, pts[a][0], pts[a][1], pts[b][0], pts[b][1],
    { color: C.background2, width: 1.25, head: false, name: `motif link ${i + 1}` }));
  pts.forEach(([x, y, col], i) => s.addShape(pres.shapes.OVAL, { x: x - 0.22, y: y - 0.22, w: 0.44, h: 0.44,
    fill: { color: col }, line: { color: C.background1, width: 2 }, objectName: `motif node ${i + 1}` }));
}
s.addNotes("这份文件把目前的 OX-ZEO 制航煤反应网络画成网络图。网络按“类别 × 碳数”集总：烯烃、烷烃、环状中间体、芳烃各按碳数一个节点，外加 CO2、H2、CO、H2O、甲醇、DME。所有通量数字来自按郭组 GaZrOx/Ga-ZSM-5（Adv. Sci. 2026）摘要数据校准的模型，320 °C、3 MPa。");

// 2. overview flow network
s = pres.addSlide({ masterName: "CONTENT", sectionTitle: "总览" });
s.addText("网络总览：氧化物生成甲醇，分子筛把甲醇变成烃", { placeholder: "title" });
s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 0.55, y: 1.45, w: 3.75, h: 5.15, rectRadius: 0.12, fill: { color: C.background2 },
  line: { color: C.background2 }, objectName: "oxide zone" });
s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 5.95, y: 1.45, w: 6.85, h: 5.15, rectRadius: 0.12, fill: { color: C.background2 },
  line: { color: C.background2 }, objectName: "zeolite zone" });
s.addText("氧化物位点（GaZrOx / ZnZrOx / In₂O₃）", { x: 0.75, y: 1.55, w: 3.4, h: 0.35, fontSize: 12, bold: true, color: C.accent6, margin: 0, isTextBox: true });
s.addText("分子筛酸位（ZSM-5）", { x: 6.15, y: 1.55, w: 4, h: 0.35, fontSize: 12, bold: true, color: C.accent6, margin: 0, isTextBox: true });
const ox = D.oxide.G1_GaZrOx_GaZSM5, ff = D.family_flux.G1_GaZrOx_GaZSM5;
const nFeed = node(s, { cx: 1.55, cy: 3.95, w: 1.6, h: 0.75, fill: C.accent5, text: "CO_2 + H_2", name: "feed CO2 H2" });
const nCO = node(s, { cx: 3.35, cy: 5.75, w: 1.2, h: 0.6, fill: C.background1, line: C.accent6, tc: C.accent6, text: "CO", name: "CO" });
const nMe = node(s, { cx: 5.12, cy: 3.95, w: 1.45, h: 0.95, fill: C.accent5, name: "methanol DME",
  text: [...runs("CH_3OH", { fontSize: 14, bold: true, color: C.background1, breakLine: false }),
         { text: "", options: { breakLine: true } }, { text: "⇌ DME", options: { fontSize: 12, color: C.background1 } }] });
const nCH4 = node(s, { cx: 5.12, cy: 5.75, w: 1.0, h: 0.5, fill: C.background1, line: C.accent6, tc: C.accent6, text: "CH_4", fs: 12, name: "CH4" });
const nOle = node(s, { cx: 7.6, cy: 3.95, w: 2.25, h: 1.15, fill: CLS.ole, name: "olefin pool",
  text: [{ text: "烯烃池 C2=–C11=", options: { fontSize: 14, bold: true, color: C.background1, breakLine: true } },
         { text: `裂解 ${fmt(ff.Z5_crack)} · 齐聚 ${fmt(ff.Z4_olig)}`, options: { fontSize: 11, color: C.background1 } }] });
const nPar = node(s, { cx: 7.55, cy: 2.3, w: 2.0, h: 0.62, fill: CLS.par, text: "烷烃 C2–C10", name: "paraffins" });
const nNap = node(s, { cx: 9.75, cy: 5.65, w: 1.8, h: 0.62, fill: CLS.nap, text: "环状中间体", name: "cyclic intermediates" });
const nAro = node(s, { cx: 11.65, cy: 3.95, w: 1.85, h: 1.15, fill: CLS.aro, name: "aromatics",
  text: [{ text: "芳烃 C6–C11", options: { fontSize: 14, bold: true, color: C.background1, breakLine: true } },
         { text: `+CH₃OH 甲基化 ${fmt(ff.Z8_ameth)}`, options: { fontSize: 11, color: C.background1 } }] });
const mx = 0.88;
edge(s, nFeed, nMe, { color: C.accent5, width: W(ox.MeOH_synthesis, mx), label: `甲醇合成 ${fmt(ox.MeOH_synthesis)}`, lw: 1.45, ly: -0.3, name: "OX1" });
edge(s, nFeed, nCO, { color: C.accent6, width: W(ox.RWGS, mx), label: `RWGS ${fmt(ox.RWGS)}`, lw: 1.05, lx: -0.55, ly: 0.05, name: "OX2" });
edge(s, nMe, nOle, { color: C.accent5, width: W(ff.Z3_meth, mx), label: `甲基化 ${fmt(ff.Z3_meth)}`, lw: 1.2, ly: -0.62, lx: -0.1, name: "Z3" });
edge(s, nMe, nCH4, { color: C.accent6, width: W(ff.Z10_ch4, mx), label: fmt(ff.Z10_ch4), lw: 0.55, lx: 0.4, name: "Z10" });
edge(s, nOle, nPar, { color: CLS.par, width: W(ff.OX3_hyd + 3 * ff.Z7_HT, mx),
  label: `加氢 ${fmt(ff.OX3_hyd)} + 氢转移 ${fmt(3 * ff.Z7_HT)}`, lw: 2.5, lx: 1.45, name: "OX3 + Z7" });
edge(s, nOle, nNap, { color: CLS.nap, width: W(ff.Z6_cyc, mx), label: `环化 ${fmt(ff.Z6_cyc)}`, lw: 1.05, lx: -0.55, ly: 0.12, name: "Z6" });
edge(s, nNap, nAro, { color: CLS.aro, width: W(ff.Z7b_dh + ff.Z7_HT, mx), label: `脱氢 ${fmt(ff.Z7b_dh)} · 氢转移 ${fmt(ff.Z7_HT)}`,
  lw: 2.2, lx: 0.6, ly: 0.22, name: "Z7b + Z7" });
{ // dealkylation aromatics -> olefins (offset parallel line)
  const A = { ...nAro, cy: nAro.cy - 0.32 }, B = { ...nOle, cy: nOle.cy - 0.32 };
  edge(s, A, B, { color: CLS.aro, width: W(ff.Z9_dealk, mx), dash: "dash", label: `脱烷基 ${fmt(ff.Z9_dealk)}`, lw: 1.2, ly: -0.22, name: "Z9" });
}
s.addText("数字 = 每转化 1 mol CO₂ 经过该步的 mol 数（G1 校准：320 °C、3 MPa）；线宽 ∝ 通量", {
  x: 0.55, y: 6.65, w: 12.2, h: 0.3, fontSize: 11, color: C.accent6, margin: 0, isTextBox: true });
s.addNotes("左边灰色区是氧化物，右边是分子筛。88% 的转化碳先在氧化物上变成甲醇，12% 经 RWGS 变成 CO 离开。甲醇（含 DME）在分子筛上主要做烯烃甲基化（0.73），烯烃池内部裂解（0.135）比齐聚（0.007）快得多。烯烃一部分被氧化物加氢或经氢转移变成烷烃，一部分环化后脱氢/氢转移成芳烃，芳烃再被甲醇甲基化到 C8–C10，少量脱烷基回到轻烯烃。");

// 3. families + scale
s = pres.addSlide({ masterName: "CONTENT", sectionTitle: "总览" });
s.addText("11 个反应族自动生成 296 个集总反应", { placeholder: "title" });
const stats = [["72", "集总物种", "烯烃 C2–C20、烷烃 C1–C20、环状 C6–C20、芳烃 C6–C16 及小分子"],
               ["299", "反应", "296 个集总反应 + 甲醇合成、RWGS、DME 3 个可逆步骤"],
               ["0.1 s", "单次反应器计算", "BDF 刚性求解，可做集合拟合与扫描"]];
stats.forEach(([big, lab, sub], i) => {
  const y = 1.5 + i * 1.75;
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 0.55, y, w: 3.6, h: 1.5, rectRadius: 0.1, fill: { color: C.background2 },
    line: { color: C.background2 }, objectName: `stat card ${i + 1}` });
  s.addText(big, { x: 0.75, y: y + 0.12, w: 3.2, h: 0.65, fontSize: 36, bold: true, color: C.text2, margin: 0, isTextBox: true });
  s.addText([{ text: lab, options: { bold: true, fontSize: 14, color: C.text1, breakLine: true } },
             { text: sub, options: { fontSize: 11, color: C.accent6 } }],
    { x: 0.75, y: y + 0.75, w: 3.25, h: 0.7, margin: 0, valign: "top", isTextBox: true });
});
const fam = [
  ["OX1", "CO₂ + 3H₂ ⇌ CH₃OH + H₂O", "1", ox.MeOH_synthesis],
  ["OX2", "CO₂ + H₂ ⇌ CO + H₂O（RWGS）", "1", ox.RWGS],
  ["OX3", "CₙH₂ₙ + H₂ → CₙH₂ₙ₊₂（加氢）", "19", ff.OX3_hyd],
  ["Z1", "2CH₃OH ⇌ DME + H₂O", "1", ox.DME_net],
  ["Z2", "2CH₃OH → C₂H₄（引发）", "1", ff.Z2_init],
  ["Z3", "烯烃 + CH₃OH/DME → 烯烃(n+1)", "36", ff.Z3_meth],
  ["Z4", "烯烃(i) + 烯烃(j) → 烯烃(i+j)", "81", ff.Z4_olig],
  ["Z5", "烯烃(n) → 烯烃(i) + 烯烃(n−i)", "80", ff.Z5_crack],
  ["Z6", "烯烃(n≥6) → 环状(n)", "15", ff.Z6_cyc],
  ["Z7", "环状 + 3 烯烃 → 芳烃 + 3 烷烃", "15", ff.Z7_HT],
  ["Z7b", "环状 → 芳烃 + H₂（Ga/Zn 位）", "11", ff.Z7b_dh],
  ["Z8", "芳烃 + CH₃OH/DME → 芳烃(n+1)", "20", ff.Z8_ameth],
  ["Z9", "芳烃(n) → 芳烃(n−k) + C₂=/C₃=", "17", ff.Z9_dealk],
  ["Z10", "CH₃OH + H₂ → CH₄ + H₂O", "1", ff.Z10_ch4],
];
const hdr = ["族", "反应", "反应数", "G1 通量"].map((t) => ({ text: t, options: { bold: true, color: C.background1, fill: { color: C.text2 } } }));
const rows = [hdr, ...fam.map((r, i) => [r[0], r[1], r[2], fmt(r[3])].map((t, j) => ({ text: t, options: {
  fill: { color: i % 2 ? C.background1 : C.background2 }, align: j >= 2 ? "right" : "left", bold: j === 0 } })))];
s.addTable(rows, { x: 4.55, y: 1.5, w: 8.25, colW: [0.75, 5.0, 1.0, 1.5], rowH: 0.33, fontSize: 12, color: C.text1,
  border: { type: "none" }, margin: [0, 0.08, 0, 0.08], valign: "middle", objectName: "family table" });
s.addText("通量单位：mol / mol CO₂ 转化；Z1 为 DME 净生成", { x: 4.55, y: 6.55, w: 8.25, h: 0.3, fontSize: 11, color: C.accent6, margin: 0, isTextBox: true });
s.addNotes("网络规模：72 个物种、296 个集总反应，加上甲醇合成、RWGS、DME 三个可逆步骤。数目最多的是齐聚（81）和裂解（80），但通量最大的是烯烃甲基化（0.73）。");

// 4. carbon-number x class grid
pres.addSection({ title: "网络细节" });
s = pres.addSlide({ masterName: "CONTENT", sectionTitle: "网络细节" });
s.addText("碳数 × 类别网络：碳沿烯烃甲基化增长，再转入芳烃", { placeholder: "title" });
const X = (n) => 1.95 + (n - 1) * 0.68;
const RY = { par: 1.95, ole: 3.3, nap: 4.6, aro: 5.85 };
const RLAB = { par: "烷烃", ole: "烯烃", nap: "环状", aro: "芳烃" };
for (let n = 1; n <= 16; n++)
  s.addText(`C${n}`, { x: X(n) - 0.3, y: 1.4, w: 0.6, h: 0.28, fontSize: 10, color: C.accent6, align: "center", margin: 0, isTextBox: true });
const maxShare = Math.max(...Object.values(D.share));
const gnodes = {};
for (const [cls, y] of Object.entries(RY)) {
  s.addText(RLAB[cls], { x: 0.55, y: y - 0.2, w: 0.9, h: 0.4, fontSize: 14, bold: true, color: CLS[cls], margin: 0, isTextBox: true });
  const [lo, hi] = { par: [1, 16], ole: [2, 16], nap: [6, 16], aro: [6, 16] }[cls];
  for (let n = lo; n <= hi; n++) {
    const sh = D.share[`${cls}${n}`] || 0;
    const d = sh > 1e-4 ? 0.16 + 0.44 * Math.sqrt(sh / maxShare) : 0.12;
    gnodes[`${cls}${n}`] = { cx: X(n), cy: y, w: d, h: d, shape: "ellipse" };
    s.addShape(pres.shapes.OVAL, { x: X(n) - d / 2, y: y - d / 2, w: d, h: d,
      fill: sh > 1e-4 ? { color: CLS[cls] } : { color: C.background1 }, line: { color: CLS[cls], width: 1 },
      objectName: `${cls}${n} share ${(sh * 100).toFixed(1)}%` });
  }
}
const maxE = Math.max(...Object.values(D.edges));
for (const [k, f] of Object.entries(D.edges)) {
  if (f < 8e-4) continue;
  const [a, b] = k.split(">");
  if (!gnodes[a] || !gnodes[b]) continue;
  const col = b.startsWith("par") ? CLS.par : b.startsWith("nap") ? CLS.nap : b.startsWith("aro") && a.startsWith("aro") ? CLS.aro
            : b.startsWith("aro") ? CLS.aro : CLS.ole;
  edge(s, gnodes[a], gnodes[b], { color: col, width: W(f, maxE, 0.6, 7), gap: 0.02, name: `flux ${k}` });
}
s.addText("每步甲基化 +CH₃OH/DME", { x: X(9) - 0.2, y: RY.ole - 0.62, w: 2.6, h: 0.28, fontSize: 11, color: C.accent6, margin: 0, isTextBox: true });
s.addText("节点面积 ∝ 出口碳选择性（G1 校准）；箭头粗细 ∝ 通量。→ 甲基化  ↑ 加氢/氢转移生成烷烃  ↓ 环化、脱氢芳构化。裂解与齐聚见下页", {
  x: 0.55, y: 6.5, w: 12.2, h: 0.4, fontSize: 11, color: C.accent6, margin: 0, isTextBox: true });
s.addNotes("横轴是碳数，纵向四行是四类物种。节点越大，出口碳选择性越高：最大的是 C8–C10 芳烃。主干通量是烯烃行上从 C3 到 C9 的甲基化箭头；C6–C8 烯烃向下环化，再脱氢成芳烃，芳烃沿底行继续甲基化到 C9–C10。C11 以上几乎没有通量：重烯烃在长到孔径上限前就被裂解（下一页）。");

// 5. cracking / oligomerization arcs (image: network form has no native PowerPoint chart)
s = pres.addSlide({ masterName: "CONTENT", sectionTitle: "网络细节" });
s.addText("烯烃池内部：裂解把重烯烃打回 C3–C6，齐聚几乎不发生", { placeholder: "title" });
s.addImage({ path: __dirname + "/arcs.png", x: 0.55, y: 1.75, w: 8.9, h: 8.9 * ARC_H / ARC_W,
  altText: "烯烃 C2=–C16= 之间的裂解（上方橙色弧）与齐聚（下方蓝色弧）通量弧线图，线宽与通量成正比", objectName: "arc diagram" });
const legendY = 1.75 + 8.9 * ARC_H / ARC_W + 0.4;
[[CLS.ole, "上方弧：裂解 CₙH₂ₙ → CᵢH₂ᵢ + Cₙ₋ᵢH₂₍ₙ₋ᵢ₎（80 个反应）"],
 [CLS.par, "下方弧：齐聚 CᵢH₂ᵢ + CⱼH₂ⱼ → Cᵢ₊ⱼH₂₍ᵢ₊ⱼ₎（81 个反应）"]].forEach(([col, t], i) => {
  line(s, 0.7, legendY + i * 0.42 + 0.14, 1.3, legendY + i * 0.42 + 0.14, { color: col, width: 4, head: false, name: `legend ${i + 1}` });
  s.addText(t, { x: 1.45, y: legendY + i * 0.42, w: 7.8, h: 0.3, fontSize: 12, color: C.text1, margin: 0, isTextBox: true });
});
s.addText("线宽 ∝ 积分通量（G1 校准）", { x: 1.45, y: legendY + 0.84, w: 7, h: 0.3, fontSize: 11, color: C.accent6, margin: 0, isTextBox: true });
s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 9.75, y: 1.75, w: 3.05, h: 4.3, rectRadius: 0.1, fill: { color: C.background2 },
  line: { color: C.background2 }, objectName: "takeaway card" });
s.addText([
  { text: `${fmt(ff.Z5_crack)} vs ${fmt(ff.Z4_olig)}`, options: { fontSize: 26, bold: true, color: C.text2, breakLine: true } },
  { text: "裂解通量约为齐聚的 18 倍", options: { fontSize: 13, color: C.text1, breakLine: true } },
  { text: " ", options: { fontSize: 8, breakLine: true } },
  { text: "C7–C11 烯烃主要裂解成 C3/C4/C5，再重新进入甲基化，形成烯烃循环。", options: { fontSize: 13, color: C.text1, breakLine: true } },
  { text: " ", options: { fontSize: 8, breakLine: true } },
  { text: "C12 以上没有通量：碳数上端由裂解速率决定，而不是孔径。", options: { fontSize: 13, color: C.text1 } },
], { x: 9.95, y: 1.95, w: 2.7, h: 3.9, valign: "top", margin: 0, isTextBox: true });
s.addNotes("这页把 161 个两两反应画成弧线图：上方橙色是裂解，从母体烯烃连到每个碎片；下方蓝色是齐聚，从两个单体连到产物。裂解通量 0.135 mol/mol，齐聚只有 0.007。重烯烃在 C7–C11 就被打回 C3–C5，这就是模型里加大孔径也拿不到 C12–C16 的原因。");

// 6. explicit C1 entry network
s = pres.addSlide({ masterName: "CONTENT", sectionTitle: "网络细节" });
s.addText("C1 入口网络（显式，含自由基）：CO₂ 如何一步步加氢成甲醇", { placeholder: "title" });
const R = (cx, cy, t) => ({ cx, cy, w: 1.18, h: 0.46, fill: C.background1, line: C.accent5, tc: C.accent5, text: t, fs: 13, name: t });
const M = (cx, cy, t) => ({ cx, cy, w: 1.18, h: 0.46, fill: C.accent5, text: t, fs: 13, name: t });
const n1 = {
  CO2: M(1.15, 3.75, "CO_2"), HOCO: R(2.85, 2.9, "HOCO•"), HCOO: R(2.85, 4.6, "HCOO•"),
  CO: M(4.7, 1.8, "CO"), HCOOH: M(4.7, 3.75, "HCOOH"), HCO: R(6.55, 1.8, "HCO•"),
  HOCHOH: R(6.55, 3.2, "HOC•HOH"), H2COHO: R(6.55, 4.6, "H_2C(OH)O•"), CH2OH2: M(8.4, 5.4, "CH_2(OH)_2"),
  CH2O: M(8.4, 2.65, "CH_2O"), CH3O: R(10.1, 1.9, "CH_3O•"), CH2OH: R(10.1, 3.6, "•CH_2OH"),
  CH3OH: M(11.75, 2.75, "CH_3OH"), CH3: R(11.75, 4.4, "•CH_3"), CH4: M(11.75, 5.6, "CH_4"),
};
const ET = { H: { color: C.accent1 }, mol: { color: C.accent6, dash: "dash" }, oh: { color: C.accent2 }, sh: { color: C.accent3, dash: "sysDot", both: true } };
const E1 = [
  ["CO2", "HOCO", "H", "+H"], ["CO2", "HCOO", "H", "+H"], ["HOCO", "HCOO", "sh", ""], ["HOCO", "HCOOH", "H", "+H"],
  ["HCOO", "HCOOH", "H", "+H"], ["HOCO", "CO", "oh", "−OH"], ["CO2", "HCOOH", "mol", "+H_2"], ["HCOOH", "CO", "mol", "−H_2O"],
  ["HCOOH", "HOCHOH", "H", "+H"], ["HCOOH", "H2COHO", "H", "+H"], ["HOCHOH", "H2COHO", "sh", ""], ["HOCHOH", "CH2OH2", "H", "+H"],
  ["H2COHO", "CH2OH2", "H", "+H"], ["H2COHO", "CH2O", "oh", "−OH"], ["CH2OH2", "CH2O", "mol", "−H_2O"], ["HCO", "HCOOH", "oh", "+OH"],
  ["CO", "HCO", "H", "+H"], ["HCO", "CH2O", "H", "+H"], ["CH2O", "CO", "mol", "−H_2"], ["CH2O", "CH3O", "H", "+H"],
  ["CH2O", "CH2OH", "H", "+H"], ["CH3O", "CH2OH", "sh", ""], ["CH3O", "CH3OH", "H", "+H"], ["CH2OH", "CH3OH", "H", "+H"],
  ["CH2O", "CH3OH", "mol", "+H_2"], ["CH3OH", "CH3", "oh", "−OH"], ["CH3", "CH4", "H", "+H"], ["HCOOH", "CH2OH2", "mol", "+H_2"],
  ["CH2OH", "CH2OH2", "oh", "+OH"],
];
const lblPos = { "CO2>HCOOH": { lt: 0.25 }, "CH2O>CO": { lt: 0.3 }, "CH2O>CH3OH": { lt: 0.8 }, "HCOOH>CH2OH2": { lt: 0.3 },
  "HCO>HCOOH": { lt: 0.3 }, "H2COHO>CH2O": { lt: 0.6 }, "HOCHOH>CH2OH2": { lt: 0.65 }, "HCOOH>H2COHO": { lt: 0.55, perp: -0.17 }, "HCOOH>CO": { perp: 0.32 }, "CH3OH>CH3": { perp: 0.32 },
  "CH2OH2>CH2O": { perp: 0.32 }, "CH3>CH4": { perp: 0.3 } };
E1.forEach(([a, b, t, lab], i) => edge(s, n1[a], n1[b], { ...ET[t], width: 1.6, label: lab || undefined, lw: 0.5, lh: 0.24, lfs: 11, perp: 0.17, lfill: null,
  ...(lblPos[`${a}>${b}`] || {}), name: `C1 edge ${i + 1} ${a}-${b}` }));
Object.values(n1).forEach((n) => node(s, n));
const LG = [["H", "加 H（自由基加成/复合）"], ["mol", "分子步骤（四中心/1,1-消除）"], ["oh", "OH 参与（复合/β-断裂/CO+OH）"], ["sh", "分子内 H 迁移"]];
const LGX = [0.6, 3.55, 6.85, 10.45], LGW = [2.3, 2.6, 2.9, 2.0];
LG.forEach(([t, lab], i) => {
  const x = LGX[i];
  line(s, x, 6.3, x + 0.5, 6.3, { ...ET[t], width: 2, head: t !== "sh" ? undefined : false, both: false, name: `C1 legend ${i + 1}` });
  s.addText(lab, { x: x + 0.6, y: 6.15, w: LGW[i], h: 0.3, fontSize: 11, color: C.text1, margin: 0, isTextBox: true });
});
s.addText("显式网络 19 个物种、141 个反应；图中 29 步，未画 111 个双分子 H 转移（抽提 45、歧化 57、取代 9）。实心 = 分子，空心 = 自由基", {
  x: 0.6, y: 6.55, w: 12.2, h: 0.4, fontSize: 11, color: C.accent6, margin: 0, isTextBox: true });
s.addNotes("这是上一阶段显式枚举得到的 C1 子网络，对应氧化物上 CO2 活化的部分（在表面上，自由基对应吸附态）。两条入口：CO2 加 H 得到 HOCO（羧基）或 HCOO（甲酸根）。羧基路线经 HOCO → CO + OH 走向 CO，再加氢到 HCO、CH2O；甲酸根路线经 HCOOH、H2C(OH)O、CH2O 到 CH3O/CH2OH，最后到甲醇。集总模型里这整块被合并成 OX1（甲醇合成）和 OX2（RWGS）两个步骤。");

// 7. flux comparison chart
pres.addSection({ title: "通量" });
s = pres.addSlide({ masterName: "CONTENT", sectionTitle: "通量" });
s.addText("三个催化剂的网络通量：RWGS 份额决定 CO，芳烃甲基化决定芳烃", { placeholder: "title" });
const labs = D.fam_labels.slice().reverse();
const ser = [["G1 GaZrOx/Ga-ZSM-5（郭组）", "G1_GaZrOx_GaZSM5"], ["G2 ZnO/ZrO₂–ZSM-5（郭组）", "G2_ZnZrO_ZSM5"], ["B1 In₂O₃/HZSM-5（对照）", "B1_In2O3_HZSM5"]]
  .map(([name, key]) => ({ name, labels: labs, values: D.compare[key].slice().reverse().map((v) => Math.round(v * 1000) / 1000) }));
s.addChart(pres.charts.BAR, ser, {
  x: 0.55, y: 1.4, w: 8.1, h: 5.45, barDir: "bar", barGrouping: "clustered", barGapWidthPct: 60,
  chartColors: [HEX.dk2, HEX.accent5, "A9B4B6"],
  showLegend: true, legendPos: "b", legendFontSize: 11, legendFontFace: "+mn-lt",
  catAxisLabelFontSize: 11, valAxisLabelFontSize: 10, catAxisLabelFontFace: "+mn-lt", valAxisLabelFontFace: "+mn-lt",
  catAxisLabelColor: HEX.accent6, valAxisLabelColor: HEX.accent6, valAxisMinVal: 0, valAxisMaxVal: 1, valAxisMajorUnit: 0.2,
  valGridLine: { color: "E2E7E6", size: 0.75 }, catGridLine: { style: "none" },
  showValue: false, showTitle: true, title: "通量（mol / mol CO₂ 转化）", titleFontSize: 12, titleColor: HEX.accent6, titleFontFace: "+mn-lt",
  objectName: "flux comparison chart",
});
const bul = [
  ["RWGS：0.12 / 0.23 / 0.45", "In₂O₃ 的 RWGS 份额最高，对应 44.8 % CO；GaZrOx 最低。"],
  ["芳烃甲基化：0.13 / 0.21 / 0.08", "ZnZrO 体系把更多甲醇送进芳烃循环，对应 70 % 芳烃。"],
  ["脱氢芳构化 > 氢转移", "三组拟合都更依赖放 H₂ 的芳构化；这一点摘要数据还不能唯一确定。"],
];
bul.forEach(([h, t], i) => {
  const y = 1.55 + i * 1.7;
  s.addShape(pres.shapes.OVAL, { x: 8.95, y: y + 0.05, w: 0.32, h: 0.32, fill: { color: C.accent5 }, line: { color: C.accent5 }, objectName: `bullet dot ${i + 1}` });
  s.addText([{ text: h, options: { bold: true, fontSize: 14, color: C.text2, breakLine: true } }, { text: t, options: { fontSize: 13, color: C.text1 } }],
    { x: 9.4, y, w: 3.4, h: 1.5, valign: "top", margin: 0, isTextBox: true });
});
s.addNotes("同一张网络，三套校准参数。差别集中在三处：氧化物上 RWGS 与甲醇合成的比例（决定 CO 选择性），甲醇有多少进入芳烃甲基化（决定芳烃），以及芳构化走脱氢还是氢转移。最后一点在集合拟合中并不唯一，需要完整的产物族组成数据来判别。");

// 8. closing
pres.addSection({ title: "结论" });
s = pres.addSlide({ masterName: "CLOSING_DARK", sectionTitle: "结论" });
s.addText("怎么读这张网络", { placeholder: "title" });
const take = [
  ["1", "主干：CO₂ → CH₃OH → 烯烃甲基化", "73 % 的转化碳走烯烃甲基化；转化率由氧化物上的甲醇合成控制。"],
  ["2", "CO 来自氧化物上的 RWGS", "分子筛把甲醇抽走后，CO 从单独氧化物的 68–94 % 降到 12–45 %。"],
  ["3", "碳数上端由裂解决定", "裂解通量约为齐聚的 18 倍，C12 以上几乎没有通量。"],
  ["4", "芳烃 = 环化 + 脱氢 + 甲基化", "仅靠氢转移，芳烃上限约 50 %；70 % 需要放 H₂ 的通道。"],
];
take.forEach(([n, h, t], i) => {
  const col = i % 2, row = Math.floor(i / 2);
  const x = 0.7 + col * 6.1, y = 2.2 + row * 1.85;
  s.addShape(pres.shapes.OVAL, { x, y, w: 0.65, h: 0.65, fill: { color: C.accent4 }, line: { color: C.accent4 }, objectName: `num ${n}` });
  s.addText(n, { x, y, w: 0.65, h: 0.65, fontSize: 20, bold: true, color: C.text2, align: "center", valign: "middle", margin: 0, isTextBox: true });
  s.addText([{ text: h, options: { bold: true, fontSize: 17, color: C.background1, breakLine: true } },
             { text: t, options: { fontSize: 14, color: C.background2 } }],
    { x: x + 0.85, y: y - 0.05, w: 4.95, h: 1.8, valign: "top", margin: 0, isTextBox: true });
});
s.addText("文件：tandem_kinetics/network_lumped.csv（296 个集总反应）· data/network_full_C1.csv（C1 显式网络）· calibration.json（参数）", {
  x: 0.7, y: 6.55, w: 12, h: 0.4, fontSize: 11, color: C.background2, margin: 0, isTextBox: true });
s.addNotes("四点读法。注意：所有通量都来自摘要级数据校准的模型；芳烃、C8–C16 等未测量量在集合拟合中范围很宽，需要完整的产物族组成数据才能收窄。");

(async () => {
  await pres.writeFile({ fileName: OUT });
  await applyTheme(OUT, THEME);
  console.log("wrote", OUT);
})();
