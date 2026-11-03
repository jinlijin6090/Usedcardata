// One-page, isomer-resolved reaction network for thermal one-step CO2 hydrogenation
// C1-C4 fully explicit (surface species, radicals, carbenium ions, every isomer), C6-C10 alkylbenzenes all isomers,
// C5+ aliphatics lumped by class. Background blocks = catalyst form.
const pptxgen = require("pptxgenjs");
const { applyTheme } = require("/root/.claude/skills/synced/62cc770a-8c2c-45f8-8ae0-ca258c8631d4_456331ad-9fc1-4922-9170-b4c06151a9b5/pptx/scripts/apply_theme.js");
const OUT = process.argv[2] || "isomer.pptx";

const THEME = {
  name: "CO2 Isomer Network", headFontFace: "Microsoft YaHei", bodyFontFace: "Microsoft YaHei",
  colors: { dk1: "14232B", lt1: "FFFFFF", dk2: "3B4A54", lt2: "F3F5F6", accent1: "E07B24", accent2: "D63B3B",
            accent3: "2A78D6", accent4: "1E9E5A", accent5: "8B5A2B", accent6: "C2185B", hlink: "2A78D6", folHlink: "8B5A2B" },
};
const pres = new pptxgen();
const SW = 44, SH = 24.75;
pres.defineLayout({ name: "BIG", width: SW, height: SH });
pres.layout = "BIG";
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
pres.title = "CO2 一步热催化加氢：同分异构体分开的完整反应网络";
const C = pres.SchemeColor;
const CAT = {
  R1: { c: C.accent1, n: "① 重整 / 合成气 / RWGS" }, R2: { c: C.accent2, n: "② 费托 / C–C 链增长" }, R3: { c: C.accent3, n: "③ 加氢还原" },
  R4: { c: C.accent4, n: "④ 酸催化 / 分子筛" }, R5: { c: C.accent5, n: "⑤ 烃类升级 / 芳烃烷基化" }, R6: { c: C.accent6, n: "⑥ 自由基链反应", dash: "dash" },
  SH: { c: C.text2, n: "共享步骤（吸附 / 脱附 / 异构）" },
};
const CATAL = {
  OX:   { fill: "DCE7F6", c: "2A5DA8", s: "还原性氧化物", n: "还原性氧化物 · 氧空位（In₂O₃ / ZnZrOₓ / GaZrOₓ）" },
  CU:   { fill: "D5EEEE", c: "1C7272", s: "Cu 基 / 界面", n: "Cu 基 / Cu–氧化物界面" },
  MET:  { fill: "FBE7D2", c: "A85810", s: "过渡金属", n: "过渡金属 Ni / Co / Ru / Fe" },
  CARB: { fill: "F6DCDC", c: "A12F2F", s: "金属碳化物", n: "金属碳化物 χ-Fe₅C₂ / Co₂C（+K / Na）" },
  DUAL: { fill: "F2DDEC", c: "86336F", s: "CO 插入双位", n: "CO 插入双位（Co⁰–Co₂C / CuCo / Cu–Fe）" },
  SUP:  { fill: "E9E6DC", c: "655D47", s: "载体 / 酸碱氧化物", n: "ZrO₂ / CeO₂（酮化）" },
  BA:   { fill: "DDEFE1", c: "1D7342", s: "B 酸分子筛", n: "Brønsted 酸分子筛（H-ZSM-5 / SAPO-34）" },
  LA:   { fill: "E9F0CC", c: "5A7410", s: "Lewis 酸改性分子筛", n: "Ga / Zn 改性分子筛（脱氢）" },
  BIF:  { fill: "F2E6D6", c: "7E5226", s: "双功能 金属–酸", n: "双功能 金属–酸（Pt / Ni + 分子筛）" },
  GAS:  { fill: "FBE4EC", c: "B0175A", s: "无催化剂（气相）", n: "无需催化剂：气相自由基（高温 / 表面引发后脱附）" },
};

pres.defineSlideMaster({ title: "BIG", background: { color: C.background1 },
  objects: [{ placeholder: { options: { name: "title", type: "title", x: 0.4, y: 0.12, w: 30, h: 0.7, fontSize: 30, bold: true,
    color: C.text1, align: "left", valign: "middle", margin: 0 }, text: "" } }] });

// ------------------------------------------------------------------ helpers (inches, scale 1)
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
function textW(s, pt) {
  let em = 0;
  const plain = s.replace(/_\{([^}]*)\}|_(\d+|[a-z])|\^\{([^}]*)\}/g, (m, a, b, c) => "\u0001".repeat((a ?? b ?? c).length));
  for (const ch of plain) {
    const c = ch.codePointAt(0);
    if (ch === "\u0001") em += 0.46; else if (c > 0x2e80) em += 1.02; else if (/[A-Z0-9]/.test(ch)) em += 0.7;
    else if (/[a-z]/.test(ch)) em += 0.58; else if (ch === " ") em += 0.3; else em += 0.55;
  }
  return em * pt / 72;
}
function seg(slide, x1, y1, x2, y2, o) {
  slide.addShape(pres.shapes.LINE, {
    x: Math.min(x1, x2), y: Math.min(y1, y2), w: Math.max(Math.abs(x2 - x1), 0.001), h: Math.max(Math.abs(y2 - y1), 0.001),
    flipH: x2 < x1, flipV: y2 < y1,
    line: { color: o.color, width: o.width || 1.25, dashType: o.dash || "solid",
            endArrowType: o.end ? "triangle" : undefined, beginArrowType: o.begin ? "triangle" : undefined },
    objectName: o.name,
  });
}
function txt(slide, x, y, w, h, text, o = {}) {
  slide.addText(Array.isArray(text) ? text : runs(text, { fontSize: o.fs || 8, color: o.color || C.text2, bold: !!o.bold, italic: !!o.italic }), {
    x, y, w, h, align: o.align || "left", valign: o.valign || "middle", margin: 0, isTextBox: true, objectName: o.name || "text",
    ...(o.fill ? { fill: { color: o.fill } } : {}) });
}
const NFS = 9, LFS = 7, NH = 0.34;
const N = {}, E = [];
function n(id, type, x, y, label, o = {}) {
  if (N[id]) throw new Error("dup node " + id);
  N[id] = { id, type, x, y, label, ...o };
}
const e = (a, b, cat, label = "", o = {}) => E.push({ a, b, cat, label, ...o });
function size(nd) {
  const fs = nd.fs || NFS;
  return [nd.w || Math.max(0.5, textW(nd.label, fs) + (nd.type === "rad" ? 0.36 : 0.24)), nd.h || NH];
}
function clip(nd, ux, uy) {
  const [w, h] = size(nd);
  if (nd.type === "rad") { const a = w / 2, b = h / 2, t = 1 / Math.sqrt((ux * ux) / (a * a) + (uy * uy) / (b * b)); return [nd.x + ux * t, nd.y + uy * t]; }
  const tx = ux ? w / 2 / Math.abs(ux) : Infinity, ty = uy ? h / 2 / Math.abs(uy) : Infinity, t = Math.min(tx, ty);
  return [nd.x + ux * t, nd.y + uy * t];
}
function label(slide, x, y, t, color, fs = LFS, name = "edge label") {
  const w = textW(t, fs) + 0.1, h = fs / 72 * 1.5;
  txt(slide, x - w / 2, y - h / 2, w, h, runs(t, { fontSize: fs, color }), { align: "center", fill: C.background1, name });
}
function drawEdge(slide, ed, i) {
  const A = N[ed.a], B = N[ed.b], cat = CAT[ed.cat];
  if (!A || !B) throw new Error(`edge ${ed.a}->${ed.b} missing`);
  const st = { color: cat.c, width: ed.width || 1.25, dash: ed.dash || cat.dash };
  const nm = `edge ${i + 1} ${ed.a}→${ed.b}`;
  if (ed.pts) {
    const P = ed.pts;
    for (let k = 0; k < P.length - 1; k++)
      seg(slide, P[k][0], P[k][1], P[k + 1][0], P[k + 1][1], { ...st, end: k === P.length - 2, begin: ed.both && k === 0, name: `${nm} seg ${k + 1}` });
    if (ed.label && ed.lpos) label(slide, ed.lpos[0], ed.lpos[1], ed.label, cat.c, LFS, nm + " label");
    return;
  }
  const ax = A.x + (ed.ax || 0), ay = A.y + (ed.ay || 0), bx = B.x + (ed.bx || 0), by = B.y + (ed.by || 0);
  const dx = bx - ax, dy = by - ay, L = Math.hypot(dx, dy), ux = dx / L, uy = dy / L;
  let [x1, y1] = clip({ ...A, x: ax, y: ay }, ux, uy);
  let [x2, y2] = clip({ ...B, x: bx, y: by }, -ux, -uy);
  const g = 0.02; x1 += ux * g; y1 += uy * g; x2 -= ux * g; y2 -= uy * g;
  seg(slide, x1, y1, x2, y2, { ...st, end: true, begin: !!ed.both, name: nm });
  if (ed.label) {
    let nx = -uy, ny = ux;
    if (ny > 0 || (Math.abs(ny) < 1e-9 && nx > 0)) { nx = -nx; ny = -ny; }
    const t = ed.lt ?? 0.5, pp = ed.perp ?? 0;
    label(slide, x1 + (x2 - x1) * t + nx * pp + (ed.lx || 0), y1 + (y2 - y1) * t + ny * pp + (ed.ly || 0), ed.label, cat.c, LFS, nm + " label");
  }
}
function drawNode(slide, nd) {
  const [w, h] = size(nd);
  const x = nd.x - w / 2, y = nd.y - h / 2, fs = nd.fs || NFS;
  if (nd.type === "txt") {
    slide.addText(runs(nd.label, { fontSize: fs, color: nd.tc || C.text2, bold: !!nd.bold, italic: nd.italic !== false }),
      { x, y, w, h, align: nd.align || "center", valign: "middle", margin: 0, isTextBox: true, objectName: `node ${nd.id}` });
    return;
  }
  const S = {
    mol: { shape: pres.shapes.ROUNDED_RECTANGLE, fill: nd.fill || C.text2, line: nd.fill || C.text2, tc: C.background1, r: 0.05 },
    ref: { shape: pres.shapes.ROUNDED_RECTANGLE, fill: "D9DEE1", line: "8A979E", tc: C.text1, r: 0.05, dash: "dash" },
    surf: { shape: pres.shapes.ROUNDED_RECTANGLE, fill: C.background1, line: C.text1, tc: C.text1, r: "pill" },
    rad: { shape: pres.shapes.OVAL, fill: C.background1, line: C.accent6, tc: C.text1, dash: "dash" },
    ion: { shape: pres.shapes.RECTANGLE, fill: C.background1, line: C.text1, tc: C.text1 },
    hub: { shape: pres.shapes.ROUNDED_RECTANGLE, fill: C.accent1, line: C.accent1, tc: C.background1, r: 0.08 },
  }[nd.type];
  if (nd.type === "ion") slide.addShape(pres.shapes.RECTANGLE, { x: x - 0.04, y: y - 0.04, w: w + 0.08, h: h + 0.08,
    fill: { color: C.background1 }, line: { color: C.text1, width: 0.75 }, objectName: `node ${nd.id} outer` });
  const opt = { x, y, w, h, shape: S.shape, fill: { color: S.fill }, line: { color: S.line, width: nd.type === "rad" ? 1.1 : 0.9, dashType: S.dash || "solid" },
    align: "center", valign: "middle", margin: 0, objectName: `node ${nd.id}` };
  if (S.r === "pill") opt.rectRadius = h / 2; else if (S.r !== undefined) opt.rectRadius = S.r;
  slide.addText(runs(nd.label, { fontSize: fs, color: S.tc, bold: true }), opt);
}
const REG = [], BLK = [], NOTES = [];
const region = (x, y, w, h, color, title, dash) => REG.push({ x, y, w, h, color, title, dash });
const block = (k, x, y, w, h, chip, lab) => BLK.push({ k, x, y, w, h, chip, lab });
const note = (x, y, w, h, t, o = {}) => NOTES.push({ x, y, w, h, t, o });

// ================================================================== REGION A: C1 surface core (x 0.4-14.0)
region(0.35, 1.0, 13.75, 11.25, C.accent1, "A · C_1 表面中间体：CO_2 活化 · 甲醇 · RWGS · CO 活化 · 甲烷化");
block("OX", 2.05, 1.5, 11.95, 4.35, [2.1, 1.52], "还原性氧化物 · 氧空位 / Cu–ZnO（甲醇 · RWGS）");
block("MET", 2.05, 6.3, 11.95, 5.85, [10.55, 11.82], "过渡金属 Ni / Co / Ru / Fe（CO 解离 · 甲烷化）");
n("fCO2", "mol", 1.2, 3.1, "CO_2(g)"); n("fH2", "mol", 1.2, 5.0, "H_2(g)"); n("Hs", "surf", 1.2, 6.7, "H*");
n("HH", "surf", 3.4, 2.1, "M–H^{−} / O–H^{+}"); n("VO", "surf", 6.4, 2.1, "V_{O} 氧空位"); n("OL", "surf", 8.6, 2.1, "O_{lat} 晶格氧");
n("CO2s", "surf", 3.2, 3.15, "CO_2^{δ−}*"); n("HCOO", "surf", 4.85, 3.15, "HCOO*"); n("H2COO", "surf", 6.5, 3.15, "H_2COO*");
n("H2COOH", "surf", 8.2, 3.15, "H_2COOH*"); n("H2CO", "surf", 9.9, 3.15, "H_2CO*"); n("CH3O", "surf", 11.4, 3.15, "CH_3O*");
n("CH3OHs", "surf", 12.9, 3.15, "CH_3OH*");
n("HCOOHs", "surf", 4.85, 4.35, "HCOOH*"); n("HCOOHg", "mol", 6.65, 4.35, "HCOOH(g)"); n("MF", "mol", 10.15, 4.6, "HCOOCH_3 甲酸甲酯");
n("MeOHg", "mol", 13.25, 4.35, "CH_3OH(g)");
n("COOH", "surf", 3.2, 5.35, "COOH*"); n("hubCO", "hub", 5.4, 6.05, "CO(g)");
n("CO", "surf", 5.4, 7.0, "CO*"); n("HCO", "surf", 7.4, 7.0, "HCO*"); n("CH2OH", "surf", 11.6, 7.0, "CH_2OH*");
n("COH", "surf", 2.9, 8.3, "COH*"); n("HCOH", "surf", 8.9, 8.3, "HCOH*");
n("Cs", "surf", 3.4, 9.8, "C*"); n("CHs", "surf", 5.4, 9.8, "CH*"); n("CH2s", "surf", 7.4, 9.8, "CH_2*"); n("CH3s", "surf", 9.4, 9.8, "CH_3*");
n("CH4g", "mol", 12.4, 9.8, "CH_4(g)");
n("coke", "mol", 2.7, 11.3, "积碳 C"); n("FeC", "surf", 4.45, 11.3, "χ-Fe_5C_2"); n("Os", "surf", 7.0, 11.3, "O*"); n("OHs", "surf", 8.8, 11.3, "OH*");
n("H2Og", "mol", 12.4, 11.3, "H_2O(g)");
e("fCO2", "CO2s", "SH", "吸附"); e("fH2", "HH", "R3", "异裂", { lt: 0.82, perp: 0.18 }); e("fH2", "Hs", "SH", "解离", { perp: 0.25 });
e("OL", "VO", "R3", "+2H → H_2O(g)", { ly: -0.2 }); e("VO", "CO2s", "SH", "吸附于 V_{O}", { perp: 0.2, lt: 0.4 });
e("CH3OHs", "OL", "SH", "O 回填晶格", { dash: "dash", perp: 0.15, lt: 0.7 });
e("CO2s", "HCOO", "R3", "+H", { both: true, ly: -0.18 }); e("HCOO", "H2COO", "R3", "+H", { ly: -0.18 }); e("H2COO", "H2COOH", "R3", "+H", { ly: -0.18 });
e("H2COOH", "H2CO", "R3", "−OH*", { ly: -0.18 }); e("H2CO", "CH3O", "R3", "+H", { ly: -0.18 }); e("CH3O", "CH3OHs", "R3", "+H", { ly: -0.18 });
e("CH3OHs", "MeOHg", "SH", ""); e("HCOO", "HCOOHs", "R3", "+H", { perp: 0.18 }); e("HCOOHs", "HCOOHg", "SH", "", { both: true });
e("HCOOHs", "H2COOH", "R3", "+H", { perp: -0.15, lt: 0.6 }); e("CH3O", "MF", "R3", "+HCOO*", { perp: 0.2 });
e("CO2s", "COOH", "R3", "+H", { perp: 0.18 }); e("COOH", "CO", "R1", "−OH*（RWGS）", { perp: 0.22, lt: 0.35 });
e("hubCO", "CO", "SH", "", { both: true }); e("CO", "HCO", "R3", "+H", { ly: -0.18 }); e("HCO", "H2CO", "R3", "+H", { perp: 0.18, lt: 0.6 });
e("CO", "COH", "R3", "+H", { perp: 0.18 }); e("HCO", "HCOH", "R3", "+H", { perp: -0.18 }); e("COH", "Cs", "SH", "−OH*", { perp: 0.18 });
e("HCOH", "CHs", "SH", "−OH*", { perp: -0.18, lt: 0.4 }); e("CO", "Cs", "SH", "−O*", { both: true, perp: -0.2 }); e("HCO", "CHs", "SH", "−O*", { both: true, perp: 0.2, lt: 0.6 });
e("H2CO", "CH2s", "SH", "−O*", { both: true, perp: 0.2, lt: 0.75 }); e("CH3O", "CH3s", "SH", "−O*", { both: true, perp: -0.2, lt: 0.8 });
e("HCOH", "CH2OH", "R3", "+H", { perp: 0.18 }); e("CH2OH", "CH3OHs", "R3", "+H", { perp: 0.2, lt: 0.35 });
e("Cs", "CHs", "R3", "+H", { both: true, ly: -0.18 }); e("CHs", "CH2s", "R3", "+H", { both: true, ly: -0.18 }); e("CH2s", "CH3s", "R3", "+H", { both: true, ly: -0.18 });
e("CH3s", "CH4g", "R3", "+H ⇌ C–H 活化", { both: true, ly: -0.2 });
e("Cs", "coke", "R1", "析碳", { perp: 0.18 }); e("Cs", "FeC", "R2", "渗碳", { both: true, perp: -0.2 }); e("Os", "OHs", "SH", "+H", { ly: -0.18 }); e("OHs", "H2Og", "SH", "+H", { ly: -0.18 });

// ================================================================== REGION R: radicals (x 0.4-14.0)
region(0.35, 12.45, 13.75, 11.2, C.accent6, "R · 气相自由基（R• 与 A 中 R* 一一对应）：C_1 加氢链 · C_1–C_4 烷基 / 烯丙基 · H / O 池", "dash");
block("GAS", 0.45, 12.85, 13.55, 10.7, [9.2, 12.9], "无需催化剂：气相自由基（高温 / 表面引发后脱附）");
n("rCO2", "ref", 1.2, 14.4, "CO_2"); n("HOCO", "rad", 3.0, 13.6, "HOCO•"); n("HCOOr", "rad", 3.0, 15.2, "HCOO•");
n("rCO", "ref", 5.0, 13.6, "CO"); n("rHCOOH", "ref", 5.0, 15.2, "HCOOH"); n("HCOr", "rad", 6.8, 13.6, "HCO•"); n("HOCHOH", "rad", 6.8, 15.2, "HOC•HOH");
n("rCH2O", "ref", 8.7, 14.4, "CH_2O"); n("CH3Or", "rad", 10.6, 13.75, "CH_3O•"); n("CH2OHr", "rad", 10.6, 15.05, "•CH_2OH"); n("rMeOH", "ref", 12.7, 14.4, "CH_3OH");
e("rCO2", "HOCO", "R6", "+H•", { perp: 0.15 }); e("rCO2", "HCOOr", "R6", "+H•", { perp: -0.15 }); e("HOCO", "HCOOr", "R6", "H 迁移", { both: true, lx: 0.4 });
e("HOCO", "rCO", "R6", "−OH• / +OH•", { both: true, ly: -0.2 }); e("HCOOr", "rHCOOH", "R6", "+H•", { ly: -0.18 }); e("rCO", "HCOr", "R6", "+H•", { both: true, ly: -0.18 });
e("rHCOOH", "HOCHOH", "R6", "+H•", { ly: -0.18 }); e("HCOr", "rCH2O", "R6", "+H• / 抽氢", { both: true, perp: 0.18 }); e("HOCHOH", "rCH2O", "R6", "−OH•", { perp: -0.18 });
e("rCH2O", "CH3Or", "R6", "+H•", { perp: 0.15 }); e("rCH2O", "CH2OHr", "R6", "+H•", { perp: -0.15 }); e("CH3Or", "CH2OHr", "R6", "", { both: true });
e("CH3Or", "rMeOH", "R6", "+H•", { perp: 0.15 }); e("CH2OHr", "rMeOH", "R6", "+H• / 抽氢", { both: true, perp: -0.18 });
// alkyl radicals by carbon number
n("rCH4", "ref", 1.6, 17.6, "CH_4"); n("CH3r", "rad", 1.6, 19.4, "•CH_3");
n("rC2H6", "ref", 4.3, 17.6, "C_2H_6"); n("C2H5r", "rad", 4.3, 19.4, "C_2H_5•"); n("rC2H4", "ref", 4.3, 21.2, "C_2H_4"); n("C2H3r", "rad", 5.9, 21.2, "C_2H_3•");
n("rC3H8", "ref", 7.6, 17.6, "C_3H_8"); n("nC3r", "rad", 6.8, 19.4, "n-C_3H_7•"); n("iC3r", "rad", 8.5, 19.4, "i-C_3H_7•"); n("rC3H6", "ref", 7.6, 21.2, "C_3H_6");
n("allyl", "rad", 8.6, 22.3, "C_3H_5• 烯丙基");
n("rnC4", "ref", 10.6, 17.6, "n-C_4H_{10}"); n("riC4", "ref", 12.6, 17.6, "i-C_4H_{10}");
n("nC4r", "rad", 10.6, 19.0, "n-C_4H_9•"); n("sC4r", "rad", 12.6, 19.0, "sec-C_4H_9•"); n("iC4r", "rad", 10.6, 20.1, "i-C_4H_9•"); n("tC4r", "rad", 12.6, 20.1, "t-C_4H_9•");
n("r1C4", "ref", 10.3, 21.5, "1-C_4H_8"); n("r2C4", "ref", 11.7, 21.5, "2-C_4H_8"); n("riC4e", "ref", 13.15, 21.5, "i-C_4H_8");
n("rH2", "ref", 1.3, 23.0, "H_2"); n("Hr", "rad", 2.9, 23.0, "H•"); n("OHr", "rad", 4.5, 23.0, "OH•"); n("rH2O", "ref", 6.1, 23.0, "H_2O");
e("rCH4", "CH3r", "R6", "抽氢 ⇌ +H•", { both: true, perp: 0.25 }); e("CH3r", "rC2H6", "R6", "2•CH_3 偶联", { perp: 0.2, lt: 0.45 });
e("rC2H6", "C2H5r", "R6", "抽氢", { both: true, perp: 0.18 }); e("C2H5r", "rC2H4", "R6", "−H• / +H•", { both: true, perp: 0.25 }); e("rC2H4", "C2H3r", "R6", "抽氢", { ly: -0.18 });
e("C2H5r", "rC3H8", "R6", "+•CH_3", { perp: 0.15, lt: 0.35 }); e("rC3H8", "nC3r", "R6", "抽氢（伯）", { both: true, perp: 0.2 }); e("rC3H8", "iC3r", "R6", "抽氢（仲）", { both: true, perp: -0.2 });
e("nC3r", "iC3r", "R6", "H 迁移", { both: true, ly: 0.2 }); e("nC3r", "rC3H6", "R6", "", { both: true }); e("iC3r", "rC3H6", "R6", "−H•", { both: true, perp: -0.15 });
e("nC3r", "rC2H4", "R6", "β-断裂 −•CH_3", { perp: 0.2, lt: 0.55 }); e("rC3H6", "allyl", "R6", "抽氢（烯丙位）", { perp: 0.2 });
e("C2H5r", "rnC4", "R6", "+C_2H_5•", { perp: 0.15, lt: 0.7 }); e("nC3r", "rnC4", "R6", "+•CH_3", { perp: -0.15, lt: 0.75 }); e("iC3r", "riC4", "R6", "+•CH_3", { perp: -0.15, lt: 0.7 });
e("rnC4", "nC4r", "R6", "", { both: true }); e("rnC4", "sC4r", "R6", "抽氢", { both: true, perp: 0.15, lt: 0.3 }); e("riC4", "iC4r", "R6", "抽氢", { both: true, perp: -0.15, lt: 0.3 });
e("riC4", "tC4r", "R6", "", { both: true }); e("nC4r", "sC4r", "R6", "H 迁移", { both: true, ly: -0.18 }); e("iC4r", "tC4r", "R6", "H 迁移", { both: true, ly: 0.2 });
e("nC4r", "r1C4", "R6", "", { both: true, ax: -0.35, bx: -0.15 }); e("sC4r", "r1C4", "R6", "", { both: true, ax: -0.2 }); e("sC4r", "r2C4", "R6", "", { both: true, ax: -0.1 });
e("iC4r", "riC4e", "R6", "", { both: true }); e("tC4r", "riC4e", "R6", "−H•", { both: true, perp: -0.15 });
e("iC4r", "rC3H6", "R6", "β-断裂 −•CH_3", { pts: [[10.0, 20.1], [9.45, 20.1], [9.45, 21.2], [8.05, 21.2]], lpos: [9.45, 20.65] });
e("rH2", "Hr", "R6", "均裂 / 复合", { both: true, ly: -0.2 }); e("OHr", "Hr", "R6", "+H_2", { ly: -0.18 }); e("OHr", "rH2O", "R6", "+H_2 / +H•", { ly: -0.18 });
note(10.0, 22.55, 4.0, 0.9, ["sec-C_4H_9• → C_3H_6 + •CH_3（β-断裂）", "n-C_4H_9• → C_2H_4 + C_2H_5•", "抽氢试剂：H•、OH•、•CH_3"], { color: C.accent6 });

// ================================================================== REGION F: chain growth, C2-C4 isomers, oxygenates, upgrading (x 14.3-29.0)
region(14.3, 1.0, 14.75, 22.65, C.accent2, "F · 金属位链增长（C_2–C_4 异构体逐个显式）· CO 插入含氧化合物 · 双功能升级");
block("CARB", 14.4, 1.4, 14.55, 7.5, [14.45, 8.5], "金属碳化物 χ-Fe₅C₂ / Co₂C（+K / Na）：CHₓ 插入 · β-H 消除 · 终止");
block("DUAL", 14.4, 9.05, 11.3, 5.6, [14.45, 14.3], "CO 插入双位（Co⁰–Co₂C / CuCo / Cu–Fe）：酰基 → 醛 → 醇");
block("SUP", 17.45, 9.6, 1.6, 3.85, null);
block("BIF", 14.4, 15.4, 14.55, 8.15, [14.45, 15.45], "双功能 金属–酸（Pt / Ni + 分子筛）：加氢异构 / 裂化 · 烯烃齐聚 → 航煤");
// paraffins / alkyls / pi-olefins / olefins / dienes
n("C2H6", "mol", 16.4, 2.0, "C_2H_6"); n("C3H8", "mol", 20.0, 2.0, "C_3H_8"); n("nC4H10", "mol", 23.0, 2.0, "n-C_4H_{10}"); n("iC4H10", "mol", 26.6, 2.0, "i-C_4H_{10}");
n("C5P", "mol", 28.35, 2.0, "C_{5+} 正 / 异构烷", { fs: 8 });
n("C2H5s", "surf", 16.4, 3.6, "C_2H_5*"); n("nC3s", "surf", 19.4, 3.6, "n-C_3H_7*"); n("nC4s", "surf", 23.0, 3.6, "n-C_4H_9*"); n("C5s", "surf", 28.25, 3.6, "C_{5+}H_{2n+1}*", { fs: 8 });
n("vinyl", "surf", 18.1, 4.4, "CH=CH_2* 乙烯基", { fs: 8 }); n("iC3s", "surf", 20.9, 4.4, "i-C_3H_7*"); n("sC4s", "surf", 24.4, 4.4, "sec-C_4H_9*");
n("iC4s", "surf", 25.95, 4.4, "i-C_4H_9*"); n("tC4s", "surf", 27.45, 4.4, "t-C_4H_9*");
n("C2H4s", "surf", 16.4, 5.2, "C_2H_4*"); n("C3H6s", "surf", 20.0, 5.2, "C_3H_6*"); n("allylS", "surf", 21.75, 5.2, "C_3H_5* 烯丙基", { fs: 8 });
n("C4H8s", "surf", 23.7, 5.2, "1- / 2-C_4H_8*"); n("iC4H8s", "surf", 26.6, 5.2, "i-C_4H_8*"); n("C5O", "mol", 28.35, 5.9, "C_{5+} α- / 内烯", { fs: 8 });
n("CCH3", "surf", 18.1, 6.05, "CCH_3* 次乙基", { fs: 8 });
n("C2H4", "mol", 16.4, 6.8, "C_2H_4"); n("C3H6", "mol", 20.0, 6.8, "C_3H_6"); n("b1", "mol", 22.85, 6.8, "1-C_4H_8"); n("bc", "mol", 24.25, 6.8, "cis-2-C_4H_8");
n("bt", "mol", 25.75, 6.8, "trans-2-C_4H_8"); n("ib", "mol", 27.25, 6.8, "i-C_4H_8");
n("C2H2", "mol", 16.4, 8.15, "C_2H_2"); n("BD", "mol", 24.25, 8.15, "1,3-C_4H_6 丁二烯");
e("C2H5s", "C2H6", "R2", "+H*", { perp: 0.18 }); e("C2H5s", "nC3s", "R2", "+CH_2*", { ly: -0.18, lt: 0.4 }); e("nC3s", "nC4s", "R2", "+CH_2*", { ly: -0.18 });
e("nC4s", "C5s", "R2", "+CH_2* …", { ly: -0.18, lt: 0.3 }); e("C5s", "C5P", "R2", "", {}); e("C5s", "C5O", "R2", "", { ax: 0.2, bx: 0.2 });
e("C2H5s", "C2H4s", "R2", "β-H", { both: true, perp: 0.18 }); e("C2H4s", "C2H4", "SH", "", { both: true }); e("C2H4s", "CCH3", "SH", "", { both: true });
e("C2H4s", "C2H2", "R2", "−2H", { perp: 0.2, ax: -0.25, bx: -0.25 });
e("vinyl", "C2H4s", "R3", "+H", { perp: 0.15 }); e("vinyl", "allylS", "R2", "+CH_2*（烯基机理）", { perp: 0.17, lt: 0.3 });
e("nC3s", "C3H8", "R2", "+H*", { perp: 0.18 }); e("iC3s", "C3H8", "R2", "+H*", { perp: -0.15, lt: 0.62 }); e("nC3s", "C3H6s", "R2", "", { both: true });
e("iC3s", "C3H6s", "R2", "", { both: true }); e("allylS", "C3H6s", "R3", "+H", { ly: -0.18 }); e("C3H6s", "C3H6", "SH", "", { both: true });
e("iC3s", "iC4s", "R2", "+CH_2*（支化）", { perp: -0.17, lt: 0.33 });
e("nC4s", "nC4H10", "R2", "+H*", { perp: 0.18 }); e("sC4s", "nC4H10", "R2", "", {}); e("iC4s", "iC4H10", "R2", "", {}); e("tC4s", "iC4H10", "R2", "+H*", { perp: -0.15 });
e("nC4s", "C4H8s", "R2", "", { both: true }); e("sC4s", "C4H8s", "R2", "", { both: true }); e("iC4s", "iC4H8s", "R2", "", { both: true }); e("tC4s", "iC4H8s", "R2", "", { both: true });
e("nC4s", "sC4s", "SH", "", { both: true, ax: 0.3 }); e("iC4s", "tC4s", "SH", "", { both: true });
e("C4H8s", "b1", "SH", "", { both: true }); e("C4H8s", "bc", "SH", "", { both: true }); e("C4H8s", "bt", "SH", "", { both: true }); e("iC4H8s", "ib", "SH", "", { both: true });
e("b1", "bc", "SH", "", { both: true, ay: 0.12, by: 0.12 }); e("bc", "bt", "SH", "顺反", { both: true, ly: 0.25 });
e("b1", "BD", "R3", "−H_2", { perp: 0.15 });
// oxygenates (CO insertion)
const OXC = [["C2", 16.4, "CH_3CO*", "CH_3CHO", "C_2H_5O*", "C_2H_5OH"], ["C3", 20.0, "C_2H_5CO*", "C_2H_5CHO", "n-C_3H_7O*", "1-C_3H_7OH"],
             ["iC4", 22.55, "i-C_3H_7CO*", "i-C_3H_7CHO", "i-C_4H_9O*", "i-C_4H_9OH"], ["nC4", 24.75, "n-C_3H_7CO*", "n-C_3H_7CHO", "n-C_4H_9O*", "1-C_4H_9OH"]];
OXC.forEach(([k, x, ac, al, ao, oh]) => {
  n(`ac${k}`, "surf", x, 10.0, ac); n(`al${k}`, "mol", x, 11.4, al); n(`ao${k}`, "surf", x, 12.8, ao); n(`oh${k}`, "mol", x, 14.2, oh);
  e(`ac${k}`, `al${k}`, "R3", k === "C2" ? "+H" : "", { perp: 0.15 }); e(`al${k}`, `ao${k}`, "R3", k === "C2" ? "+H" : "", { perp: 0.15 });
  e(`ao${k}`, `oh${k}`, "R3", k === "C2" ? "+H" : "", { perp: 0.15 });
});
n("AcOH", "mol", 18.25, 10.0, "CH_3COOH"); n("acetone", "mol", 18.25, 11.4, "丙酮"); n("iPrOH", "mol", 18.25, 12.8, "2-C_3H_7OH"); n("MeOAc", "mol", 18.25, 14.2, "CH_3COOCH_3");
n("C5OH", "mol", 27.3, 12.1, "C_{5+} 醇（同理）", { fs: 8 });
e("acC2", "AcOH", "R3", "+OH*", { ly: -0.18 }); e("AcOH", "acetone", "R5", "酮化", { perp: 0.15 }); e("acetone", "iPrOH", "R3", "+H_2", { perp: 0.15 });
e("acC2", "MeOAc", "R3", "+CH_3O*", { perp: 0.17, lt: 0.75 });
e("C2H5s", "acC3", "R2", "+CO*", { pts: [[16.75, 3.77], [16.95, 3.97], [16.95, 9.55], [19.8, 9.55], [19.8, 9.83]], lpos: [18.6, 9.55] });
e("nC3s", "acnC4", "R2", "+CO*", { pts: [[19.4, 3.77], [19.4, 9.4], [24.75, 9.4], [24.75, 9.83]], lpos: [23.8, 9.4] });
e("iC3s", "aciC4", "R2", "+CO*", { pts: [[20.85, 4.57], [20.85, 9.25], [22.55, 9.25], [22.55, 9.83]], lpos: [21.7, 9.25] });
e("C5s", "C5OH", "R2", "+CO* → 醛 → 醇", { pts: [[28.85, 3.6], [28.95, 3.6], [28.95, 12.1], [28.35, 12.1]], lpos: [28.3, 10.6], dash: "dash" });
// upgrading (bifunctional) + olefin oligomerization
n("wax", "mol", 27.6, 16.6, "C_{20+} 蜡"); n("nC12", "mol", 25.3, 16.6, "n-C_{12}H_{26}（代表）"); n("nC12e", "mol", 22.9, 16.6, "n-C_{12}H_{24}"); n("sec12", "ion", 20.7, 16.6, "sec-C_{12}H_{25}^{+}");
n("mono12", "ion", 20.7, 18.1, "单支链 C_{12}^{+}"); n("di12", "ion", 20.7, 19.6, "双支链 C_{12}^{+}"); n("tri12", "ion", 20.7, 21.1, "三支链 C_{12}^{+}");
n("im", "mol", 24.4, 18.1, "单支链 i-C_{12}"); n("id", "mol", 24.4, 19.6, "双支链 i-C_{12}"); n("it", "mol", 24.4, 21.1, "三支链 i-C_{12}");
n("frag", "ion", 17.3, 19.6, "C_5–C_7^{+} + C_5–C_7^{=}"); n("lightP", "mol", 15.6, 21.1, "C_5–C_8 异构烷");
n("jet", "mol", 24.4, 22.85, "C_8–C_{16} 异构烷烃航煤 SPK", { fill: C.accent5 }); n("fOle", "ref", 15.6, 22.85, "C_3^{=} / C_4^{=}"); n("olig", "mol", 19.4, 22.85, "C_8–C_{16}^{=} 低聚物");
e("C5s", "wax", "R2", "链增长 → 蜡", { pts: [[28.85, 3.75], [28.85, 3.75], [28.85, 16.6], [28.17, 16.6]], lpos: [28.3, 15.9] });
e("wax", "nC12", "R5", "加氢裂化", { ly: -0.2 }); e("nC12", "nC12e", "R5", "−H_2 / +H_2", { both: true, ly: -0.2 }); e("nC12e", "sec12", "R5", "+H^{+}", { both: true, ly: -0.2 });
e("sec12", "mono12", "R5", "PCP", { perp: 0.2 }); e("mono12", "di12", "R5", "PCP", { perp: 0.2 }); e("di12", "tri12", "R5", "PCP", { perp: 0.2 });
e("mono12", "im", "R5", "−H^{+} +H_2", { ly: -0.2 }); e("di12", "id", "R5", "", {}); e("tri12", "it", "R5", "", {});
e("mono12", "frag", "R5", "C 型", { perp: -0.18, lt: 0.6 }); e("di12", "frag", "R5", "B 型", { ly: -0.2 }); e("tri12", "frag", "R5", "A 型（最快）", { perp: 0.2, lt: 0.6 });
e("frag", "lightP", "R5", "−H^{+} +H_2", { perp: 0.25 }); e("it", "jet", "R5", "→ 航煤", { perp: 0.2 }); e("im", "id", "SH", "", { both: true, ax: 0.65, bx: 0.65 }); e("id", "it", "SH", "", { both: true, ax: 0.65, bx: 0.65 });
e("fOle", "olig", "R4", "酸催化齐聚", { ly: -0.2 }); e("olig", "jet", "R3", "+H_2", { ly: -0.2 });

// ================================================================== REGION Z: zeolite, methanol -> olefins, explicit carbenium ions (x 29.3-43.6)
region(29.3, 1.0, 14.35, 11.25, C.accent4, "Z · 分子筛：甲醇活化 · 第一个 C–C 键 · 烯烃循环（C_2–C_4 碳正离子与异构体显式）");
block("BA", 29.4, 1.4, 14.15, 10.75, [29.45, 11.82], "Brønsted 酸分子筛（H-ZSM-5 / SAPO-34）");
block("OX", 37.7, 1.4, 5.85, 2.55, [41.55, 3.6], "氧化物（OEO）");
block("OX", 29.4, 8.45, 12.95, 0.75, [29.45, 9.25], "氧化物 H 溢流加氢 / 氢转移");
n("zMeOH", "ref", 30.2, 2.0, "CH_3OH"); n("MeOH2", "ion", 31.8, 2.0, "CH_3OH_2^{+}"); n("SMS", "surf", 33.5, 2.0, "Z–OCH_3"); n("DME", "mol", 35.2, 2.0, "DME");
n("TMO", "ion", 36.75, 2.0, "(CH_3)_3O^{+}"); n("acZ", "surf", 38.75, 2.0, "CH_3CO–Z"); n("ket", "surf", 40.55, 2.0, "CH_2=C=O"); n("zMeOAc", "ref", 42.45, 2.0, "CH_3COOCH_3");
n("zCO", "ref", 38.75, 3.3, "CO"); n("HCHO", "mol", 33.5, 3.3, "HCHO + CH_4"); n("MMC", "ion", 35.2, 3.3, "CH_3OCH_2^{+}");
n("zEtOH", "ref", 30.2, 4.6, "C_2H_5OH"); n("ethoxy", "surf", 32.0, 4.6, "Z–OC_2H_5 (≡C_2H_5^{+})");
n("zC2H4", "ref", 30.4, 6.0, "C_2H_4"); n("zC3H6", "ref", 33.0, 6.0, "C_3H_6"); n("zb1", "ref", 35.6, 6.0, "1-C_4H_8"); n("zb2", "ref", 37.2, 6.0, "2-C_4H_8（顺 / 反）", { fs: 8 });
n("zib", "ref", 38.8, 6.0, "i-C_4H_8"); n("zC5O", "ref", 41.2, 6.0, "C_{5+} 烯烃");
n("nC3p", "ion", 31.1, 7.4, "n-C_3H_7^{+}"); n("iC3p", "ion", 32.45, 7.4, "i-C_3H_7^{+}"); n("sC4p", "ion", 34.6, 7.4, "sec-C_4H_9^{+}");
n("iC4p", "ion", 36.3, 7.4, "i-C_4H_9^{+}"); n("tC4p", "ion", 37.75, 7.4, "t-C_4H_9^{+}"); n("C5p", "ion", 40.0, 7.4, "C_{5+}^{+}");
n("zC2H6", "ref", 30.4, 8.8, "C_2H_6"); n("zC3H8", "ref", 33.0, 8.8, "C_3H_8"); n("znC4", "ref", 35.6, 8.8, "n-C_4H_{10}"); n("ziC4", "ref", 38.8, 8.8, "i-C_4H_{10}");
n("zC5P", "ref", 41.2, 8.8, "C_{5+} 烷烃"); n("GAS", "mol", 42.05, 10.25, "C_5–C_{11} 汽油");
e("zMeOH", "MeOH2", "R4", "+H^{+}", { ly: -0.2 }); e("MeOH2", "SMS", "R4", "−H_2O", { ly: -0.2 }); e("SMS", "DME", "R4", "+CH_3OH", { both: true, ly: -0.2 });
e("DME", "TMO", "R4", "+SMS", { ly: -0.2 }); e("SMS", "HCHO", "R4", "+CH_3OH 氢转移", { perp: 0.2 }); e("DME", "MMC", "R4", "+SMS −CH_4", { perp: -0.2 });
e("TMO", "MMC", "R4", "−H^{+}", { perp: 0.15 }); e("zCO", "acZ", "R4", "+SMS（Koch）", { perp: 0.2 }); e("acZ", "ket", "R4", "−H^{+}", { both: true, ly: -0.2 });
e("ket", "zMeOAc", "R4", "+CH_3OH", { ly: -0.2 }); e("zEtOH", "ethoxy", "R4", "−H_2O", { ly: -0.2 }); e("ethoxy", "zC2H4", "R4", "−H^{+}", { perp: 0.15 });
e("MMC", "zC2H4", "R4", "第一个 C–C 键", { pts: [[35.2, 3.47], [35.2, 5.35], [30.55, 5.35], [30.55, 5.83]], lpos: [34.0, 5.35] });
e("ket", "zC2H4", "R4", "乙烯酮 → C_2H_4（OEO）", { pts: [[41.25, 2.0], [41.25, 2.0], [41.3, 5.15], [30.25, 5.15], [30.25, 5.83]], lpos: [39.6, 5.15], dash: "dash" });
e("zC2H4", "nC3p", "R4", "+CH_3^{+}", { perp: 0.15 }); e("nC3p", "iC3p", "SH", "H 迁移", { both: true, ly: 0.25 }); e("iC3p", "zC3H6", "R4", "−H^{+}", { both: true, perp: -0.15 });
e("zC3H6", "sC4p", "R4", "+CH_3^{+}", { perp: 0.15 }); e("zC3H6", "iC4p", "R4", "+CH_3^{+}（中碳）", { perp: 0.15, lt: 0.6 }); e("iC4p", "tC4p", "SH", "", { both: true });
e("sC4p", "zb1", "R4", "", { both: true }); e("sC4p", "zb2", "R4", "−H^{+}", { both: true, perp: 0.15, lt: 0.7 }); e("tC4p", "zib", "R4", "−H^{+}", { both: true, perp: -0.15 });
e("sC4p", "tC4p", "SH", "PCP 骨架异构", { both: true, pts: [[34.6, 7.61], [34.6, 8.12], [37.75, 8.12], [37.75, 7.61]], lpos: [36.2, 8.12] });
e("zb1", "C5p", "R4", "", { lt: 0.5 }); e("zb2", "C5p", "R4", "+CH_3^{+} / 齐聚", { perp: 0.15, lt: 0.55 }); e("zib", "C5p", "R4", "", {});
e("C5p", "zC5O", "R4", "−H^{+}", { both: true, perp: -0.15 });
e("C5p", "zC3H6", "R4", "β-断裂 → C_3^{=} / C_4^{=}", { pts: [[39.6, 7.21], [39.6, 5.45], [33.15, 5.45], [33.15, 5.83]], lpos: [38.4, 5.45] });
e("zC2H4", "zC2H6", "R3", "", {}); e("zC3H6", "zC3H8", "R3", "", {}); e("zb1", "znC4", "R3", "", {}); e("zib", "ziC4", "R3", "", {}); e("zC5O", "zC5P", "R3", "+H", { perp: 0.15 });
e("zC5P", "GAS", "R5", "", {});

// ================================================================== REGION AR: aromatics, all C6-C10 alkylbenzene isomers (x 29.3-43.6)
region(29.3, 12.45, 14.35, 11.2, C.accent5, "AR · 芳烃循环：C_6–C_{10} 全部烷基苯异构体 · 七甲基苯鎓^{+} · 萘 / 积碳");
block("BA", 29.4, 12.85, 14.15, 10.7, [37.9, 23.2], "Brønsted 酸分子筛：甲基化 · 烷基化 · 异构 · paring");
block("LA", 29.45, 12.95, 2.35, 2.05, [29.5, 15.08], "Ga / Zn（脱氢）");
n("cyc", "ion", 30.6, 13.35, "环烯基^{+} / 环二烯", { fs: 8 }); n("Bz", "mol", 30.6, 14.6, "苯");
n("Tol", "mol", 32.4, 14.6, "甲苯"); n("oX", "mol", 34.6, 13.6, "o-二甲苯"); n("mX", "mol", 34.6, 14.4, "m-二甲苯"); n("pX", "mol", 34.6, 15.2, "p-二甲苯"); n("EB", "mol", 34.6, 16.4, "乙苯");
const C9 = [["T123", "1,2,3-三甲苯"], ["T124", "1,2,4-三甲苯"], ["T135", "1,3,5-三甲苯"], ["oET", "o-乙基甲苯"], ["mET", "m-乙基甲苯"], ["pET", "p-乙基甲苯"], ["cum", "异丙苯"], ["nPB", "正丙苯"]];
C9.forEach(([id, l], i) => n(id, "mol", 37.15, 13.3 + 0.62 * i, l, { fs: 8 }));
const TE = [["Te1234", "1,2,3,4-四甲苯"], ["Te1235", "1,2,3,5-四甲苯"], ["Te1245", "1,2,4,5-四甲苯"]];
const DMEB = [["d34", "3-乙基-邻二甲苯"], ["d45", "4-乙基-邻二甲苯"], ["d2m", "2-乙基-间二甲苯"], ["d4m", "4-乙基-间二甲苯"], ["d5m", "5-乙基-间二甲苯"], ["d2p", "2-乙基-对二甲苯"]];
TE.forEach(([id, l], i) => n(id, "mol", 40.3, 13.25 + 0.5 * i, l, { fs: 8 }));
DMEB.forEach(([id, l], i) => n(id, "mol", 40.3, 15.2 + 0.5 * i, l, { fs: 8 }));
n("PMB", "mol", 42.85, 13.45, "五甲苯", { fs: 8 }); n("HMB", "mol", 42.85, 14.3, "六甲苯", { fs: 8 }); n("hepta", "ion", 42.85, 15.2, "七甲基苯鎓^{+}", { fs: 8 });
n("CP", "ion", 42.85, 16.45, "环戊烯基^{+}", { fs: 8 }); n("tCP", "txt", 42.85, 17.0, "→ C_3^{=} / i-C_4^{=}", { fs: 8, tc: C.accent4, bold: true });
n("side", "ion", 42.85, 18.0, "外亚甲基 / 乙基侧链^{+}", { fs: 7 }); n("tSide", "txt", 42.85, 18.55, "→ C_2H_4 / C_3H_6", { fs: 8, tc: C.accent4, bold: true });
const DEB = [["oDEB", "o-二乙苯"], ["mDEB", "m-二乙苯"], ["pDEB", "p-二乙苯"]];
DEB.forEach(([id, l], i) => n(id, "mol", 33.35 + 1.25 * i, 19.9, l, { fs: 8 }));
const CYM = [["oCy", "o-伞花烃"], ["mCy", "m-伞花烃"], ["pCy", "p-伞花烃"]], PT = [["oPT", "o-丙基甲苯"], ["mPT", "m-丙基甲苯"], ["pPT", "p-丙基甲苯"]];
CYM.forEach(([id, l], i) => n(id, "mol", 31.25 + 1.2 * i, 21.25, l, { fs: 8 }));
PT.forEach(([id, l], i) => n(id, "mol", 35.95 + 1.3 * i, 21.25, l, { fs: 8 }));
const BB = [["nBB", "正丁苯"], ["sBB", "仲丁苯"], ["iBB", "异丁苯"], ["tBB", "叔丁苯"]];
BB.forEach(([id, l], i) => n(id, "mol", 30.45 + 1.15 * i, 22.6, l, { fs: 8 }));
n("Naph", "mol", 40.3, 18.95, "萘", { fs: 8 }); n("MN1", "mol", 39.65, 19.85, "1-甲基萘", { fs: 8 }); n("MN2", "mol", 41.0, 19.85, "2-甲基萘", { fs: 8 }); n("PAH", "mol", 40.3, 20.85, "多环芳烃 / 积碳", { fs: 8 });
// family frames
const FAM = [["四甲苯", 40.3, 13.25, 14.25], ["二甲基乙苯（6）", 40.3, 15.2, 17.7], ["二乙苯", null], ["伞花烃 / 丙基甲苯", null], ["丁苯（4）", null]];
// edges: cyclization lane from C5+ cation
e("C5p", "cyc", "R4", "环化（C_6–C_8 碳正离子）", { pts: [[40.0, 7.61], [40.0, 7.85], [39.55, 7.85], [39.55, 12.93], [30.6, 12.93], [30.6, 13.17]], lpos: [34.4, 12.93] });
e("cyc", "Bz", "R4", "−3H：氢转移 / −H_2", { perp: -0.2 });
e("Bz", "Tol", "R5", "+CH_3^{+}", { ly: -0.2 }); e("Tol", "oX", "R5", "", {}); e("Tol", "mX", "R5", "", {}); e("Tol", "pX", "R5", "", {});
e("oX", "T123", "R5", "", {}); e("oX", "T124", "R5", "", {}); e("mX", "T123", "R5", "", {}); e("mX", "T124", "R5", "", {}); e("mX", "T135", "R5", "", {}); e("pX", "T124", "R5", "", {});
e("oX", "mX", "SH", "", { both: true, ax: -0.55, bx: -0.55 }); e("mX", "pX", "SH", "", { both: true, ax: -0.55, bx: -0.55 });
e("T123", "Te1234", "R5", "", {}); e("T123", "Te1235", "R5", "", {}); e("T124", "Te1234", "R5", "", {}); e("T124", "Te1235", "R5", "", {}); e("T124", "Te1245", "R5", "", {}); e("T135", "Te1235", "R5", "", {});
e("T123", "T124", "SH", "", { both: true, ax: -0.75, bx: -0.75 }); e("T124", "T135", "SH", "", { both: true, ax: -0.75, bx: -0.75 });
e("Te1234", "Te1235", "SH", "", { both: true, ax: -0.85, bx: -0.85 }); e("Te1235", "Te1245", "SH", "", { both: true, ax: -0.85, bx: -0.85 });
e("Te1234", "PMB", "R5", "", {}); e("Te1235", "PMB", "R5", "", {}); e("Te1245", "PMB", "R5", "", {});
e("PMB", "HMB", "R5", "", {}); e("HMB", "hepta", "R5", "+CH_3^{+}", { perp: -0.2 }); e("hepta", "CP", "R4", "paring", { perp: -0.25 }); e("CP", "tCP", "R4", "", {});
e("hepta", "side", "R4", "侧链", { pts: [[43.42, 15.2], [43.6, 15.2], [43.6, 18.0], [43.57, 18.0]], lpos: [43.3, 17.45] });
e("side", "tSide", "R4", "", {});
// alkylation (+C2H4 / +C3H6 / +C4H8)
e("Bz", "EB", "R5", "+C_2H_4", { pts: [[30.6, 14.77], [30.6, 16.4], [34.12, 16.4]], lpos: [32.4, 16.4], both: true });
e("Tol", "oET", "R5", "+C_2H_4", { pts: [[32.4, 14.77], [32.4, 15.8], [35.95, 15.8], [35.95, 15.16], [36.5, 15.16]], lpos: [33.4, 15.8] }); e("EB", "mET", "R5", "+CH_3^{+}", { perp: 0.15, lt: 0.6 }); e("EB", "pET", "R5", "", {});
e("EB", "oET", "R5", "", {});
e("Bz", "cum", "R5", "+C_3H_6", { pts: [[30.75, 14.77], [30.75, 17.02], [36.7, 17.02]], lpos: [33.0, 17.02] });
e("Bz", "nPB", "R5", "+C_3H_6（少）", { pts: [[30.45, 14.77], [30.45, 17.64], [36.7, 17.64]], lpos: [33.0, 17.64] });
// C10 families
e("oET", "d34", "R5", "+CH_3^{+}", { perp: 0.15, lt: 0.7 }); e("oET", "d45", "R5", "", {}); e("mET", "d2m", "R5", "", {}); e("mET", "d4m", "R5", "", {}); e("mET", "d5m", "R5", "", {}); e("pET", "d2p", "R5", "", {});
e("EB", "mDEB", "R5", "+C_2H_4", { pts: [[34.6, 16.57], [34.6, 19.73]], lpos: [34.6, 18.9] }); e("mDEB", "oDEB", "SH", "", { both: true }); e("mDEB", "pDEB", "SH", "", { both: true });
e("Tol", "mCy", "R5", "+C_3H_6", { pts: [[32.25, 14.77], [32.25, 21.08]], lpos: [32.25, 20.4] }); e("mCy", "oCy", "SH", "", { both: true }); e("mCy", "pCy", "SH", "", { both: true });
e("nPB", "mPT", "R5", "+CH_3^{+}", { pts: [[37.25, 17.81], [37.25, 21.08]], lpos: [37.25, 20.4] }); e("mPT", "oPT", "SH", "", { both: true }); e("mPT", "pPT", "SH", "", { both: true });
e("Bz", "sBB", "R5", "+C_4H_8", { pts: [[30.3, 14.77], [30.3, 22.0], [31.6, 22.0], [31.6, 22.43]], lpos: [30.95, 22.0] });
e("sBB", "nBB", "SH", "", {}); e("sBB", "iBB", "SH", "", {}); e("iBB", "tBB", "SH", "", { both: true });
e("Te1245", "Naph", "R4", "+烯烃 环化缩合", { pts: [[40.95, 14.25], [41.45, 14.25], [41.45, 18.95], [40.62, 18.95]], lpos: [41.45, 18.4] });
e("Naph", "MN1", "R5", "", {}); e("Naph", "MN2", "R5", "", {}); e("MN1", "PAH", "R4", "", {}); e("MN2", "PAH", "R4", "缩合", { perp: -0.15 });

// ================================================================== draw
const s = pres.addSlide({ masterName: "BIG" });
s.addText("CO₂ 一步热催化加氢：同分异构体分开的完整反应网络（C₁–C₄ 全显式 · C₆–C₁₀ 烷基苯全部异构体 · 含表面物种 / 自由基 / 碳正离子）", { placeholder: "title" });
REG.forEach((r) => {
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: r.x, y: r.y, w: r.w, h: r.h, rectRadius: 0.15, fill: { color: C.background1 },
    line: { color: r.color, width: 2.25, dashType: r.dash ? "dash" : "solid" }, objectName: `region ${r.title}` });
});
BLK.forEach((b, i) => s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: b.x, y: b.y, w: b.w, h: b.h, rectRadius: 0.1,
  fill: { color: CATAL[b.k].fill }, line: { color: CATAL[b.k].fill, width: 0.5 }, objectName: `catalyst block ${i + 1} ${b.k}` }));
// family frames in AR
[["四甲苯（3）", 39.45, 13.0, 1.75, 1.5], ["二甲基乙苯（6）", 39.35, 14.95, 1.95, 2.95], ["二乙苯（3）", 32.6, 19.62, 3.5, 0.56],
 ["伞花烃（3）", 30.55, 20.97, 3.4, 0.56], ["丙基甲苯（3）", 35.2, 20.97, 3.55, 0.56], ["丁苯（4）", 29.85, 22.32, 4.3, 0.56], ["萘系", 39.0, 18.68, 2.65, 2.45]]
  .forEach(([t, x, y, w, h], i) => {
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, rectRadius: 0.06, fill: { color: "FFFFFF", transparency: 55 },
      line: { color: "8B5A2B", width: 0.75, dashType: "sysDash" }, objectName: `family ${t}` });
  });
REG.forEach((r) => txt(s, r.x + 0.18, r.y + 0.06, r.w - 0.4, 0.32, runs(r.title, { fontSize: 13, bold: true, color: r.color }), { name: "region title" }));
E.forEach((ed, i) => drawEdge(s, ed, i));
Object.values(N).forEach((nd) => drawNode(s, nd));
// family labels
[["四甲苯（3）", 39.45, 12.8], ["二甲基乙苯（6）", 39.35, 14.75], ["二乙苯（3）", 32.6, 19.42], ["伞花烃（3）", 30.55, 20.77],
 ["丙基甲苯（3）", 35.2, 20.77], ["丁苯（4）", 33.2, 22.12], ["萘系", 39.0, 18.48]].forEach(([t, x, y]) =>
  txt(s, x, y - 0.02, 1.4, 0.22, t, { fs: 7, color: "8B5A2B", bold: true, name: "family label" }));
NOTES.forEach(({ x, y, w, h, t, o }) => {
  const r = []; t.forEach((line, k) => { const rr = runs(line, { fontSize: 7.5, color: o.color || C.text2 }); if (k < t.length - 1) rr[rr.length - 1].options.breakLine = true; r.push(...rr); });
  txt(s, x, y, w, h, r, { valign: "top", name: "note" });
});
// catalyst chips
BLK.forEach((b, i) => {
  if (!b.chip) return;
  const K = CATAL[b.k], t = b.lab || K.n, fs = 8, w = textW(t, fs) * 1.1 + 0.2, h = 0.22;
  s.addText(runs(t, { fontSize: fs, color: "FFFFFF", bold: true }), { x: b.chip[0], y: b.chip[1], w, h, shape: pres.shapes.ROUNDED_RECTANGLE, rectRadius: h / 2,
    fill: { color: K.c }, line: { color: K.c, width: 0.5 }, align: "center", valign: "middle", margin: 0, objectName: `catalyst label ${i + 1}` });
});
// legends (top-right strip + bottom)
{
  const y = 23.85;
  const types = [["mol", "分子 / 产物"], ["ref", "同一物种的引用副本"], ["surf", "表面吸附物种 *"], ["rad", "自由基 •"], ["ion", "离子 / 碳正离子 ⁺"]];
  let x = 0.4;
  types.forEach(([t, lab]) => { drawNode(s, { id: `lg-${t}`, type: t, x: x + 0.3, y: y + 0.17, label: "", w: 0.55, h: 0.24 });
    txt(s, x + 0.68, y + 0.03, 2.2, 0.28, lab, { fs: 10, color: C.text1 }); x += textW(lab, 10) + 1.05; });
  x += 0.3;
  Object.values(CAT).forEach((k) => { seg(s, x, y + 0.17, x + 0.45, y + 0.17, { color: k.c, width: 2, dash: k.dash, end: true, name: "legend line" });
    txt(s, x + 0.52, y + 0.03, 3.0, 0.28, k.n, { fs: 10, color: C.text1 }); x += textW(k.n, 10) + 0.85; });
  const y2 = 24.3; x = 0.4;
  txt(s, x, y2, 1.6, 0.3, "催化剂形式（色块）", { fs: 10, color: C.text1, bold: true }); x += 1.75;
  Object.values(CATAL).forEach((K) => { s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y: y2 + 0.04, w: 0.45, h: 0.22, rectRadius: 0.05,
    fill: { color: K.fill }, line: { color: K.c, width: 1 }, objectName: "catalyst legend swatch" });
    txt(s, x + 0.52, y2, 2.6, 0.3, K.s, { fs: 10, color: C.text1 }); x += textW(K.s, 10) + 0.85; });
  txt(s, 30.5, y2, 13.2, 0.3, "C₅₊ 脂肪烃按类集总；芳烃 C₆–C₁₀ 全部烷基苯构造异构体；未区分对映异构体", { fs: 10, color: C.text2, align: "right" });
}
s.addNotes("一页完整网络：A 区为 C1 表面中间体（氧化物上的甲酸盐 / 氧空位路线、RWGS、金属上的 CO 活化与甲烷化）；R 区为对应的气相自由基（C1 加氢链、C1–C4 烷基自由基的抽氢 / β-断裂 / 复合 / H 迁移）；F 区为金属碳化物上的链增长，C2–C4 每个烷基、烯烃、烷烃异构体单独画出，CO 插入得到各个醛、醇，以及双功能升级与烯烃齐聚到航煤；Z 区为分子筛上的甲醇活化、第一个 C–C 键和 C2–C4 碳正离子；AR 区为芳烃循环，C6–C10 所有烷基苯构造异构体（36 个）及其甲基化、烷基化、异构化关系。灰色虚线框节点是同一物种在其他区的引用副本。");

(async () => { await pres.writeFile({ fileName: OUT }); await applyTheme(OUT, THEME);
  const cnt = {}; Object.values(N).forEach((x) => { cnt[x.type] = (cnt[x.type] || 0) + 1; });
  console.log("wrote", OUT, "nodes", Object.keys(N).length, JSON.stringify(cnt), "edges", E.length); })();
