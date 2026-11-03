// Complete CO2/C1 catalytic reaction network incl. radicals (no NH3/urea part)
const pptxgen = require("pptxgenjs");
const { applyTheme } = require("/root/.claude/skills/synced/62cc770a-8c2c-45f8-8ae0-ca258c8631d4_456331ad-9fc1-4922-9170-b4c06151a9b5/pptx/scripts/apply_theme.js");
const OUT = process.argv[2] || "full_network.pptx";

const THEME = {
  name: "C1 Network", headFontFace: "Microsoft YaHei", bodyFontFace: "Microsoft YaHei",
  colors: { dk1: "14232B", lt1: "FFFFFF", dk2: "3B4A54", lt2: "F3F5F6", accent1: "E07B24", accent2: "D63B3B",
            accent3: "2A78D6", accent4: "1E9E5A", accent5: "7B4FC4", accent6: "8B5A2B", hlink: "2A78D6", folHlink: "7B4FC4" },
};
const pres = new pptxgen();
pres.defineLayout({ name: "POSTER", width: 20, height: 12.5 });
pres.layout = "POSTER";
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
pres.title = "CO2 加氢及 C1 催化完整反应网络（含自由基）";
const C = pres.SchemeColor;
const RAD = "C2185B";   // radical category colour (hex, also used where scheme colours are not accepted)
const CAT = {
  R1: { c: C.accent1, hex: "E07B24", n: "① 重整 / 合成气" },
  R2: { c: C.accent2, hex: "D63B3B", n: "② 费托 / C–C 链增长" },
  R3: { c: C.accent3, hex: "2A78D6", n: "③ 加氢还原" },
  R4: { c: C.accent4, hex: "1E9E5A", n: "④ 酸催化 / 分子筛" },
  R5: { c: C.accent5, hex: "7B4FC4", n: "⑤ 羰基化" },
  R6: { c: C.accent6, hex: "8B5A2B", n: "⑥ 烃类升级" },
  R7: { c: RAD, hex: RAD, n: "⑦ 自由基链反应", dash: "dash" },
  SH: { c: C.text2, hex: "3B4A54", n: "多类反应共享步骤" },
};
const CAT_ORDER = ["R1", "R2", "R3", "R4", "R5", "R6", "R7"];

pres.defineSlideMaster({
  title: "POSTER", background: { color: C.background1 },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.5, y: 0.12, w: 19, h: 0.6, fontSize: 28, bold: true,
        color: C.text1, align: "center", valign: "middle", margin: 0 }, text: "" } },
  ],
});
pres.defineSlideMaster({
  title: "TABLE", background: { color: C.background1 },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.6, y: 0.35, w: 18.8, h: 0.9, fontSize: 32, bold: true,
        color: C.text1, align: "left", valign: "middle", margin: 0 }, text: "" } },
    { text: { text: "CO₂ 加氢及 C1 催化完整反应网络（含自由基）", options: { x: 0.6, y: 11.9, w: 10, h: 0.35, fontSize: 12, color: C.text2, margin: 0 } } },
  ],
  slideNumber: { x: 18.9, y: 11.9, w: 0.5, h: 0.35, fontSize: 12, color: C.text2, align: "right" },
});

// ---------------------------------------------------------------- text helpers
// "C_nH_{2n+1}^+" style: _x / _{..} subscript, ^{..} superscript
function runs(s, opt = {}) {
  const out = [];
  const re = /_\{([^}]*)\}|_(\d+|[a-z])|\^\{([^}]*)\}/g;
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

// ---------------------------------------------------------------- network data
// type: mol | surf | rad | ion | hub | txt
const N = {};
function n(id, type, cx, cy, label, w, extra = {}) {
  const dw = { mol: 1.05, surf: 1.05, rad: 1.0, ion: 1.9, hub: 1.7, txt: 1.0 }[type];
  N[id] = { id, type, cx, cy, label, w: w || dw, h: type === "rad" ? 0.42 : 0.44, ...extra };
}
// feeds
n("fCO2", "mol", 0.95, 3.05, "CO_2"); n("fH2", "mol", 0.95, 4.35, "H_2"); n("fH2O", "mol", 0.95, 5.4, "H_2O");
n("fO2", "mol", 0.95, 6.45, "O_2"); n("fCH4", "mol", 0.95, 7.55, "CH_4");
// electro
n("eCO2m", "rad", 3.0, 1.75, "CO_2•^{−}", 1.05); n("mHCOOH", "mol", 4.6, 1.75, "HCOOH", 1.15);
// C1 surface core
n("sCO2", "surf", 3.35, 3.05, "CO_2*"); n("sHCOO", "surf", 4.6, 3.05, "HCOO*"); n("sH2COO", "surf", 5.85, 3.05, "H_2COO*", 1.1);
n("sH2COOH", "surf", 7.08, 3.05, "H_2COOH*", 1.2); n("sH2CO", "surf", 8.25, 3.05, "H_2CO*", 0.95); n("sCH3O", "surf", 9.25, 3.05, "CH_3O*", 0.9);
n("sH", "surf", 3.0, 4.35, "H*", 0.8); n("sOH", "surf", 3.0, 5.4, "OH*", 0.8); n("sO", "surf", 3.0, 6.45, "O*", 0.8);
n("sCOOH", "surf", 4.6, 4.35, "COOH*"); n("hubCO", "hub", 6.45, 4.25, "CO(g) 合成气"); n("sHCO", "surf", 8.8, 4.35, "HCO*", 0.9);
n("sCO", "surf", 7.9, 5.5, "CO*", 0.9); n("sCOH", "surf", 6.55, 6.3, "COH*", 0.9);
n("sC", "surf", 7.9, 7.45, "C*", 0.8); n("sCH", "surf", 6.6, 7.45, "CH*", 0.8); n("sCH2", "surf", 5.3, 7.45, "CH_2*", 0.85);
n("sCH3", "surf", 4.0, 7.45, "CH_3*", 0.85);
// radicals (gas phase / homogeneous)
n("tCO2H", "txt", 2.85, 9.45, "CO_2 + H•", 0.95);
n("rHOCO", "rad", 3.95, 9.1, "HOCO•"); n("rHCOO", "rad", 3.95, 9.8, "HCOO•");
n("rHCO", "rad", 5.25, 9.1, "HCO•"); n("rCH2O", "mol", 6.45, 9.1, "CH_2O", 0.9);
n("rCH3O", "rad", 7.65, 9.1, "CH_3O•"); n("rCH2OH", "rad", 7.65, 9.8, "•CH_2OH");
n("rMeOH", "mol", 9.0, 9.45, "CH_3OH", 1.0); n("rHCOOH", "txt", 5.25, 9.8, "→ HCOOH", 1.0);
n("rCH3", "rad", 3.1, 10.5, "•CH_3"); n("rC2", "mol", 4.55, 10.5, "C_2H_6/C_2H_4", 1.5);
n("rCH3OO", "rad", 6.45, 10.5, "CH_3OO•", 1.05);
n("rH", "rad", 7.6, 10.5, "H•", 0.8); n("rHO2", "rad", 8.5, 10.5, "HO_2•", 0.85); n("rOH", "rad", 9.35, 10.5, "OH•", 0.8);
// carbonylation
n("mMeONO", "mol", 11.4, 1.97, "CH_3ONO", 1.15); n("rMeONO", "rad", 13.0, 1.97, "CH_3O• + NO•", 1.6);
n("sPd", "surf", 14.7, 1.97, "Pd–COOCH_3", 1.45); n("mDMO", "mol", 16.2, 1.97, "DMO", 0.85);
n("mMG", "mol", 17.35, 1.97, "MG", 0.8); n("mEG", "mol", 18.5, 1.97, "EG", 0.8);
// zeolite acid sites
n("mMeOH", "mol", 10.75, 3.05, "CH_3OH", 1.05); n("zOCH3", "surf", 12.25, 3.05, "Z–OCH_3 表面甲氧基", 1.85);
n("mDME", "mol", 13.95, 3.05, "DME", 0.85); n("mARO", "mol", 15.65, 3.05, "芳烃 BTX / C_9–C_{10}", 1.95);
n("mGAS", "mol", 17.95, 5.95, "C_5–C_{11} 汽油", 1.6);
n("zKET", "surf", 11.3, 4.6, "CH_2=C=O 乙烯酮", 1.65); n("iPOOL", "ion", 14.0, 4.6, "烃池：多甲基苯^{+} / 环戊烯基^{+}", 2.75);
n("mC24", "mol", 17.95, 4.6, "C_2–C_4 烯烃", 1.45);
n("zOEt", "surf", 12.55, 5.95, "Z–OC_2H_5 乙氧基", 1.75); n("mC5O", "mol", 15.6, 5.95, "C_{5+} 烯烃", 1.25);
// metal C-C chain growth + acid upgrading
n("mEtOH", "mol", 11.0, 7.52, "C_2H_5OH", 1.1); n("sCH3CO", "surf", 11.0, 8.24, "CH_3CO*", 1.05);
n("sCnH", "surf", 13.35, 9.15, "C_nH_{2n+1}* 烷基链", 1.85);
n("mAOle", "mol", 15.6, 7.6, "α-烯烃 C_nH_{2n}（FTO）", 2.05); n("mNPar", "mol", 15.75, 9.15, "正构烷烃 C_nH_{2n+2}", 1.85);
n("mWax", "mol", 15.75, 10.4, "C_{20+} 蜡", 1.25); n("mHA", "mol", 12.3, 10.4, "C_{2+} 醇（CO 插入）", 1.85);
n("iCARB", "ion", 18.0, 8.3, "C_nH_{2n+1}^{+} 碳正离子", 2.15); n("mJET", "mol", 18.0, 10.2, "C_8–C_{16} 航煤 / 柴油", 2.05);

// edges: [from, to, cat, label, opts]
const E = [];
const e = (a, b, cat, label, o = {}) => E.push({ a, b, cat, label, ...o });
// feeds
e("fCO2", "sCO2", "SH", "吸附", { lt: 0.45, lw: 0.5 });
e("fH2", "sH", "SH", "解离", { lt: 0.45, lw: 0.5 });
e("fH2O", "sOH", "R1", "SMR", { both: true, lt: 0.45, lw: 0.5 });
e("fO2", "sO", "R1", "POM", { lt: 0.45, lw: 0.5 });
e("fCH4", "sCH3", "SH", "C–H 活化 ⇌ 甲烷化", { both: true, lt: 0.4, lw: 1.55 });
e("sO", "sOH", "SH", "+H*", { perp: 0.25 });
e("sOH", "sH", "SH", "", {});
// electro
e("fCO2", "eCO2m", "R3", "+e^{−}（电还原）", { lt: 0.45, lw: 1.35, perp: 0.25 });
e("eCO2m", "mHCOOH", "R3", "+H^{+} + e^{−}", { lw: 0.95, ly: -0.27 });
e("mHCOOH", "sHCOO", "R3", "", { both: true });
// formate route (methanol): labels replaced by one annotation
e("sCO2", "sHCOO", "R3", ""); e("sHCOO", "sH2COO", "R3", ""); e("sH2COO", "sH2COOH", "R3", "");
e("sH2COOH", "sH2CO", "R3", ""); e("sH2CO", "sCH3O", "R3", ""); e("sCH3O", "mMeOH", "R3", "");
// carboxyl / RWGS / CO
e("sCO2", "sCOOH", "R3", "+H*", { perp: -0.22, lw: 0.5 });
e("sCOOH", "sCO", "SH", "−OH*（RWGS）", { lt: 0.42, lw: 1.3, perp: -0.2 });
e("hubCO", "sCO", "SH", "", { both: true });
e("sCO", "sHCO", "R3", "+H*", { perp: -0.2 }); e("sHCO", "sH2CO", "R3", "+H*", { perp: -0.2 });
e("sCO", "sCOH", "R3", "+H*", { perp: 0.2 }); e("sCOH", "sC", "SH", "H 助 C–O 解离", { lt: 0.55, lw: 1.5, perp: -0.2 });
e("sCO", "sC", "SH", "C* + O* ⇌ CO*", { both: true, lw: 1.4, perp: -0.75 });
e("sC", "sCH", "R3", "", { both: true }); e("sCH", "sCH2", "R3", "", { both: true }); e("sCH2", "sCH3", "R3", "", { both: true });
// chain growth links (lanes between the core and the radical zone)
e("sCH2", "sCnH", "R2", "CH_2* 插入 → 链增长（F–T / FTO）", { pts: [[5.3, 7.67], [5.3, 8.08], [10.0, 8.08], [10.0, 9.15], [12.4, 9.15]],
  lpos: [6.55, 8.08], lw: 2.3 });
e("sCH3", "sCH3CO", "R2", "+CO* 插入（还原偶联 / 同系化）", { pts: [[4.0, 7.67], [4.0, 8.24], [10.45, 8.24]],
  lpos: [8.7, 8.24], lw: 2.25 });
// methanol -> zeolite
e("mMeOH", "zOCH3", "R4", "脱水", { perp: 0.21, lw: 0.5 });
e("zOCH3", "mDME", "R4", "", {});
e("zOCH3", "iPOOL", "R4", "甲基化", { perp: 0.22, lt: 0.45 });
e("mDME", "iPOOL", "R4", "MTH", { perp: -0.22 });
e("sCO", "zKET", "R4", "OEO：CO* + CH_3* → 乙烯酮", { dash: "dash", lt: 0.55, lw: 2.25, perp: -0.2 });
e("zKET", "mC24", "R4", "乙烯酮 → 乙烯", { pts: [[11.3, 4.82], [11.3, 5.35], [17.45, 5.35], [17.45, 4.82]], lpos: [12.9, 5.35], dash: "dash", lw: 1.25 });
e("iPOOL", "mC24", "R4", "侧链消除 / 脱烷基", { lt: 0.5, lw: 1.65, perp: 0.2 });
e("iPOOL", "mARO", "R4", "氢转移 / 脱氢", { lw: 1.3, perp: 0.22 });
e("mC24", "mC5O", "R4", "甲基化 / 齐聚 ⇄ 裂解", { both: true, lw: 1.85, perp: 0.3, lt: 0.35 });
e("mC5O", "iPOOL", "R4", "环化", { perp: -0.2, lt: 0.4 });
e("mC5O", "mGAS", "R3", "+H（氧化物加氢）", { ly: -0.28, lw: 1.45 });
e("mEtOH", "zOEt", "R4", "脱水", { perp: 0.22, lw: 0.5 });
e("zOEt", "mC24", "R4", "", {});
// carbonylation (labels below the edges)
e("mMeOH", "mMeONO", "R5", "+NO + ½O_2", { lw: 0.95, lx: -0.5, ly: 0.1 });
e("mMeONO", "rMeONO", "R5", "Pd 上 O–N 断裂", { lw: 1.35, ly: -0.31 });
e("rMeONO", "sPd", "R5", "+CO", { ly: -0.31, lw: 0.5 });
e("sPd", "mDMO", "R5", "偶联", { ly: -0.31, lw: 0.5 });
e("mDMO", "mMG", "R3", "+H_2", { ly: -0.31, lw: 0.5 }); e("mMG", "mEG", "R3", "+H_2", { ly: -0.31, lw: 0.5 });
// metal chain growth
e("sCH3CO", "mEtOH", "R3", "+4H*（经 CH_3CHO）", { lw: 1.7, perp: -1.4 });
e("sCH3CO", "sCnH", "R2", "链增长", { perp: 0.2 });
e("sCnH", "mAOle", "R2", "β-H 消除", { perp: -0.25 });
e("sCnH", "mNPar", "R2", "+H* 终止", { ly: -0.25 });
e("sCnH", "mWax", "R2", "链增长", { perp: -0.2, lt: 0.6 });
e("sCnH", "mHA", "R2", "CO* 插入", { perp: 0.22 });
e("mAOle", "iCARB", "R6", "质子化", { perp: 0.2 });
e("mWax", "iCARB", "R6", "金属脱氢 → 质子化", { perp: -0.22, lt: 0.45, lw: 1.65 });
e("iCARB", "mJET", "R6", "β-断裂 / 异构 / 加氢", { lw: 1.75, perp: 0.25 });
e("iCARB", "mC24", "R6", "裂解 / 低聚", { pts: [[19.08, 8.3], [19.35, 8.3], [19.35, 4.6], [18.68, 4.6]], lpos: [18.85, 7.55], lw: 1.0 });
// radicals
e("fCH4", "rCH3", "R7", "高温均裂 / +O_2", { pts: [[0.95, 7.77], [0.95, 10.5], [2.6, 10.5]], lpos: [1.75, 10.25], lw: 1.4 });
e("tCO2H", "rHOCO", "R7", ""); e("tCO2H", "rHCOO", "R7", "");
e("rHOCO", "rHCO", "R7", "−OH• → CO；+H•", { ly: 0.27, lw: 1.5 });
e("rHCOO", "rHCOOH", "R7", "", {});
e("rHCO", "rCH2O", "R7", "+H•", { ly: 0.27, lw: 0.5 });
e("rCH2O", "rCH3O", "R7", "+H•", { ly: 0.27, lw: 0.5 }); e("rCH2O", "rCH2OH", "R7", "+H•", { perp: -0.22, lw: 0.5 });
e("rCH3O", "rMeOH", "R7", "+H•", { perp: 0.2, lw: 0.5 }); e("rCH2OH", "rMeOH", "R7", "", {});
e("rCH3O", "rCH2OH", "R7", "", { both: true });
e("rCH3", "rC2", "R7", "OCM", { ly: -0.25, lw: 0.5 });
e("rCH3", "rCH3OO", "R7", "+O_2", { pts: [[3.1, 10.71], [3.1, 10.88], [6.45, 10.88], [6.45, 10.71]], lpos: [5.75, 10.88], lw: 0.45 });
e("rCH3OO", "rCH2O", "R7", "−OH•", { perp: 0.25, lw: 0.5 });
e("rH", "rHO2", "R7", "+O_2", { ly: -0.25, lw: 0.5 });
e("rHO2", "rOH", "R7", "+H•", { ly: -0.25, lw: 0.5 });
e("rOH", "rH", "R7", "+H_2 → H_2O", { pts: [[9.35, 10.71], [9.35, 10.88], [7.6, 10.88], [7.6, 10.71]], lpos: [8.5, 10.88], lw: 1.1 });

// categories passing through each node (by chemistry; dots and slide 3)
const NODE_CATS = {
  fCO2: "R1 R3", fH2: "R2 R3", fH2O: "R1", fO2: "R1 R7", fCH4: "R1 R3 R7", eCO2m: "R3 R7", mHCOOH: "R3",
  sCO2: "R1 R3", sHCOO: "R3", sH2COO: "R3", sH2COOH: "R3", sH2CO: "R3", sCH3O: "R3 R5", sH: "R1 R2 R3", sOH: "R1 R3", sO: "R1 R3",
  sCOOH: "R3", hubCO: "R1 R2 R3 R4 R5", sHCO: "R2 R3", sCO: "R1 R2 R3 R4", sCOH: "R2 R3", sC: "R1 R2 R3", sCH: "R1 R2 R3",
  sCH2: "R2 R3", sCH3: "R1 R2 R3",
  rHOCO: "R7", rHCOO: "R7", rHCO: "R7", rCH2O: "R7", rCH3O: "R7", rCH2OH: "R7", rMeOH: "R7", rCH3: "R1 R7", rC2: "R7",
  rCH3OO: "R7", rH: "R7", rHO2: "R7", rOH: "R7",
  mMeONO: "R5", rMeONO: "R5 R7", sPd: "R5", mDMO: "R3 R5", mMG: "R3", mEG: "R3",
  mMeOH: "R3 R4 R5", zOCH3: "R4", mDME: "R4", mARO: "R4 R6", mGAS: "R4 R6", zKET: "R4", iPOOL: "R4", mC24: "R2 R4 R6",
  zOEt: "R4", mC5O: "R4 R6", mEtOH: "R2 R3 R4", sCH3CO: "R2 R3", sCnH: "R2", mAOle: "R2 R6", mNPar: "R2", mWax: "R2 R6",
  mHA: "R2 R3", iCARB: "R4 R6", mJET: "R6",
};

// ---------------------------------------------------------------- drawing helpers
function clip(nd, ux, uy) {
  if (nd.type === "rad") {
    const a = nd.w / 2, b = nd.h / 2, t = 1 / Math.sqrt((ux * ux) / (a * a) + (uy * uy) / (b * b));
    return [nd.cx + ux * t, nd.cy + uy * t];
  }
  if (nd.type === "txt") { const t = Math.min(nd.w / 2 / Math.max(Math.abs(ux), 1e-9), 0.16 / Math.max(Math.abs(uy), 1e-9)); return [nd.cx + ux * t, nd.cy + uy * t]; }
  const tx = ux ? nd.w / 2 / Math.abs(ux) : Infinity, ty = uy ? nd.h / 2 / Math.abs(uy) : Infinity;
  const t = Math.min(tx, ty);
  return [nd.cx + ux * t, nd.cy + uy * t];
}
function seg(slide, x1, y1, x2, y2, o) {
  slide.addShape(pres.shapes.LINE, {
    x: Math.min(x1, x2), y: Math.min(y1, y2), w: Math.max(Math.abs(x2 - x1), 0.001), h: Math.max(Math.abs(y2 - y1), 0.001),
    flipH: x2 < x1, flipV: y2 < y1,
    line: { color: o.color, width: o.width || 1.5, dashType: o.dash || "solid",
            endArrowType: o.end ? "triangle" : undefined, beginArrowType: o.begin ? "triangle" : undefined },
    objectName: o.name,
  });
}
function label(slide, x, y, text, color, o = {}) {
  const lw = o.lw || 0.75, lh = 0.26;
  slide.addText(runs(text, { fontSize: o.fs || 10, color }), { x: x - lw / 2, y: y - lh / 2, w: lw, h: lh, align: "center",
    valign: "middle", margin: 0, isTextBox: true, fill: { color: C.background1 }, objectName: o.name || "edge label" });
}
function drawEdge(slide, ed, i) {
  const A = N[ed.a], B = N[ed.b], cat = CAT[ed.cat];
  const st = { color: cat.c, width: ed.width || 1.6, dash: ed.dash || cat.dash };
  const nm = `edge ${i + 1} ${ed.a}→${ed.b}`;
  if (ed.pts) {
    const P = ed.pts;
    for (let k = 0; k < P.length - 1; k++)
      seg(slide, P[k][0], P[k][1], P[k + 1][0], P[k + 1][1], { ...st, end: k === P.length - 2, begin: ed.both && k === 0, name: `${nm} seg ${k + 1}` });
    if (ed.label) label(slide, ed.lpos[0], ed.lpos[1], ed.label, cat.c, { lw: ed.lw, name: `${nm} label` });
    return;
  }
  const ax = A.cx + (ed.ax || 0), ay = A.cy + (ed.ay || 0), bx = B.cx + (ed.bx || 0), by = B.cy + (ed.by || 0);
  const dx = bx - ax, dy = by - ay, L = Math.hypot(dx, dy), ux = dx / L, uy = dy / L;
  let [x1, y1] = clip({ ...A, cx: ax, cy: ay }, ux, uy);
  let [x2, y2] = clip({ ...B, cx: bx, cy: by }, -ux, -uy);
  const g = 0.03;
  x1 += ux * g; y1 += uy * g; x2 -= ux * g; y2 -= uy * g;
  seg(slide, x1, y1, x2, y2, { ...st, end: true, begin: !!ed.both, name: nm });
  if (ed.label) {
    let nx = -uy, ny = ux;
    if (ny > 0 || (Math.abs(ny) < 1e-9 && nx > 0)) { nx = -nx; ny = -ny; }
    const t = ed.lt ?? 0.5, pp = ed.perp ?? 0;
    label(slide, x1 + (x2 - x1) * t + nx * pp + (ed.lx || 0), y1 + (y2 - y1) * t + ny * pp + (ed.ly || 0), ed.label, cat.c,
      { lw: ed.lw, name: `${nm} label` });
  }
}
// categories passing through each node (dots)
const nodeCats = {};
for (const [id, cs] of Object.entries(NODE_CATS)) nodeCats[id] = new Set(cs.split(" "));
function drawNode(slide, nd) {
  const x = nd.cx - nd.w / 2, y = nd.cy - nd.h / 2;
  const t = nd.type;
  if (t === "txt") {
    slide.addText(runs(nd.label, { fontSize: 11, color: RAD, bold: true }), { x, y, w: nd.w, h: nd.h, align: "center", valign: "middle",
      margin: 0, isTextBox: true, objectName: `node ${nd.id}` });
    return;
  }
  const style = {
    mol: { shape: pres.shapes.ROUNDED_RECTANGLE, fill: C.text2, line: C.text2, tc: C.background1, bold: true, r: 0.06 },
    surf: { shape: pres.shapes.ROUNDED_RECTANGLE, fill: C.background1, line: C.text1, tc: C.text1, bold: true, r: 0.22 },
    rad: { shape: pres.shapes.OVAL, fill: C.background1, line: RAD, tc: C.text1, bold: true, dash: "dash" },
    ion: { shape: pres.shapes.RECTANGLE, fill: C.background1, line: C.text1, tc: C.text1, bold: true },
    hub: { shape: pres.shapes.ROUNDED_RECTANGLE, fill: C.accent1, line: C.accent1, tc: C.background1, bold: true, r: 0.08 },
  }[t];
  if (t === "ion") slide.addShape(pres.shapes.RECTANGLE, { x: x - 0.05, y: y - 0.05, w: nd.w + 0.1, h: nd.h + 0.1,
    fill: { color: C.background1 }, line: { color: C.text1, width: 1 }, objectName: `node ${nd.id} outer border` });
  const opt = { x, y, w: nd.w, h: nd.h, shape: style.shape, fill: { color: style.fill },
    line: { color: style.line, width: t === "rad" ? 1.5 : 1.25, dashType: style.dash || "solid" },
    align: "center", valign: "middle", margin: 0.02, objectName: `node ${nd.id}` };
  if (style.r !== undefined) opt.rectRadius = style.r;
  slide.addText(runs(nd.label, { fontSize: t === "hub" ? 13 : 12, color: style.tc, bold: style.bold }), opt);
  // category dots
  const cats = CAT_ORDER.filter((k) => nodeCats[nd.id] && nodeCats[nd.id].has(k));
  // a lone ⑦ dot on a radical-zone node carries no information (the zone already says so)
  if (nd.noDots || !cats.length || (cats.length === 1 && cats[0] === "R7")) return;
  const d = 0.085, gap = 0.035, tot = cats.length * d + (cats.length - 1) * gap;
  cats.forEach((k, i) => slide.addShape(pres.shapes.OVAL, { x: nd.cx + nd.w / 2 - tot - 0.02 + i * (d + gap),
    y: nd.cy + nd.h / 2 + (t === "ion" ? 0.09 : 0.04), w: d, h: d, fill: { color: CAT[k].c }, line: { color: CAT[k].c, width: 0.5 },
    objectName: `node ${nd.id} dot ${k}` }));
}
function zone(slide, x, y, w, h, color, title, o = {}) {
  slide.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, rectRadius: 0.12, fill: { color: o.fill || C.background2 },
    line: { color, width: 1.75, dashType: o.dash || "solid" }, objectName: `zone ${title}` });
  slide.addText(runs(title, { fontSize: 12, bold: true, color }), { x: o.tx ?? x + 0.15, y: o.ty ?? y + 0.06, w: o.tw ?? w - 0.3, h: 0.3,
    align: o.talign || "left", margin: 0, isTextBox: true, objectName: `zone title ${title}` });
}

// ================================================================= slide 1: full network
let s = pres.addSlide({ masterName: "POSTER" });
s.addText("CO₂ 加氢及 C1 催化完整反应网络 · 中间体与自由基", { placeholder: "title" });
s.addText("节点形状 = 物种类型；连线颜色 = 反应类别；深灰 = 多类反应共享步骤；节点右下圆点 = 经过该中间体的反应类别。表面吸附态 R* 与气相自由基 R• 一一对应（上下两区）。", {
  x: 0.5, y: 0.74, w: 19, h: 0.3, fontSize: 11, color: C.text2, align: "center", margin: 0, isTextBox: true });
zone(s, 2.45, 2.32, 7.35, 5.68, CAT.R1.c, "C_1 表面中间体核心（金属 / 金属–氧化物界面）", { tx: 4.95, tw: 4.7, talign: "right" });
s.addText(runs("甲酸盐路线：每步 +H*（H_2COOH* → H_2CO* 为 −OH*）", { fontSize: 10, color: CAT.R3.c }), { x: 4.25, y: 3.42, w: 3.95, h: 0.26,
  align: "center", margin: 0, isTextBox: true, objectName: "formate note" });

zone(s, 2.45, 8.42, 7.35, 2.63, RAD, "⑦ 气相 / 均相自由基（POM、OCM、高温）：R• 与上方 R* 一一对应", { dash: "dash", fill: "FFF7F9" });
zone(s, 10.1, 1.1, 9.5, 1.3, CAT.R5.c, "⑤ 羰基化（亚硝酸甲酯路线，Pd）→ 乙二醇", { ty: 1.13 });
zone(s, 10.1, 2.48, 9.5, 4.17, CAT.R4.c, "④ 分子筛酸位：表面甲氧基 / 烃池 / 碳正离子（OX-ZEO 甲醇中间体路线）", { ty: 2.5, tx: 11.5, tw: 7.95 });
zone(s, 10.1, 6.88, 9.5, 4.17, CAT.R2.c, "② 金属位 C–C 链增长（Fe_5C_2 / Co_2C，费托 / FTO）· ⑥ 酸位烃类升级", { ty: 6.92, tx: 11.75, tw: 7.7 });
// pull HCOOH, CO2•- into the top band label
s.addText("电化学 / 光化学", { x: 0.45, y: 1.62, w: 1.9, h: 0.28, fontSize: 11, bold: true, color: CAT.R3.c, margin: 0, isTextBox: true });
E.forEach((ed, i) => drawEdge(s, ed, i));
Object.values(N).forEach((nd) => drawNode(s, nd));

// legend
const LY = 11.2;
s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 0.45, y: LY, w: 19.15, h: 1.1, rectRadius: 0.1, fill: { color: C.background2 },
  line: { color: C.background2 }, objectName: "legend panel" });
s.addText("物种类型", { x: 0.65, y: LY + 0.08, w: 1.5, h: 0.28, fontSize: 12, bold: true, color: C.text1, margin: 0, isTextBox: true });
const LT = [["mol", "原料 / 稳定分子 / 产物"], ["surf", "表面吸附物种 *"], ["rad", "自由基 •"], ["ion", "离子 / 碳正离子 ⁺"], ["hub", "合成气枢纽"]];
LT.forEach(([t, lab], i) => {
  const x = 0.75 + (i % 3) * 2.55, y = LY + 0.5 + Math.floor(i / 3) * 0.33;
  drawNode(s, { id: `legend-${t}`, type: t, cx: x + 0.3, cy: y, w: 0.55, h: 0.24, label: "", noDots: true });
  s.addText(lab, { x: x + 0.7, y: y - 0.13, w: 1.8, h: 0.26, fontSize: 10, color: C.text1, margin: 0, isTextBox: true });
});
s.addText("反应类别（连线颜色 / 节点圆点）", { x: 8.25, y: LY + 0.08, w: 4, h: 0.28, fontSize: 12, bold: true, color: C.text1, margin: 0, isTextBox: true });
[...CAT_ORDER, "SH"].forEach((k, i) => {
  const x = 8.3 + (i % 4) * 2.25, y = LY + 0.5 + Math.floor(i / 4) * 0.33;
  seg(s, x, y, x + 0.45, y, { color: CAT[k].c, width: 2, dash: CAT[k].dash, end: true, name: `legend line ${k}` });
  s.addText(CAT[k].n, { x: x + 0.55, y: y - 0.13, w: 1.8, h: 0.26, fontSize: 10, color: C.text1, margin: 0, isTextBox: true });
});
s.addText([
  { text: "DMO 草酸二甲酯 · MG 乙醇酸甲酯", options: { breakLine: true } },
  { text: "EG 乙二醇 · Z 分子筛骨架", options: { breakLine: true } },
  { text: "OEO 氧化物–乙烯酮–分子筛", options: { breakLine: true } },
  { text: "FTO 费托制烯烃", options: { breakLine: true } },
  { text: "不含氨 / 尿素（C–N）部分" },
], { x: 17.35, y: LY + 0.08, w: 2.2, h: 0.98, fontSize: 9, color: C.text2, margin: 0, valign: "top", isTextBox: true });
s.addNotes("完整网络图：左上为电化学还原 CO2；左侧为原料；中部橙框为 C1 表面中间体核心（甲酸盐路线 CO2*→HCOO*→H2COO*→H2COOH*→H2CO*→CH3O*→CH3OH；羧基/RWGS 路线 CO2*→COOH*→CO*；CO* 加氢与解离到 CHx*，甲烷化与重整）；左下红色虚线框为与之对应的气相自由基网络（HOCO•、HCOO•、HCO•、CH3O•、•CH2OH，以及 POM/OCM 的 •CH3、CH3OO•、H•、HO2•、OH•）。右侧：紫框为亚硝酸甲酯羰基化到乙二醇；绿框为分子筛（OX-ZEO）甲醇到烃，包括表面甲氧基、乙烯酮、烃池、烯烃循环和芳烃；红框为金属位链增长（费托/FTO，CO 插入到醇）和酸位升级到航煤/柴油。");

// ================================================================= slide 2: radical inventory
s = pres.addSlide({ masterName: "TABLE" });
s.addText("网络中的全部自由基：生成、消耗与对应的表面物种", { placeholder: "title" });
const RADS = [
  ["H•", "H₂ 均裂；OH• + H₂ → H₂O + H•", "加成到 CO₂ / CO / CH₂O；+O₂ → HO₂•", "H*", "通用"],
  ["OH•", "HOCO• → CO + OH•；HO₂• + H•；CH₃OO• 分解", "+H₂ → H₂O + H•；H 抽提 RH → R• + H₂O", "OH*", "通用"],
  ["O•", "O₂ 高温解离", "+CH₄ → •CH₃ + OH•", "O*", "POM / 高温"],
  ["HO₂•", "H• + O₂", "+H• → 2OH•", "—", "POM"],
  ["CO₂•⁻", "CO₂ + e⁻", "+H⁺ + e⁻ → HCOOH（或 → CO + OH⁻）", "CO₂(δ−)*", "电化学"],
  ["HOCO•", "CO₂ + H•", "→ CO + OH•；+H• → HCOOH；⇌ HCOO•", "COOH*", "RWGS"],
  ["HCOO•", "CO₂ + H•", "+H• → HCOOH；⇌ HOCO•", "HCOO*", "甲酸盐路线"],
  ["HOC•HOH / H₂C(OH)O•", "HCOOH + H•", "+H• → CH₂(OH)₂ → CH₂O + H₂O", "H₂COOH* / H₂COO*", "甲酸盐路线"],
  ["HCO•", "CO + H•", "+H• → CH₂O；→ CO + H•", "HCO*", "CO 加氢"],
  ["CH₃O•", "CH₂O + H•；CH₃ONO → CH₃O• + NO•", "+H• → CH₃OH；+CO（Pd）→ COOCH₃", "CH₃O*", "甲醇 / 羰基化"],
  ["•CH₂OH", "CH₂O + H•；CH₃OH 被 H 抽提", "+H• → CH₃OH；⇌ CH₃O•", "CH₂OH*", "甲醇"],
  ["•CH₃", "CH₄ 均裂或 H 抽提；CH₃OH → •CH₃ + OH•", "2•CH₃ → C₂H₆；+O₂ → CH₃OO•；+H• → CH₄", "CH₃*", "OCM / POM / 甲烷化"],
  ["CH₃OO•", "•CH₃ + O₂", "→ CH₂O + OH•", "—", "POM"],
  ["C₂H₅• / CₙH₂ₙ₊₁•", "烷烃 + H• / OH•（H 抽提）；烯烃 + H•", "β-断裂、歧化、复合", "CₙH₂ₙ₊₁*", "OCM / 裂解"],
  ["NO•", "CH₃ONO 在 Pd 上分解", "+CH₃O• + ½O₂ 再生 CH₃ONO", "—", "羰基化"],
];
const th = ["自由基", "生成", "消耗", "对应表面物种", "场景"].map((t) => ({ text: t, options: { bold: true, color: C.background1, fill: { color: C.text2 } } }));
s.addTable([th, ...RADS.map((r, i) => r.map((t, j) => ({ text: t, options: { fill: { color: i % 2 ? C.background1 : C.background2 },
  bold: j === 0, color: j === 0 ? RAD : C.text1 } })))],
  { x: 0.6, y: 1.5, w: 18.8, colW: [2.9, 5.6, 5.6, 2.6, 2.1], rowH: 0.58, fontSize: 13, border: { type: "none" },
    margin: [0, 0.1, 0, 0.1], valign: "middle", objectName: "radical table" });
s.addText("显式枚举（本仓库 data/network_full_C1.csv）：C1 子网络含 10 个自由基、141 个反应；C≤12 全空间约 3×10¹¹ 个自由基，按产物骨架剪枝后为 10²–10⁶ 个。", {
  x: 0.6, y: 11.35, w: 18.8, h: 0.4, fontSize: 12, color: C.text2, margin: 0, isTextBox: true });
s.addNotes("这一页列出总图中出现的全部自由基，以及每个自由基的主要生成与消耗反应和它在催化剂表面上的对应吸附态。");

// ================================================================= slide 3: shared intermediates
s = pres.addSlide({ masterName: "TABLE" });
s.addText("共享中间体：同一个物种贯穿多类反应", { placeholder: "title" });
const ROLE = {
  sCO: "重整终点、费托 / 甲烷化起点、OEO 乙烯酮来源、RWGS 产物", sH: "所有加氢、脱氢步骤的 H 源", sCH3: "甲烷化末步、C–H 活化、CO 插入起点",
  sCH2: "费托链增长单体，加氢走向甲烷", sCO2: "甲酸盐与羧基两条路线的分叉点", mMeOH: "氧化物产物 → 分子筛原料 → 羰基化原料",
  mC24: "MTO / OEO / FTO 汇合点，烯烃循环中心", iCARB: "裂解、异构、低聚的共同中间体", sOH: "SMR 起点、RWGS 与甲酸盐路线的 O 移除",
  sHCOO: "甲醇前驱体；与电化学 HCOOH 相连", sCnH: "β-H 消除、终止、CO 插入的分叉点", iPOOL: "芳烃与轻烯烃的共同来源",
  hubCO: "合成气枢纽", sCH3O: "甲醇前驱体", mHCOOH: "电化学产物，与表面甲酸盐互通", sHCO: "CO 加氢第一步",
  rCH2O: "自由基链中 HCO• → CH₃O•/•CH₂OH 的必经分子", rCH3: "OCM 偶联、POM 氧化、甲烷化的气相对应物",
};
const SHARED = [
  ["hubCO", "合成气枢纽：重整产物，费托、甲醇、OEO、羰基化原料"], ["sCO", "重整终点、费托 / 甲烷化起点、OEO 乙烯酮来源、RWGS 产物"],
  ["sH", "所有加氢、脱氢步骤的 H 源"], ["sCH3", "甲烷化末步、C–H 活化、CO 插入起点"], ["sC", "重整积碳、甲烷化与费托链的起点"],
  ["sCH", "CHₓ 加氢–脱氢链"], ["mMeOH", "氧化物产物 → 分子筛原料 → 羰基化原料"], ["mC24", "MTO / OEO / FTO 汇合点，烯烃循环中心"],
  ["mEtOH", "CO 插入产物，脱水成乙烯"], ["fCH4", "重整原料、甲烷化产物、自由基链起点"], ["sCO2", "甲酸盐与羧基两条路线的分叉点；DRM 原料"],
  ["sOH", "SMR 起点；RWGS 与甲酸盐路线移除 O"], ["sCH2", "费托链增长单体，加氢走向甲烷"], ["iCARB", "裂解、异构、低聚的共同中间体"],
];
const SUB = { 0: "₀", 1: "₁", 2: "₂", 3: "₃", 4: "₄", 5: "₅", 6: "₆", 7: "₇", 8: "₈", 9: "₉", n: "ₙ", "+": "₊", x: "ₓ" };
const SUP = { "+": "⁺", "−": "⁻" };
const label2 = (id) => N[id].label
  .replace(/_\{([^}]*)\}|_(\d+|[a-z])/g, (m, a, b) => (a ?? b).split("").map((c) => SUB[c] ?? c).join(""))
  .replace(/\^\{([^}]*)\}/g, (m, a) => a.split("").map((c) => SUP[c] ?? c).join(""));
const rows3 = SHARED.map(([id, role], i) => {
  const cats = CAT_ORDER.filter((k) => nodeCats[id].has(k));
  const fill = { color: i % 2 ? C.background1 : C.background2 };
  return [{ text: label2(id), options: { bold: true, fill } },
          { text: cats.map((k) => ({ text: "● ", options: { color: CAT[k].c, fontSize: 16 } })), options: { fill } },
          { text: cats.map((k) => CAT[k].n.split(" ")[0]).join(" "), options: { fill, color: C.text2 } },
          { text: role, options: { fill } }];
});
const th3 = ["中间体", "经过的反应类别", "", "作用"].map((t) => ({ text: t, options: { bold: true, color: C.background1, fill: { color: C.text2 } } }));
s.addTable([th3, ...rows3], { x: 0.6, y: 1.5, w: 13.4, colW: [3.0, 2.2, 2.4, 5.8], rowH: 0.6, fontSize: 13, border: { type: "none" },
  margin: [0, 0.1, 0, 0.1], valign: "middle", objectName: "shared intermediates table" });
// counts card
const cnt = {}; Object.values(N).forEach((nd) => { cnt[nd.type] = (cnt[nd.type] || 0) + 1; });
const ecnt = {}; E.forEach((ed) => { ecnt[ed.cat] = (ecnt[ed.cat] || 0) + 1; });
s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 14.5, y: 1.5, w: 4.9, h: 6.3, rectRadius: 0.12, fill: { color: C.background2 },
  line: { color: C.background2 }, objectName: "count card" });
const nTot = Object.values(N).filter((x) => x.type !== "txt").length;
s.addText([
  { text: `${nTot}`, options: { fontSize: 44, bold: true, color: C.text1, breakLine: true } },
  { text: `个节点：分子 ${cnt.mol} · 表面物种 ${cnt.surf} · 自由基 ${cnt.rad} · 离子 ${cnt.ion} · 枢纽 ${cnt.hub}`, options: { fontSize: 14, color: C.text1, breakLine: true } },
  { text: " ", options: { fontSize: 10, breakLine: true } },
  { text: `${E.length}`, options: { fontSize: 44, bold: true, color: C.text1, breakLine: true } },
  { text: "条反应边：", options: { fontSize: 14, color: C.text1, breakLine: true } },
  ...[...CAT_ORDER, "SH"].map((k, i, arr) => ({ text: `● ${CAT[k].n}  ${ecnt[k] || 0}`, options: { fontSize: 13, color: CAT[k].c, breakLine: i < arr.length - 1 } })),
], { x: 14.8, y: 1.75, w: 4.4, h: 5.9, valign: "top", margin: 0, isTextBox: true, objectName: "count text" });
s.addNotes("每个中间体经过的反应类别按化学含义标注（与总图节点圆点一致）；右侧为总图节点与反应边的统计。");

(async () => { await pres.writeFile({ fileName: OUT }); await applyTheme(OUT, THEME); console.log("wrote", OUT, Object.keys(N).length, "nodes", E.length, "edges"); })();
