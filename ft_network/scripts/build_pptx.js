// 由 data/c8/*.csv 生成 C1–C8 费托反应网络 PPTX（全部 578 步 + 258 物种，逐条列出）
// 用法: node build_pptx.js <pptx skill 目录>
const fs = require("fs");
const path = require("path");
const pptxgen = require("pptxgenjs");
const { applyTheme } = require(path.join(process.argv[2], "scripts/apply_theme.js"));

const ROOT = path.resolve(__dirname, "..");
const OUT = path.join(ROOT, "FT_network_C1-C8.pptx");

function readCsv(file) {
  const txt = fs.readFileSync(file, "utf8").replace(/^﻿/, "");
  const rows = [];
  for (const line of txt.split(/\r?\n/).filter(Boolean)) {
    const cells = []; let cur = "", q = false;
    for (let i = 0; i < line.length; i++) {
      const ch = line[i];
      if (q) { if (ch === '"' && line[i + 1] === '"') { cur += '"'; i++; } else if (ch === '"') q = false; else cur += ch; }
      else if (ch === '"') q = true; else if (ch === ",") { cells.push(cur); cur = ""; } else cur += ch;
    }
    cells.push(cur); rows.push(cells);
  }
  const head = rows.shift();
  return rows.map(r => Object.fromEntries(head.map((h, i) => [h, r[i]])));
}
const RX = readCsv(path.join(ROOT, "data/c8/reactions.csv"));
const SP = readCsv(path.join(ROOT, "data/c8/species.csv"));

const THEME = {
  name: "FT Syngas", headFontFace: "Arial", bodyFontFace: "Arial",
  colors: { dk1: "1B2A33", lt1: "FFFFFF", dk2: "0F4C5C", lt2: "EEF3F4", accent1: "0F4C5C", accent2: "D9580D",
            accent3: "6B2D5C", accent4: "2E8B57", accent5: "3B6FB6", accent6: "8A9AA3", hlink: "3B6FB6", folHlink: "6B2D5C" },
};
const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";            // 13.33 x 7.5 in
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
pres.title = "费托合成直链烃完整基元反应网络 (C1–C8)";
const C = pres.SchemeColor;
const W = 13.33;

pres.defineSlideMaster({
  title: "FT_TITLE", background: { color: C.text2 },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.8, y: 2.3, w: 11.7, h: 1.5, fontSize: 40, bold: true, color: C.background1, align: "left", valign: "bottom" }, text: "" } },
    { placeholder: { options: { name: "body", type: "body", x: 0.8, y: 4.0, w: 11.7, h: 1.6, fontSize: 20, color: C.background2, align: "left", valign: "top" }, text: "" } },
  ],
});
pres.defineSlideMaster({
  title: "FT_CONTENT", background: { color: C.background1 },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.5, y: 0.3, w: 12.33, h: 0.8, fontSize: 28, bold: true, color: C.text2, align: "left", valign: "middle" }, text: "" } },
    { text: { text: "费托合成直链烃基元反应网络 · C1–C8", options: { x: 0.5, y: 7.05, w: 8, h: 0.3, fontSize: 10, color: C.accent6 } } },
  ],
  slideNumber: { x: 12.2, y: 7.05, w: 0.63, h: 0.3, fontSize: 10, color: C.accent6, align: "right" },
});

let section = "";
function newSection(t) { section = t; pres.addSection({ title: t }); }
function content(title) {
  const s = pres.addSlide({ masterName: "FT_CONTENT", sectionTitle: section });
  s.addText(title, { placeholder: "title" });
  return s;
}

// ---------------------------------------------------------------- 1 标题
newSection("概览");
{
  const s = pres.addSlide({ masterName: "FT_TITLE", sectionTitle: section });
  s.addText("费托合成生成直链烃的完整基元反应网络", { placeholder: "title" });
  s.addText("C1–C8 · 258 个物种 · 578 个可逆基元步骤 · 含全部自由基型表面中间体", { placeholder: "body" });
}

// ---------------------------------------------------------------- 2 规模
{
  const s = content("C1–C8 网络共 258 个物种、578 个可逆基元步骤");
  const surf = SP.filter(r => r.phase === "surface");
  const open = surf.filter(r => +r.free_valence > 0).length;
  const stats = [["578", "可逆基元步骤", C.accent1], [String(SP.length), "物种（表面 + 气相）", C.accent5],
                 [String(open), "开壳层(自由基型)表面中间体", C.accent3], ["15", "烃产物：C1–C8 烷烃 + C2–C8 α-烯烃", C.accent2]];
  stats.forEach(([n, l, col], i) => {
    const x = 0.5 + i * 3.11;
    s.addShape(pres.ShapeType.roundRect, { x, y: 1.4, w: 2.96, h: 2.0, rectRadius: 0.08, fill: { color: C.background2 }, line: { type: "none" }, objectName: `stat${i}` });
    s.addText(n, { x, y: 1.5, w: 2.96, h: 1.1, fontSize: 60, bold: true, color: col, align: "center", valign: "middle", isTextBox: true });
    s.addText(l, { x: x + 0.1, y: 2.6, w: 2.76, h: 0.7, fontSize: 14, color: C.text1, align: "center", valign: "top", isTextBox: true });
  });
  s.addText("总反应", { x: 0.5, y: 3.75, w: 6, h: 0.4, fontSize: 20, bold: true, color: C.text2, isTextBox: true, margin: 0 });
  s.addText([
    { text: "正构烷烃    n CO + (2n+1) H2 → CnH2n+2 + n H2O", options: { breakLine: true } },
    { text: "直链 α-烯烃  n CO + 2n H2 → CnH2n + n H2O", options: { breakLine: true } },
    { text: "水煤气变换   CO + H2O ⇌ CO2 + H2" },
  ], { x: 0.5, y: 4.2, w: 6.2, h: 1.5, fontSize: 16, color: C.text1, isTextBox: true, margin: 0, paraSpaceAfter: 8, valign: "top" });
  s.addText("每个碳数的规模", { x: 7.1, y: 3.75, w: 5.7, h: 0.4, fontSize: 20, bold: true, color: C.text2, isTextBox: true, margin: 0 });
  const hdr = o => ({ bold: true, color: C.background1, fill: { color: C.text2 }, ...o });
  s.addTable([
    [{ text: "部分", options: hdr({}) }, { text: "物种", options: hdr({ align: "center" }) }, { text: "基元步骤", options: hdr({ align: "center" }) }],
    ["通用小分子 + C1", { text: "24", options: { align: "center" } }, { text: "35", options: { align: "center" } }],
    ["C1 → C2", { text: "36", options: { align: "center" } }, { text: "81", options: { align: "center" } }],
    ["C3 … C8（每个碳数）", { text: "33", options: { align: "center" } }, { text: "77", options: { align: "center" } }],
    [{ text: "合计", options: { bold: true } }, { text: "258", options: { bold: true, align: "center" } }, { text: "578", options: { bold: true, align: "center" } }],
  ], { x: 7.1, y: 4.2, w: 5.73, colW: [3.13, 1.2, 1.4], rowH: 0.4, fontSize: 14, color: C.text1, valign: "middle",
       border: { type: "solid", pt: 0.5, color: THEME.colors.lt2 }, objectName: "sizeTable" });
  s.addText("每一步均通过 C/H/O 原子守恒与表面位点守恒校验；* 表示表面空位，所有步骤可逆。",
    { x: 0.5, y: 6.4, w: 12.33, h: 0.4, fontSize: 14, italic: true, color: C.accent6, isTextBox: true, margin: 0 });
}

// ---------------------------------------------------------------- 3 示意图
{
  const s = content("网络骨架：CO 活化生成 CHx 单体，再逐碳增长与终止");
  const iw = 12.33, ih = iw * 1480 / 3100;
  s.addImage({ path: path.join(ROOT, "figures/ft_network_scheme.png"), x: 0.5, y: 1.15, w: iw, h: ih, altText: "费托反应网络示意图：C1 网络与 Cn 中间体格子", objectName: "scheme" });
}

// ---------------------------------------------------------------- 4 格子
{
  const s = content("每个碳数有 12 个烃类中间体，按 α、β 碳的含氢数排成格子");
  const cells = {
    "0,0": ["R–C≡C*", "炔基 · 自由价 5"], "0,1": ["R–C≡CH*", "1-炔(吸附) · 4"], "0,2": ["R–C(*)=CH2", "1-烯-2-基 · 3"], "0,3": ["R–C(*)–CH3", "2-亚烷基(卡宾) · 2"],
    "1,0": ["R–CH=C*", "亚乙烯基型 · 4"], "1,1": ["R–CH=CH*", "烯基 · 3"], "1,2": ["R–CH=CH2*", "吸附 α-烯烃 · 2"], "1,3": ["R–CH(*)–CH3", "仲烷基自由基 · 1"],
    "2,0": ["R–CH2–C*", "次烷基(卡拜) · 3"], "2,1": ["R–CH2–CH*", "亚烷基(卡宾) · 2"], "2,2": ["R–CH2–CH2*", "烷基自由基 · 1"], "2,3": ["R–CH2–CH3*", "正构烷烃 · 0"],
  };
  const x0 = 2.0, y0 = 1.95, cw = 2.2, ch = 1.15, gx = 0.2, gy = 0.2;
  for (let b = 0; b < 4; b++)
    s.addText(`α-碳 H 数 b = ${b}`, { x: x0 + b * (cw + gx), y: 1.45, w: cw, h: 0.4, fontSize: 14, bold: true, color: C.text2, align: "center", isTextBox: true });
  for (let a = 0; a < 3; a++) {
    s.addText(`β-碳 H 数\na = ${a}`, { x: 0.5, y: y0 + a * (ch + gy), w: 1.4, h: ch, fontSize: 14, bold: true, color: C.text2, align: "center", valign: "middle", isTextBox: true });
    for (let b = 0; b < 4; b++) {
      const [f, k] = cells[`${a},${b}`];
      const prod = (a === 2 && b === 3) || (a === 1 && b === 2), grow = a === 2 && b < 3;
      s.addShape(pres.ShapeType.roundRect, { x: x0 + b * (cw + gx), y: y0 + a * (ch + gy), w: cw, h: ch, rectRadius: 0.08,
        fill: { color: prod ? C.accent2 : grow ? C.accent4 : C.background2 }, line: { type: "none" }, objectName: `cell${a}${b}` });
      const tc = prod || grow ? C.background1 : C.text1;
      s.addText([{ text: f, options: { bold: true, fontSize: 16, breakLine: true } }, { text: k, options: { fontSize: 12 } }],
        { x: x0 + b * (cw + gx), y: y0 + a * (ch + gy), w: cw, h: ch, color: tc, align: "center", valign: "middle", isTextBox: true });
    }
  }
  const lx = 11.75;
  [["链增长片段", C.accent4, "再接 CHx 或 CO"], ["脱附得产物", C.accent2, "烷烃 / α-烯烃"], ["其余中间体", C.background2, "加氢/脱氢互变"]].forEach(([t, col, d], i) => {
    s.addShape(pres.ShapeType.roundRect, { x: lx, y: 2.0 + i * 1.35, w: 0.3, h: 0.3, rectRadius: 0.05, fill: { color: col }, line: { type: "none" }, objectName: `leg${i}` });
    s.addText([{ text: t, options: { bold: true, breakLine: true } }, { text: d }], { x: lx - 0.05, y: 2.35 + i * 1.35, w: 1.5, h: 0.7, fontSize: 12, color: C.text1, isTextBox: true, margin: 0, valign: "top" });
  });
  s.addText([
    { text: "相邻格点之间是 ±H* 的基元步骤（17 步/碳数）。R = CH3(CH2)n-3–；C2 因两端对称只有 10 个状态。", options: { breakLine: true } },
    { text: "“自由价”= 片段脱离表面后未饱和的价键数：1 为单自由基，2 为卡宾/双自由基，3 为卡拜型。" },
  ], { x: 0.5, y: 6.1, w: 12.33, h: 0.8, fontSize: 14, color: C.text1, isTextBox: true, margin: 0, paraSpaceAfter: 4, valign: "top" });
}

// ---------------------------------------------------------------- 5 机理
{
  const s = content("每增加一个碳重复 77 步：四种链增长机理加两条终止路径");
  const cards = [
    ["碳化物机理 · 12 步", "R'–CH_a* + CH_b* ⇌ R'–CH_a–CH_b* + *", "a = 0,1,2；b = 0…3。含烷基机理（CH2 插入）、亚烷基与次烷基路径", C.accent4],
    ["CO 插入机理 · 42 步", "R'–CH_a* + CO* ⇌ R'–CH_a–CO* + *", "插入 3 步，含氧中间体加氢 21 步，C–O 断裂 18 步", C.accent3],
    ["烯醇缩合机理 · 1 步", "R–C(OH)* + HCOH* ⇌ R–C–COH* + H2O*", "羟基卡宾脱水缩合成 C–C 键", C.accent5],
    ["烯基机理 · 3 步", "R'CH=CH* + CH2* ⇌ η3-烯丙基 ⇌ R'CH2CH=CH*", "烯丙基经 1,3-H 迁移成高一级烯基，或加氢成 α-烯烃", C.accent1],
    ["加氢/脱氢 · 17 步", "R–CH_a–CH_b* + H* ⇌ 相邻格点 + *", "含烷基加氢 → 正构烷烃，烷基 β-H 消除 → α-烯烃", C.accent2],
    ["脱附/再吸附 · 2 步", "RCH2CH3* ⇌ CnH2n+2(g)；RCH=CH2* ⇌ CnH2n(g)", "逆向即烯烃再吸附，引出二次加氢与再引发", C.accent2],
  ];
  cards.forEach(([t, eq, d, col], i) => {
    const x = 0.5 + (i % 3) * 4.18, y = 1.4 + Math.floor(i / 3) * 2.75;
    s.addShape(pres.ShapeType.roundRect, { x, y, w: 3.97, h: 2.5, rectRadius: 0.06, fill: { color: C.background2 }, line: { type: "none" }, objectName: `mech${i}` });
    s.addShape(pres.ShapeType.ellipse, { x: x + 0.25, y: y + 0.25, w: 0.5, h: 0.5, fill: { color: col }, line: { type: "none" }, objectName: `mechDot${i}` });
    s.addText(String(i + 1), { x: x + 0.25, y: y + 0.25, w: 0.5, h: 0.5, fontSize: 16, bold: true, color: C.background1, align: "center", valign: "middle", isTextBox: true, margin: 0 });
    s.addText(t, { x: x + 0.9, y: y + 0.25, w: 2.9, h: 0.5, fontSize: 18, bold: true, color: C.text2, valign: "middle", isTextBox: true, margin: 0 });
    s.addText(eq, { x: x + 0.25, y: y + 0.95, w: 3.5, h: 0.65, fontSize: 14, bold: true, color: C.text1, valign: "top", isTextBox: true, margin: 0 });
    s.addText(d, { x: x + 0.25, y: y + 1.65, w: 3.5, h: 0.7, fontSize: 14, color: C.text1, valign: "top", isTextBox: true, margin: 0 });
  });
}

// ---------------------------------------------------------------- 反应表 / 物种表
const COLX = [0.5, 6.84], COLW = 5.99, TOP = 1.25;
const RROW = 0.24, RCAP = 23;       // 反应表: 每列 23 行
const SROW = 0.285, SCAP = 19;      // 物种表: 表头 + 18 行

function reactionSlides(title, cats, rows) {
  // 把 (类别表头 + 反应行) 流式排进两列, 满了换页
  const items = [];
  for (const cat of cats) {
    const rs = rows.filter(r => r.category === cat);
    if (!rs.length) continue;
    items.push({ head: `${cat}（${rs.length} 步）`, cat });
    rs.forEach(r => items.push({ r, cat }));
  }
  const cols = [[]];
  for (let i = 0; i < items.length; i++) {
    let col = cols[cols.length - 1];
    const it = items[i];
    if (col.length >= RCAP || (it.head && col.length >= RCAP - 1)) {
      col = []; cols.push(col);
      if (!it.head) col.push({ head: `${it.cat}（续）` });
    }
    col.push(it);
  }
  for (let p = 0; p * 2 < cols.length; p++) {
    const s = content(title + (p ? "（续）" : ""));
    [cols[2 * p], cols[2 * p + 1]].forEach((col, ci) => {
      if (!col) return;
      const body = col.map(it => {
        if (it.head) return [{ text: it.head.replace(/^\S+\s/, ""), options: { colspan: 2, bold: true, color: C.background1, fill: { color: C.text2 } } }];
        const key = /链终止|脱附|甲烷化/.test(it.r.note);
        const o = key ? { bold: true, color: C.accent2 } : {};
        return [{ text: it.r.id, options: { color: C.accent6, ...(key ? { bold: true } : {}) } }, { text: it.r.equation, options: o }];
      });
      s.addTable(body, { x: COLX[ci], y: TOP, w: COLW, colW: [0.72, COLW - 0.72], rowH: RROW, fontSize: 12, color: C.text1, valign: "middle",
        margin: [0, 0.06, 0, 0.06], border: { type: "none" }, objectName: `rxn${ci}` });
    });
  }
}

function speciesSlides(title, rows) {
  const per = SCAP - 1;
  for (let p = 0; p * 2 * per < rows.length; p++) {
    const s = content(title + (p ? "（续）" : ""));
    for (let ci = 0; ci < 2; ci++) {
      const chunk = rows.slice((2 * p + ci) * per, (2 * p + ci + 1) * per);
      if (!chunk.length) continue;
      const h = t => ({ text: t, options: { bold: true, color: C.background1, fill: { color: C.text2 } } });
      const body = [[h("物种"), h("自由价"), h("类型 / 游离态对应物")]];
      chunk.forEach((r, i) => {
        const gasP = r.phase === "gas", rad = !gasP && +r.free_valence > 0;
        const fill = i % 2 ? { fill: { color: C.background2 } } : {};
        let desc = gasP ? r.type : (r.radical_form.startsWith("含氧端基") || !r.type ? (r.radical_form || r.type) : `${r.type.replace(/ [a-z]+$/, "")} · ${r.radical_form}`);
        if (!gasP && !rad) desc = r.type || "闭壳层分子";
        if (desc.length > 30) desc = r.radical_form;
        body.push([
          { text: r.species, options: { bold: true, color: gasP ? C.accent2 : C.text1, ...fill } },
          { text: gasP ? "–" : r.free_valence, options: { align: "center", color: rad ? C.accent3 : C.accent6, bold: rad, ...fill } },
          { text: desc, options: { ...fill } },
        ]);
      });
      s.addTable(body, { x: COLX[ci], y: TOP, w: COLW, colW: [2.15, 0.6, COLW - 2.75], rowH: SROW, fontSize: 11, color: C.text1, valign: "middle",
        margin: [0, 0.06, 0, 0.06], border: { type: "none" }, objectName: `sp${ci}` });
    }
  }
}

const catsOf = n => [...new Set(RX.filter(r => +r.n === n).map(r => r.category))];

// C1 + 通用
newSection("通用步骤与 C1");
reactionSlides("通用步骤与 C1 化学：吸附、CO 活化、除氧、CHx 加氢（35 步）",
  [...catsOf(0).filter(c => c.startsWith("1")), ...catsOf(1).filter(c => c.startsWith("2")), ...catsOf(0).filter(c => c.startsWith("3a")),
   ...catsOf(1).filter(c => c.startsWith("3b")), ...catsOf(1).filter(c => c.startsWith("4"))],
  RX.filter(r => +r.n <= 1));
speciesSlides("通用小分子与 C1 物种（24 个）", SP.filter(r => +r.carbon_number <= 1));

for (let n = 2; n <= 8; n++) {
  newSection(`C${n - 1} → C${n}`);
  const rows = RX.filter(r => +r.n === n);
  const A = ["5a", "5b", "5e", "5f", "5g", "5h"], B = ["5c", "5d"];
  const pick = keys => catsOf(n).filter(c => keys.some(k => c.startsWith(k))).sort((x, y) => keys.indexOf(x.slice(0, 2)) - keys.indexOf(y.slice(0, 2)));
  const nA = rows.filter(r => A.some(k => r.category.startsWith(k))).length;
  reactionSlides(`C${n - 1} → C${n}：链增长、加氢/脱氢与链终止（${nA} 步）`, pick(A), rows);
  reactionSlides(`C${n - 1} → C${n}：CO 插入中间体的加氢与 C–O 断裂（${rows.length - nA} 步）`, pick(B), rows);
  const sp = SP.filter(r => +r.carbon_number === n);
  speciesSlides(`C${n} 物种（${sp.length} 个）`, sp);
}

// ---------------------------------------------------------------- 范围
newSection("范围");
{
  const s = content("网络的范围：只含直链骨架，产物为正构烷烃与 α-烯烃");
  const items = [
    ["包含", C.accent4, ["CO 直接解离与氢助解离的全部 C1 路径", "碳化物、CO 插入、烯醇缩合、烯基四种链增长机理", "烷基加氢与 β-H 消除两条终止路径", "除氧生成 H2O / CO2 与水煤气变换", "可逆性自带的二次反应：烯烃再吸附、二次加氢、氢解"]],
    ["不包含", C.accent2, ["支链化（CHx 加到链内碳）", "双键位移生成内烯烃", "醇、醛、酸等含氧产物的脱附", "气相自由基链反应", "能垒、速率常数及各催化剂上的主导路径判断"]],
  ];
  items.forEach(([t, col, list], i) => {
    const x = 0.5 + i * 6.34;
    s.addShape(pres.ShapeType.roundRect, { x, y: 1.4, w: 5.99, h: 4.6, rectRadius: 0.05, fill: { color: C.background2 }, line: { type: "none" }, objectName: `scope${i}` });
    s.addShape(pres.ShapeType.ellipse, { x: x + 0.3, y: 1.7, w: 0.4, h: 0.4, fill: { color: col }, line: { type: "none" }, objectName: `scopeDot${i}` });
    s.addText(t, { x: x + 0.85, y: 1.65, w: 4.8, h: 0.5, fontSize: 22, bold: true, color: C.text2, valign: "middle", isTextBox: true, margin: 0 });
    s.addText(list.map((l, k) => ({ text: l, options: { bullet: true, breakLine: k < list.length - 1 } })),
      { x: x + 0.3, y: 2.4, w: 5.4, h: 3.4, fontSize: 16, color: C.text1, valign: "top", isTextBox: true, margin: 0, paraSpaceAfter: 10 });
  });
  s.addText("同一套规则可延伸到任意碳数：python scripts/build_network.py --nmax N",
    { x: 0.5, y: 6.3, w: 12.33, h: 0.4, fontSize: 14, italic: true, color: C.accent6, isTextBox: true, margin: 0 });
}

(async () => {
  await pres.writeFile({ fileName: OUT });
  await applyTheme(OUT, THEME);
  console.log(OUT);
})();
