// 单页版: 把 C1–C8 全部 578 步 + 258 物种排进一张 48×27 英寸 (16:9) 的幻灯片
// 用法: node build_onepage.js <pptx skill 目录>
const fs = require("fs");
const path = require("path");
const pptxgen = require("pptxgenjs");
const { applyTheme } = require(path.join(process.argv[2], "scripts/apply_theme.js"));
const ROOT = path.resolve(__dirname, "..");
const OUT = path.join(ROOT, "FT_network_C1-C8_onepage.pptx");

function readCsv(file) {
  const txt = fs.readFileSync(file, "utf8").replace(/^﻿/, "");
  const rows = txt.split(/\r?\n/).filter(Boolean).map(line => {
    const cells = []; let cur = "", q = false;
    for (let i = 0; i < line.length; i++) {
      const ch = line[i];
      if (q) { if (ch === '"' && line[i + 1] === '"') { cur += '"'; i++; } else if (ch === '"') q = false; else cur += ch; }
      else if (ch === '"') q = true; else if (ch === ",") { cells.push(cur); cur = ""; } else cur += ch;
    }
    cells.push(cur); return cells;
  });
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
const W = 48, H = 27;
pres.defineLayout({ name: "POSTER", width: W, height: H });
pres.layout = "POSTER";
pres.theme = { headFontFace: "Arial", bodyFontFace: "Arial" };
pres.title = "费托合成直链烃完整基元反应网络 (C1–C8) 单页版";
const C = pres.SchemeColor;
const s = pres.addSlide();
s.background = { color: C.background1 };

s.addText("费托合成生成直链烃的完整基元反应网络（C1–C8）",
  { x: 0.7, y: 0.45, w: 30, h: 1.0, fontSize: 54, bold: true, color: C.text2, margin: 0, valign: "middle", isTextBox: true });
s.addText("578 个可逆基元步骤 · 258 个物种 · 含全部自由基型表面中间体 · * 为表面空位 · 橙色粗体为链终止与产物脱附",
  { x: 0.7, y: 1.45, w: 34, h: 0.5, fontSize: 20, color: C.accent6, margin: 0, valign: "middle", isTextBox: true });
s.addText([
  { text: "正构烷烃  n CO + (2n+1) H2 → CnH2n+2 + n H2O", options: { breakLine: true } },
  { text: "直链 α-烯烃  n CO + 2n H2 → CnH2n + n H2O", options: { breakLine: true } },
  { text: "水煤气变换  CO + H2O ⇌ CO2 + H2" },
], { x: 35.3, y: 0.45, w: 12, h: 1.5, fontSize: 18, color: C.text1, align: "right", valign: "middle", margin: 0, isTextBox: true });

const NCOL = 11, GAP = 0.1, X0 = 0.7, CW = (W - 2 * X0 - (NCOL - 1) * GAP) / NCOL, TOP = 2.25, RH = 0.262, FS = 9.5;
const colX = i => X0 + i * (CW + GAP);
const M = [0, 0.05, 0, 0.05];
const top = (t, col) => [{ text: t, options: { colspan: 2, bold: true, fontSize: 14, color: C.background1, fill: { color: col } } }];
const head = t => [{ text: t, options: { colspan: 2, bold: true, color: C.text2, fill: { color: C.background2 } } }];

function rxnColumn(i, title, rows) {
  const body = [top(`${title} · ${rows.length} 步`, C.text2)];
  let cat = null;
  for (const r of rows) {
    if (r.category !== cat) {
      cat = r.category;
      body.push(head(`${cat.replace(/^\S+\s/, "")}（${rows.filter(q => q.category === cat).length}）`));
    }
    const key = /链终止|脱附|甲烷化/.test(r.note);
    body.push([{ text: r.id, options: { color: C.accent6 } }, { text: r.equation, options: key ? { bold: true, color: C.accent2 } : {} }]);
  }
  s.addTable(body, { x: colX(i), y: TOP, w: CW, colW: [0.52, CW - 0.52], rowH: RH, fontSize: FS, color: C.text1, valign: "middle", margin: M, border: { type: "none" }, objectName: `rxn${i}` });
  return body.length;
}

const order = ["1", "2a", "2b", "2c", "3a", "3b", "4"];
const c1 = RX.filter(r => +r.n <= 1).sort((a, b) => order.indexOf(a.category.split(" ")[0]) - order.indexOf(b.category.split(" ")[0]) || a.id.localeCompare(b.id));
const n1 = rxnColumn(0, "通用步骤与 C1", c1);
for (let n = 2; n <= 8; n++) rxnColumn(n - 1, `C${n - 1} → C${n}`, RX.filter(r => +r.n === n));

// 第 1 列下方: 读图说明
const ny = TOP + n1 * RH + 0.5;
s.addShape(pres.ShapeType.roundRect, { x: colX(0), y: ny, w: CW, h: 11.4, rectRadius: 0.03, fill: { color: C.background2 }, line: { type: "none" }, objectName: "guide" });
const g = (t, o = {}) => ({ text: t, options: { breakLine: true, ...o } });
const gh = t => g(t, { bold: true, fontSize: 16, color: C.text2 });
s.addText([
  gh("怎么读这张表"),
  g("第 1 列是与链长无关的步骤和 C1 化学；第 2–8 列每列是增加一个碳的全部步骤；最右 3 列是物种表。"),
  g(" "),
  gh("每增加一个碳的 77 步"),
  g("C–C 偶联（碳化物机理）12 步"),
  g("CO 插入 3 步 + 中间体加氢 21 步 + C–O 断裂 18 步"),
  g("烯醇缩合 1 步；烯基机理 3 步"),
  g("中间体加氢/脱氢 17 步；脱附/再吸附 2 步"),
  g("C1 → C2 因 CHx 两两组合为 81 步。"),
  g(" "),
  gh("自由价"),
  g("片段脱离表面后未饱和的价键数：1 = 单自由基，2 = 卡宾/双自由基，3 = 卡拜型，0 = 闭壳层分子。"),
  g(" "),
  gh("范围"),
  g("只含直链骨架；产物为正构烷烃和直链 α-烯烃。"),
  { text: "不含支链化、内烯烃、含氧产物脱附、气相自由基链反应，也不含能垒与速率常数。" },
], { x: colX(0) + 0.25, y: ny + 0.25, w: CW - 0.5, h: 10.9, fontSize: 13, color: C.text1, valign: "top", margin: 0, paraSpaceAfter: 5, isTextBox: true });

// 物种表: 3 列
const items = [];
let cur = null;
for (const r of SP) {
  if (r.carbon_number !== cur) {
    cur = r.carbon_number;
    const k = SP.filter(q => q.carbon_number === cur).length;
    items.push({ head: `${cur === "0" ? "通用小分子" : "C" + cur} 物种（${k}）` });
  }
  items.push({ r });
}
const per = Math.ceil(items.length / 3);
for (let c = 0; c < 3; c++) {
  const chunk = items.slice(c * per, (c + 1) * per);
  const hd = t => ({ text: t, options: { bold: true, fontSize: 14, color: C.background1, fill: { color: C.accent3 } } });
  const body = [[hd(c ? "物种（续）" : "物种 · 258 个"), hd("自由价"), hd("游离态对应物")]];
  for (const it of chunk) {
    if (it.head) { body.push([{ text: it.head, options: { colspan: 3, bold: true, color: C.text2, fill: { color: C.background2 } } }]); continue; }
    const r = it.r, gasP = r.phase === "gas", rad = !gasP && +r.free_valence > 0;
    let d = gasP ? r.type.replace("产物:", "产物 ") : rad ? r.radical_form.replace("含氧端基: ", "").replace(" (单自由基)", "").replace(" / 1,2-双自由基", "") : (r.type || "闭壳层分子");
    body.push([
      { text: r.species, options: { bold: true, color: gasP ? C.accent2 : C.text1 } },
      { text: gasP ? "–" : r.free_valence, options: { align: "center", bold: rad, color: rad ? C.accent3 : C.accent6 } },
      { text: d },
    ]);
  }
  s.addTable(body, { x: colX(8 + c), y: TOP, w: CW, colW: [1.82, 0.55, CW - 2.37], rowH: RH, fontSize: 10, color: C.text1, valign: "middle", margin: M, border: { type: "none" }, objectName: `sp${c}` });
}

(async () => { await pres.writeFile({ fileName: OUT }); await applyTheme(OUT, THEME); console.log(OUT); })();
