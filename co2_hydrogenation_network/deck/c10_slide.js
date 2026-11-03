// One-slide C1–C10 fully isomer-resolved network: rendered graph image + native labels/legend.
// Inputs: c10_net.jpg, c10_layout.json (from c10_layout.py). Usage: node c10_slide.js out.pptx
const pptxgen = require("pptxgenjs");
const fs = require("fs");
const L = JSON.parse(fs.readFileSync("c10_layout.json", "utf8"));
const out = process.argv[2] || "c10.pptx";
const SW = 44, SH = 24.75, X0 = 2.6, Y0 = 1.55;
const INK = "1B1B1B", INK2 = "55534E";
const BLK = { MET: "D64545", GAS: "8A8A8A", OX: "D99A00", BA: "2A6FD6", BA2: "7FA8E6", AR: "1E9C4B" };
const BLKNAME = {
  MET: "金属/碳化物 + 氧化物：链增长 · CO 插入 · 含氧化合物", GAS: "气相自由基", OX: "氧化物 H 溢流加氢 / 脱氢",
  BA: "分子筛 Brønsted 酸：碳正离子", BA2: "分子筛：氢转移 / 成环", AR: "芳烃池（BA/LA）",
};

const pres = new pptxgen();
pres.defineLayout({ name: "BIG", width: SW, height: SH });
pres.layout = "BIG";
const s = pres.addSlide();
s.background = { color: "FFFFFF" };
const fmt = (n) => n.toLocaleString("en-US");

s.addText([
  { text: "CO₂ 一步热催化加氢 · C1–C10 全同分异构体反应网络", options: { bold: true, fontSize: 30, color: INK } },
  { text: `   每个点 = 一个构造异构体（含表面物种 / 自由基 / 碳正离子），共 ${fmt(L.n_nodes)} 个节点 · ${fmt(L.n_edges)} 条基元反应边`, options: { fontSize: 16, color: INK2 } },
], { x: 0.25, y: 0.18, w: SW - 0.5, h: 0.75, margin: 0, valign: "middle" });

// column headers
for (const [c, [a, b]] of Object.entries(L.colx)) {
  const label = c === "1" ? "C1 表面核心" : `C${c}`;
  const sub = c === "1" ? "" : `  ${fmt(L.cols[c])}`;
  s.addText([
    { text: label, options: { bold: true, fontSize: c === "1" ? 13 : 15, color: INK } },
    { text: sub, options: { fontSize: 11, color: INK2 } },
  ], { x: X0 + a, y: 1.08, w: b - a, h: 0.42, align: "center", valign: "middle", margin: 0 });
}
// row labels
L.rows.forEach(([key, lab, blk, n], r) => {
  const [y0, y1] = L.rowy[String(r)];
  s.addShape(pres.shapes.RECTANGLE, { x: 0.15, y: Y0 + y0 + 0.04, w: 0.09, h: y1 - y0 - 0.08, fill: { color: BLK[blk] }, line: { color: BLK[blk], width: 0 } });
  s.addText([
    { text: lab, options: { bold: true, fontSize: 12.5, color: INK, breakLine: true } },
    { text: `${fmt(n)} 个`, options: { fontSize: 11, color: INK2 } },
  ], { x: 0.32, y: Y0 + y0, w: X0 - 0.4, h: y1 - y0, valign: "middle", margin: 0, fit: "shrink" });
});

s.addImage({ path: "c10_net.jpg", x: X0, y: Y0, w: L.W, h: L.H });

// legend: node classes
const NODES = [["醇", "8E5BB5"], ["醛", "B5487F"], ["酮", "C77A2E"], ["酸", "6A3D9A"], ["R*", "D64545"], ["R•", "6E6E6E"],
  ["烷烃", "3A3A3A"], ["烯烃", "E07B00"], ["R⁺", "1F5FBF"], ["二烯", "C9A100"], ["烯基⁺", "3D8FE0"], ["环烷", "5C6B73"],
  ["环烯", "2F9C95"], ["环状⁺", "6C4FD1"], ["环二烯", "0E7C66"], ["芳烃", "1E8C3A"]];
let x = 0.25, y = 23.72;
s.addText("节点", { x, y, w: 0.6, h: 0.32, fontSize: 12, bold: true, color: INK, margin: 0, valign: "middle" });
x += 0.62;
for (const [t, c] of NODES) {
  s.addShape(pres.shapes.OVAL, { x, y: y + 0.08, w: 0.17, h: 0.17, fill: { color: c }, line: { color: c, width: 0 } });
  s.addText(t, { x: x + 0.22, y, w: 0.85, h: 0.32, fontSize: 11, color: INK, margin: 0, valign: "middle" });
  x += 0.3 + 0.17 * t.length + 0.35;
}
const ef = L.edge_fams, sum = (ks) => ks.reduce((a, k) => a + (ef[k] || 0), 0);
const EDGES = [
  ["质子化/去质子", "5B9BE8", sum(["Z:prot"])],
  ["H 迁移 · 烷基迁移 · PCP · β-断裂/齐聚 · 成环 · 甲基化", "2A6FD6", sum(["Z:hshift", "Z:alkylshift", "Z:pcp", "Z:beta", "Z:cyclise", "Z:meth"])],
  ["氢转移 HT", "13A39A", sum(["Z:HT"])],
  ["加氢/脱氢", "D99A00", sum(["H:hyd"])],
  ["金属位 β-H · CH₂/CO 插入 · 加氢", "D64545", sum(["M:betaH", "M:CH2", "M:CO", "M:hyd", "M:hub"])],
  ["含氧化合物 加氢/脱水", "8E5BB5", sum(["O:dehyd", "O:hyd", "O:ox"])],
  ["自由基 β-断裂 · H 抽提", "777777", sum(["R:beta", "R:abs"])],
  ["芳烃 甲基化 · 异构 · 烷基化", "1E9C4B", sum(["A:meth", "A:iso", "A:FC", "A:hub"])],
];
x = 0.25; y = 24.08;
s.addText("边", { x, y, w: 0.6, h: 0.32, fontSize: 12, bold: true, color: INK, margin: 0, valign: "middle" });
x += 0.62;
for (const [t, c, n] of EDGES) {
  s.addShape(pres.shapes.LINE, { x, y: y + 0.16, w: 0.45, h: 0, line: { color: c, width: 2.5 } });
  const lab = `${t}  ${fmt(n)}`;
  const w = 0.13 * lab.length + 0.4;
  s.addText(lab, { x: x + 0.52, y, w, h: 0.32, fontSize: 11, color: INK, margin: 0, valign: "middle" });
  x += 0.52 + w + 0.25;
}
// block legend (right side of node row)
x = 19.0; y = 23.72;
for (const k of ["MET", "GAS", "OX", "BA", "BA2", "AR"]) {
  s.addShape(pres.shapes.RECTANGLE, { x, y: y + 0.07, w: 0.09, h: 0.2, fill: { color: BLK[k] }, line: { color: BLK[k], width: 0 } });
  const w = 0.135 * BLKNAME[k].length + 0.2;
  s.addText(BLKNAME[k], { x: x + 0.15, y, w, h: 0.32, fontSize: 10, color: INK2, margin: 0, valign: "middle" });
  x += w + 0.35;
  if (x > 42.5) { x = 26.2; y += 0.0; }
}
s.addText("列 = 碳数，行 = 物种类别，格内数字 = 该格异构体数；格内按骨架排序（直链在左）。不含伯碳正离子；双分子负氢转移（~10⁷ 对）未画。每个点下方有 SMILES，放大可读（另附矢量 PDF）。",
  { x: 0.25, y: 24.44, w: SW - 0.5, h: 0.25, fontSize: 10, color: INK2, margin: 0, valign: "middle" });

pres.writeFile({ fileName: out }).then(() => console.log("wrote", out));
