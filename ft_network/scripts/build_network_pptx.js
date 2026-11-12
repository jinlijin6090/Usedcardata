// 网络图版: 一页 36.4×32.6 英寸幻灯片, 放 figures/ft_network_C1-C8_full.png
// 用法: node build_network_pptx.js <pptx skill 目录>   (先运行 plot_full_network.py)
const path = require("path");
const pptxgen = require("pptxgenjs");
const { applyTheme } = require(path.join(process.argv[2], "scripts/apply_theme.js"));
const ROOT = path.resolve(__dirname, "..");
const N = Number(process.argv[3] || 8);              // 最大碳数: 8 或 12
const PW = 36.4, PH = N <= 8 ? 32.6 : 47.0;
const NSP = 24 + 36 + 33 * (N - 2), NRX = 35 + 81 + 77 * (N - 2);
const OUT = path.join(ROOT, `FT_network_C1-C${N}_graph.pptx`);
const THEME = {
  name: "FT Syngas", headFontFace: "Arial", bodyFontFace: "Arial",
  colors: { dk1: "1B2A33", lt1: "FFFFFF", dk2: "0F4C5C", lt2: "EEF3F4", accent1: "0F4C5C", accent2: "D9580D",
            accent3: "6B2D5C", accent4: "2E8B57", accent5: "3B6FB6", accent6: "8A9AA3", hlink: "3B6FB6", folHlink: "6B2D5C" },
};
const pres = new pptxgen();
pres.defineLayout({ name: "POSTER", width: PW, height: PH });
pres.layout = "POSTER";
pres.theme = { headFontFace: "Arial", bodyFontFace: "Arial" };
pres.title = `费托合成直链烃完整反应网络图 (C1–C${N})`;
const C = pres.SchemeColor;
const s = pres.addSlide();
s.background = { color: C.background1 };
s.addImage({ path: path.join(ROOT, `figures/ft_network_C1-C${N}_full.png`), x: 0, y: 0, w: PW, h: PH, objectName: "network",
  altText: `费托合成 C1–C${N} 反应网络图：${NSP} 个物种节点、${NRX} 条基元步骤箭头，按碳数分块、分段阶梯排列` });
s.addText(`费托合成生成直链烃的完整反应网络（C1–C${N}）`,
  { x: 14.0, y: 0.7, w: 21.8, h: 1.1, fontSize: 48, bold: true, color: C.text2, margin: 0, valign: "middle", isTextBox: true });
s.addText(`${NSP} 个物种 · ${NRX} 个可逆基元步骤 · 含全部自由基型表面中间体`,
  { x: 14.0, y: 1.85, w: 21.8, h: 0.6, fontSize: 24, color: C.accent6, margin: 0, valign: "middle", isTextBox: true });
(async () => { await pres.writeFile({ fileName: OUT }); await applyTheme(OUT, THEME); console.log(OUT); })();
