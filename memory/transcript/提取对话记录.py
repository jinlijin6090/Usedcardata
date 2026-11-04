import json, re, sys
from datetime import datetime
src, dst = sys.argv[1], sys.argv[2]
SR = re.compile(r"<system-reminder>.*?</system-reminder>", re.S)
out = ["# 对话记录（文字版）", "",
       "本文件由会话日志自动提取：保留用户消息和助手的文字回复，工具调用只保留一行摘要，图片和工具输出省略。完整原始日志见同目录 `session_full.jsonl.gz`（Claude Code 会话格式）。",
       "说明：会话很长，中途有几次上下文压缩；最早的部分在原始日志里完整保留。", "", "---", ""]
nu = 0
last = None
def ts(e):
    t = e.get("timestamp")
    try:
        return datetime.fromisoformat(t.replace("Z", "+00:00")).strftime("%Y-%m-%d %H:%M UTC")
    except Exception:
        return ""
for line in open(src, encoding="utf-8"):
    try:
        e = json.loads(line)
    except Exception:
        continue
    if e.get("type") not in ("user", "assistant") or e.get("isMeta") or e.get("isCompactSummary"):
        continue
    msg = e.get("message") or {}
    c = msg.get("content")
    if e["type"] == "user":
        texts = []
        if isinstance(c, str):
            texts = [c]
        elif isinstance(c, list):
            if any(x.get("type") == "tool_result" for x in c if isinstance(x, dict)):
                continue
            texts = [x.get("text", "") if x.get("type") == "text" else "[图片]" for x in c
                     if isinstance(x, dict) and x.get("type") in ("text", "image")]
            texts = [t for t in texts if not t.startswith("[Image: original")]
        t = SR.sub("", "\n".join(texts)).strip()
        if not t or t.startswith("This session is being continued") or t.startswith("<") or t.startswith("[Request interrupted"):
            continue
        nu += 1
        out += [f"### 用户（{ts(e)}）", "", t, ""]
        last = "user"
    else:
        parts = []
        for x in c if isinstance(c, list) else []:
            if x.get("type") == "text" and x.get("text", "").strip():
                parts.append(x["text"].strip())
            elif x.get("type") == "tool_use":
                inp = x.get("input") or {}
                d = inp.get("description") or inp.get("file_path") or inp.get("pattern") or inp.get("query") or ""
                parts.append(f"> 〔工具 {x.get('name')}〕{str(d)[:160]}")
        if parts:
            out += [f"**助手（{ts(e)}）**", "", "\n".join(parts), ""]
open(dst, "w", encoding="utf-8").write("\n".join(out) + "\n")
print("user messages:", nu)
