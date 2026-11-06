#!/usr/bin/env bash
# Offline install of the bundled environment (no internet, no root, no conda needed).
#   bash env/install_env.sh [path/to/nicoke_env.tar.gz]
# Default archive location: env/nicoke_env.tar.gz
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ARCHIVE="${1:-$HERE/nicoke_env.tar.gz}"
TARGET="$HERE/nicoke_env"

echo "== system check"
GLIBC=$(ldd --version 2>/dev/null | head -1 | grep -oE '[0-9]+\.[0-9]+$' || echo unknown)
echo "   glibc: $GLIBC (bundle needs >= 2.27)   arch: $(uname -m) (needs x86_64)"
if [ "$(uname -m)" != "x86_64" ]; then echo "ERROR: x86_64 only"; exit 1; fi
if [ "$GLIBC" != unknown ] && [ "$(printf '%s\n2.27\n' "$GLIBC" | sort -V | head -1)" != "2.27" ]; then
    echo "ERROR: glibc $GLIBC is older than 2.27 (e.g. CentOS 7). See README section 'old systems'."; exit 1
fi
if [ ! -f "$ARCHIVE" ] && ls "$HERE"/parts/nicoke_env.tar.gz.part* >/dev/null 2>&1; then
    echo "== joining split archive from env/parts/ (git checkout)"
    cat "$HERE"/parts/nicoke_env.tar.gz.part* > "$ARCHIVE"
    if [ -f "$HERE/parts/SHA256SUM" ]; then
        (cd "$HERE" && sha256sum -c parts/SHA256SUM) || { echo "ERROR: checksum mismatch - parts incomplete (git lfs/partial clone?)"; rm -f "$ARCHIVE"; exit 1; }
    fi
fi
[ -f "$ARCHIVE" ] || { echo "ERROR: $ARCHIVE not found. Copy nicoke_env.tar.gz into env/ first."; exit 1; }

echo "== unpacking $(du -h "$ARCHIVE" | cut -f1) archive into $TARGET"
rm -rf "$TARGET"
mkdir -p "$TARGET"
tar -xzf "$ARCHIVE" -C "$TARGET"
echo "== fixing prefixes (conda-unpack)"
"$TARGET/bin/python" "$TARGET/bin/conda-unpack" >/dev/null 2>&1 || "$TARGET/bin/conda-unpack"
echo "== quick test"
bash "$HERE/check_env.sh"
