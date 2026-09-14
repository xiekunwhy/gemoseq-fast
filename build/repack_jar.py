#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""repack_jar.py -- 重新打包 GeMoSeq-1.2.3-fast.jar

将官方 GeMoSeq fat jar 中的 htsjdk 替换为新版，并用本仓库 src/gemoseq 编译出的
class 覆盖同名类。

用法:
    python repack_jar.py \
        --base GeMoSeq-1.2.3.jar \
        --htsjdk htsjdk-2.24.1.jar \
        --classes ../out-fast \
        --out GeMoSeq-1.2.3-fast.jar

其中 --classes 指向先执行过如下编译的输出目录（Windows 下 classpath 用分号）:
    javac -encoding UTF-8 -cp GeMoSeq-1.2.3.jar;htsjdk-2.24.1.jar -d out-fast src/gemoseq/*.java
"""
import argparse
import os
import zipfile


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, help="官方 GeMoSeq-1.2.3.jar 路径")
    ap.add_argument("--htsjdk", required=True, help="新版 htsjdk jar 路径")
    ap.add_argument("--classes", required=True, help="编译输出目录（含 projects/gemoseq/*.class）")
    ap.add_argument("--out", required=True, help="输出 jar 路径")
    a = ap.parse_args()

    newfiles = {}
    for root, _, files in os.walk(a.classes):
        for f in files:
            if f.endswith(".class"):
                full = os.path.join(root, f)
                rel = os.path.relpath(full, a.classes).replace(os.sep, "/")
                newfiles[rel] = full
    if not newfiles:
        raise SystemExit("no .class files found under %s" % a.classes)

    zin = zipfile.ZipFile(a.base)
    zh = zipfile.ZipFile(a.htsjdk)
    seen = set()
    with zipfile.ZipFile(a.out, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            if item.filename in seen or item.filename.startswith("htsjdk/"):
                continue  # 去重 + 剔除旧版 htsjdk
            seen.add(item.filename)
            if item.filename in newfiles:
                with open(newfiles.pop(item.filename), "rb") as fh:
                    zout.writestr(item.filename, fh.read())
            else:
                zout.writestr(item, zin.read(item.filename))
        for item in zh.infolist():
            if item.filename in seen or item.filename.endswith("/"):
                continue
            seen.add(item.filename)
            zout.writestr(item, zh.read(item.filename))
        for rel, full in newfiles.items():
            zout.write(full, rel)
    print("built %s (%d bytes, %d patched/new classes)" % (a.out, os.path.getsize(a.out), len(newfiles)))


if __name__ == "__main__":
    main()
