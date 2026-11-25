# CPU C11 优化版

以 AutoDock Vina 上游提交 `3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645` 为基线。正式验收采用 `profiled` 构建，包含 C11 原子对、AD4 与 BFGS 内核，以及 GCC 根据运行采样优化生成代码所用的 21 个 PGO 文件。[上游版本](https://github.com/forlilab/AutoDock-Vina/tree/3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645) · [冻结源码与采样清单](source.json)

## 构建

使用 Linux x86-64、支持 AVX-512F/DQ/VL 的 Ice Lake Server 兼容 CPU、GCC/G++ 11.4.0、GNU Make、ripgrep 和 Boost 1.74.0 静态库。Boost 需要 program_options、filesystem、thread、serialization、system、chrono 和 atomic。

从仓库根目录执行。

```bash
bash cpu/build.sh
```

脚本在 `cpu/rebuilds/` 新建构建目录，打印可执行文件路径及 SHA-256。Boost 默认从 `/usr` 查找，自定义安装可以传入 `BOOST_ROOT` 和 `BOOST_LIBDIR`。正式验收的可执行文件 SHA-256 为 `9468577089e8e219a288d15f09ef7c1547783beeea99d422e85bf3c19c385b24`。与该哈希一致需要匹配冻结环境的编译器、静态库及链接器。

清单记录正式构建文件的哈希。公开 Makefile 将 Boost 默认前缀设为 `/usr`，并沿用冻结构建的编译选项和 PGO 布局。

## 公开样例

将 `vina` 设为上一步打印的可执行文件路径，从仓库根目录运行 1iep。

```bash
mkdir -p .build
example="$PWD/benchmarks/examples/cpu-1iep"
"$vina" --receptor "$example/1iep_receptor.pdbqt" \
    --ligand "$example/1iep_ligand.pdbqt" \
    --config "$example/1iep_receptor.box.txt" \
    --cpu 1 --exhaustiveness 1 --seed 983 --out .build/cpu-1iep.pdbqt
```

正式批次固定四个物理 CPU 核，使用 19 个完整对接配置计时，并为另外 6 个仅评分或局部优化配置检查输出。[实验报告](../docs/report.md)
