# AccelVina

**在 Intel Xeon CPU、NVIDIA A100 和摩尔线程 S5000 上加速 AutoDock Vina 分子对接，输出与同平台基线一致。**

[实验报告](docs/report.md) · [HTML 版报告](docs/report.html) · [公开样例与复测](benchmarks/) · [上游 AutoDock Vina](https://github.com/forlilab/AutoDock-Vina/tree/3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645) · [上游 Vina-GPU 2.1](https://github.com/DeltaGroupNJUPT/Vina-GPU-2.1/tree/180272b8a5265d6ed9664178345933cebe2cd349)

本仓库收录三条优化路线正式验收时使用的源码、构建脚本、公开复测样例和实验报告。三条路线都不改评分函数和搜索逻辑：相同输入、随机种子和设置下，优化版输出的 PDBQT 与同平台基线一致。

## 效果

| 平台 | 同平台基线 | 端到端耗时 | 耗时减少 | 峰值内存或显存 | 输出校验 |
| --- | --- | ---: | ---: | ---: | --- |
| [Intel Xeon Platinum 8362 · C11](cpu/) | AutoDock Vina 上游 | 3445.941 → 1936.001 秒 | **43.82%** | 828.30 → 418.57 MiB（−49.47%） | 75/75 次通过 |
| [NVIDIA A100 · OpenCL](gpu-opencl/) | Vina-GPU 2.1 上游 | 70.042 → 45.226 秒 | **35.43%** | 19.046 → 18.499 GiB（−2.87%） | 两版合计 150/150 次通过 |
| [摩尔线程 S5000 · MUSA](gpu-musa/) | 同卡 MUSA 移植版 | 1181.236 → 128.220 秒 | **89.15%** | 3542 → 620 MiB（−82.50%） | 两版合计 150/150 次通过 |

- **计时口径**：每个用例运行三次取中位数，再把各用例相加。CPU 覆盖 19 个完整对接配置，两条 GPU 路线各覆盖 25 个最终用例。
- **内存口径**：CPU 为进程峰值驻留内存（RSS）；A100 为公开单例的整卡显存采样峰值；S5000 为正式批次的进程显存采样峰值。
- **输出校验**：比较完整 PDBQT 文件，CPU 另比较规范化日志，GPU 另比较姿态表。参考输出分别是 CPU 上游、A100 上游和 S5000 上的 MUSA 移植版。
- **S5000 的基线**是我们把 Vina-GPU 2.1 移植到 MUSA 后的未优化版本，不是官方发行版。跨平台输出存在数值差异，例如 2bm2 样例的最佳亲和力在 A100 上为 −8.6 kcal/mol，在 S5000 上为 −8.1 kcal/mol。

各版本的对照、失败尝试和测量细节见 [实验报告](docs/report.md)。

## 做了哪些优化

- **CPU**：原子对求值、AD4 网格插值和 BFGS 局部优化改用 C11 + AVX-512 内核，保持原有累加顺序；蒙特卡洛搜索仍是 C++；最后用 GCC 的剖析引导优化（PGO）编译。
- **A100**：两份 OpenCL 内核保持不变，只改主机端：只为实际用到的原子类型建网格、只上传网格头、批处理复用受体解析。公开单例的主机到设备上传量从 1115.17 MiB 降到 13.53 MiB。
- **S5000**：先把 Vina-GPU 2.1 的 OpenCL 内核移植到 MUSA 作为基线，再优化主机数据准备、按实际访问范围映射显存，并在网格内核报警时重放原内核以保留原有输出。

## 仓库结构

| 目录 | 内容 |
| --- | --- |
| [`cpu/`](cpu/) | CPU C11 优化版源码、PGO 剖析文件和构建脚本 |
| [`gpu-opencl/`](gpu-opencl/) | A100 OpenCL 优化版源码与构建文件 |
| [`gpu-musa/`](gpu-musa/) | S5000 MUSA 移植补丁、优化补丁、构建脚本和复测脚本 |
| [`upstream/vina-gpu-2.1/`](upstream/vina-gpu-2.1/) | 两条 GPU 路线共用的 Vina-GPU 2.1 上游源码及逐文件校验和 |
| [`benchmarks/`](benchmarks/) | 公开样例 1iep（CPU）和 2bm2（GPU）及复测清单 |
| [`docs/`](docs/) | 实验报告（Markdown 和 HTML） |
| `build_baseline.sh` · `run_benchmark.py` | A100 基线构建与 GPU 对照复测入口 |

## 快速开始

每个平台的依赖、构建命令和公开样例都写在各自目录的 README 里：[CPU](cpu/README.md) · [A100 OpenCL](gpu-opencl/README.md) · [S5000 MUSA](gpu-musa/README.md)。

以 CPU 为例（需要 Linux x86-64、支持 AVX-512 的 Ice Lake Server 兼容 CPU、GCC 11.4 和 Boost 1.74 静态库）：

```bash
git clone https://github.com/holamochi/AccelVina.git
cd AccelVina
bash cpu/build.sh          # 打印可执行文件路径及 SHA-256

vina=<上一步打印的路径>
example=benchmarks/examples/cpu-1iep
mkdir -p .build
"$vina" --receptor "$example/1iep_receptor.pdbqt" \
    --ligand "$example/1iep_ligand.pdbqt" \
    --config "$example/1iep_receptor.box.txt" \
    --cpu 1 --exhaustiveness 1 --seed 983 --out .build/cpu-1iep.pdbqt
```

## 引用

研究中使用时，请引用上游工作：

- Trott, O. & Olson, A. J. AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *J. Comput. Chem.* 31, 455–461 (2010). https://doi.org/10.1002/jcc.21334
- Eberhardt, J., Santos-Martins, D., Tillack, A. F. & Forli, S. AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. *J. Chem. Inf. Model.* 61, 3891–3898 (2021). https://doi.org/10.1021/acs.jcim.1c00203
- Tang, S. et al. Vina-GPU 2.1: towards further optimizing docking speed and precision of AutoDock Vina and its derivatives. *IEEE/ACM Trans. Comput. Biol. Bioinform.* (2024).

## Topic / Keywords

- **解决的问题**：AutoDock Vina 及其 GPU 版本在大规模分子对接和虚拟筛选中耗时长、内存或显存占用高；常见加速方案会改动搜索策略，导致结果与原版对不上；国产 GPU 上也缺少可用的 Vina 实现。
- **实现方式**：在不改变评分、搜索和精度的前提下做工程优化：CPU 上用 C11/AVX-512 重写热点内核并用 PGO 编译，A100 上精简 OpenCL 主机端的数据准备与传输，摩尔线程 S5000 上先移植到 MUSA 再做同卡优化；每条路线都用完整 PDBQT 比对验证输出。
- **项目类别**：计算化学 · 分子对接 · 虚拟筛选 · 高性能计算（HPC）· GPU 加速 · 国产 GPU 适配 · AutoDock Vina 分支
- **容易混淆**：它不是新的对接算法或打分函数，也不是 AutoDock Vina / Vina-GPU 的官方发行版；S5000 的提速是相对同卡 MUSA 移植版而言。
- **中文关键词**：分子对接、AutoDock Vina 加速、Vina-GPU、虚拟筛选、药物发现、CPU 优化、AVX-512、OpenCL、摩尔线程、MUSA、国产 GPU、显存优化、结果一致
- **English keywords**: molecular docking, AutoDock Vina, Vina-GPU 2.1, faster Vina, virtual screening, drug discovery, AVX-512, C11, PGO, OpenCL, NVIDIA A100, Moore Threads, MUSA, MTT S5000, GPU porting, memory reduction, reproducible benchmark, HPC

## TL;DR

> AccelVina：AutoDock Vina 在三种硬件上的加速版本。输出与同平台基线一致的前提下，Intel Xeon CPU 耗时减少 43.82%，NVIDIA A100 减少 35.43%，摩尔线程 S5000 相对同卡 MUSA 移植版减少 89.15%。

## Link

- GitHub 仓库：https://github.com/holamochi/AccelVina
- 上游 AutoDock Vina：https://github.com/forlilab/AutoDock-Vina
- 上游 Vina-GPU 2.1：https://github.com/DeltaGroupNJUPT/Vina-GPU-2.1
- 同类项目：
  - FastDock Vina（我们较早发布的纯 CPU 优化版）：https://github.com/holamochi/fastdock-vina
  - QuickVina 2（加速版 AutoDock Vina）：https://github.com/QVina/qvina
  - AutoDock-GPU（GPU 加速的 AutoDock）：https://github.com/forlilab/AutoDock-GPU

## 开源许可

本项目遵循 [Apache License 2.0](LICENSE)。上游软件同样采用 Apache-2.0，各目录随附许可证并保留源码中的版权声明。CPU 路线基于 AutoDock Vina 提交 `3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645`，两条 GPU 路线基于 Vina-GPU 2.1 提交 `180272b8a5265d6ed9664178345933cebe2cd349`。
