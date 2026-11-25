# NVIDIA A100 OpenCL 优化版

正式验收候选为 `a066bd9c7352c72c055179ed4e83bb1644b9225d`，直接基于上游 Vina-GPU 2.1 提交 `180272b8a5265d6ed9664178345933cebe2cd349`。本目录的 90 个源文件、构建文件及许可证与候选提交逐字节一致。候选沿用上游两份 OpenCL 内核，修改主机端网格准备、缓冲传输和批处理流程。[候选源码清单](SOURCE_SHA256SUMS) · [验收版本与 SHA-256](provenance.json)

## 构建

使用 Ubuntu 22.04 x86-64、GCC/G++、GNU Make、NVIDIA 驱动、OpenCL ICD 和 Python 3.10 或更新版本。在 Python 环境安装 `conda-package-handling`，并保证 `python` 和 `python3` 指向该环境。在只可见一张 NVIDIA 卡的环境中，从仓库根目录执行。

```bash
work="$PWD/.build"
archives="$work/archives"
mkdir -p "$archives"
curl -fL https://conda.anaconda.org/conda-forge/linux-64/libboost-1.84.0-h6fcfa73_0.conda \
    -o "$archives/libboost-1.84.0-h6fcfa73_0.conda"
curl -fL https://conda.anaconda.org/conda-forge/linux-64/libboost-headers-1.84.0-ha770c72_0.conda \
    -o "$archives/libboost-headers-1.84.0-ha770c72_0.conda"
curl -fL https://codeload.github.com/KhronosGroup/OpenCL-Headers/tar.gz/refs/tags/v2024.10.24 \
    -o "$archives/OpenCL-Headers-v2024.10.24.tar.gz"

base="$work/a100-base"
bash build_baseline.sh "$base" "$archives"
make -C gpu-opencl -j4 BASE_BUILD="$base" BUILD_DIR="$work/a100-opt"
```

基线脚本校验三份依赖和共享上游源码，生成当前驱动可加载的两份 OpenCL 内核，下载并校验正式基线所用的上游 Linux 可执行文件。`$base` 须尚未存在。正式环境为 A100-SXM4-80GB，NVIDIA 驱动 535.129.03，冻结内核目标为 `sm_80`。

## 公开样例

复用刚构建的运行库和内核，在同一张卡上运行 2bm2 并比较输出。

```bash
runtime="$base/runtime"
LD_LIBRARY_PATH="$runtime/lib:/usr/lib/x86_64-linux-gnu" \
OCL_ICD_VENDORS="$runtime/icd" \
python3 gpu-musa/bench/run_benchmark.py --suite benchmarks/gpu-2bm2.json \
    --backend opencl --gpu "NVIDIA A100" \
    --baseline "$runtime/AutoDock-Vina-GPU-2-1" \
    --candidate "$work/a100-opt/AutoDock-Vina-GPU-2-1" \
    --baseline-kernels "$runtime" --candidate-kernels "$runtime" \
    --repeats 1 --output "$work/a100-2bm2"
```

正式计时使用 25 个最终用例和三轮复测。[实验报告](../docs/report.md)
