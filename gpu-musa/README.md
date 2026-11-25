# 摩尔线程 S5000 MUSA 优化版

从相同的 Vina-GPU 2.1 上游源码构建同卡 MUSA 移植基线，再按正式验收顺序应用三份优化补丁。补丁涉及主机数据准备、显存映射和网格诊断。警告回退时重放原网格设备代码，搜索内核沿用移植基线。[冻结源码与补丁清单](report/evidence/freeze.json) · [正式计时汇总](report/evidence/final-summary.json)

## 构建

使用 Linux x86-64、GNU Make、GNU patch、GCC/G++ 和 MUSA SDK 4.3.8。依赖准备使用 curl、unzip、zstd 和 GNU tar。正式验证卡为 MTT S5000，目标架构为 `mp_31`。从仓库根目录执行，将 `MUSA_HOME` 指向本机的 SDK。

```bash
work="$PWD/.build"
mkdir -p "$work"
export MUSA_HOME=/usr/local/musa-4.3.8
upstream="$PWD/upstream/vina-gpu-2.1"
bash gpu-musa/prepare_boost.sh "$work/musa-boost"
bash gpu-musa/build.sh baseline "$work/musa-base" "$work/musa-boost" "$upstream"
bash gpu-musa/build.sh candidate "$work/musa-opt" "$work/musa-boost" "$upstream"
```

各构建输出目录应为新目录，`build.sh` 为各自的程序、内核对象和运行库生成 `SHA256SUMS`。

## 公开样例

复测脚本需要 Python 3.10 或更新版本。将 `MUSA_VISIBLE_DEVICES` 设为本机测试卡的编号。可以用 `VINA_EXPECT_PCI` 指定卡的实际 PCI 地址来校验选卡。

```bash
MUSA_VISIBLE_DEVICES=0 LD_LIBRARY_PATH="$MUSA_HOME/lib" \
python3 gpu-musa/bench/run_benchmark.py --suite benchmarks/gpu-2bm2.json \
    --backend musa --gpu "MTT S5000" \
    --baseline "$work/musa-base/AutoDock-Vina-GPU-2-1" \
    --candidate "$work/musa-opt/AutoDock-Vina-GPU-2-1" \
    --repeats 1 --output "$work/musa-2bm2"
```

S5000 优化版以同一张 S5000 上的 MUSA 移植版为输出参考。正式计时覆盖 25 个最终用例，各运行三次。[实验报告](../docs/report.md)
