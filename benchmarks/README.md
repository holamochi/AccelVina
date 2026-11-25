# 公开样例与复测

| 样例 | 固定来源 | 运行入口 |
| --- | --- | --- |
| [1iep](examples/cpu-1iep/) | AutoDock Vina 提交 `3c65c0b` 的公开入门样例 | [CPU C11](../cpu/#公开样例) |
| [2bm2](examples/gpu-2bm2/) | Vina-GPU 2.1 提交 `180272b` 的公开入门样例 | [A100 OpenCL](../gpu-opencl/#公开样例) · [S5000 MUSA](../gpu-musa/#公开样例) |

[2bm2 清单](gpu-2bm2.json) 固定两份输入的 SHA-256、对接盒、种子和设置。GPU 复测脚本使用 Python 3.10 或更新版本，在同一设备上比较基线和候选的完整 PDBQT 文件和姿态表。使用新输出目录，并将 `--repeats 1` 改为 `--repeats 3` 可以检查三轮输出与汇总耗时。

正式结果及其计时样本范围、峰值内存或显存的采样范围见实验报告。正式最终用例的逐例结果由评测方保管。[实验报告](../docs/report.md)
