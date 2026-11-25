# AccelVina

我们在 Intel Xeon CPU、NVIDIA A100 和摩尔线程 S5000 上优化了 AutoDock Vina 的分子对接。这里收录三条路线正式验收使用的源码、构建入口、公开复测样例和实验报告。

| 平台与构建入口 | 同平台基线 | 端到端耗时 | 耗时减少 | 峰值内存或显存 | 输出校验 |
| --- | --- | ---: | ---: | ---: | --- |
| [Intel Xeon · C11](cpu/) | AutoDock Vina 上游 | 3445.941 → 1936.001 秒 | **43.82%** | 828.30 → 418.57 MiB | 优化版 75/75 次通过 |
| [NVIDIA A100 · OpenCL](gpu-opencl/) | Vina-GPU 2.1 上游 | 70.042 → 45.226 秒 | **35.43%** | 19.046 → 18.499 GiB | 两版合计 150/150 次通过 |
| [摩尔线程 S5000 · MUSA](gpu-musa/) | 同卡 MUSA 移植版 | 1181.236 → 128.220 秒 | **89.15%** | 3542 → 620 MiB | 两版合计 150/150 次通过 |

耗时是每个计时用例运行三次的中位数之和，CPU 计时覆盖 19 个完整对接配置，两条 GPU 路线各覆盖 25 个最终用例。CPU 测进程驻留内存，A100 显存取公开单例的整卡采样峰值，S5000 显存取正式批次的进程采样峰值。通过数均以对应平台的冻结参考输出为准。

[实验报告](docs/report.md) · [HTML 文件](docs/report.html) · [公开样例与复测](benchmarks/) · [上游来源与许可证](#上游来源与许可证)

## 上游来源与许可证

- CPU 基于 AutoDock Vina 的固定提交 `3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645`。[AutoDock Vina 上游](https://github.com/forlilab/AutoDock-Vina/tree/3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645)
- OpenCL 和 MUSA 均基于 Vina-GPU 2.1 的固定提交 `180272b8a5265d6ed9664178345933cebe2cd349`。共享目录保留构建所需的上游源码及逐文件校验和。[Vina-GPU 2.1 上游](https://github.com/DeltaGroupNJUPT/Vina-GPU-2.1/tree/180272b8a5265d6ed9664178345933cebe2cd349) · [共享上游源码](upstream/vina-gpu-2.1/)

上游软件采用 Apache-2.0。各目录随附许可证并保留源码中的版权声明。[Apache-2.0 许可证](LICENSE)
