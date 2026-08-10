# MethylVI IR/NR 分析项目

本仓库保存针对 10 个 MethSCAn 样本（5 个 IR、5 个 NR）整理的 ALLCools/MethylVI 分析脚本。

正式脚本和完整运行说明位于：

- [`20260810/scripts/README.md`](20260810/scripts/README.md)

其中包括输入审计、ALLCools 5-kb 上游准备、`mc/cov` 计数重建、MethylVI 训练、UMAP/Leiden 分析和结果绘图。

本仓库不包含原始 `.cov`、`.allc`、H5AD/H5MU 数据、训练结果或日志。运行前请在 `20260810/scripts/mvi_00_config.sh` 中设置服务器上的数据路径。

## 快速入口

```bash
cd 20260810/scripts
bash mvi_05_run_pipeline.sh verify
bash mvi_05_run_pipeline.sh smoke
bash mvi_05_run_pipeline.sh build
bash mvi_05_run_pipeline.sh train
```

正式分析前请先阅读脚本目录中的 README，并确认 MethylVI API、Conda 环境和输入数据目录已经通过审计。
