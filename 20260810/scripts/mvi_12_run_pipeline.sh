#!/usr/bin/env bash

# MethylVI 主流程的 dsub 任务入口。
# CPU、内存、工作目录和日志由外层 dsub 命令设置。

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
bash "$SCRIPT_DIR/mvi_04_run_pipeline.sh" all
