# 模型与数据附件

大型文件通过 [GitHub Release](https://github.com/dawnmengsjtu/YeastBridge/releases/tag/reproducibility-assets-20261003)提供，不放入Git历史。文件均从原实验环境恢复，下载器根据 [asset_downloads.json](asset_downloads.json)校验SHA-256。

## 完整筛选所需文件

```bash
python scripts/install_assets.py --download
python main.py --mode check
```

| 文件 | SHA-256 | 用途 |
| --- | --- | --- |
| `raw/tier1_models/route_b/final_model.pt` | `112b39afc687cc59fa54c2de22b5e6c3268c082fc94625740971cdd5759f2be7` | B2投影权重与偏置 |
| `raw/tier1_response/strain_response.npz` | `d174dc0a1cc403b9ab4ce730aa07a11eb51aa290fe04196c2c903d47c486aadf` | HIP/HOP响应矩阵 |

两个文件合计约513 MiB。其他主链输入已在仓库中，共21项主链资产均按冻结哈希检查。离线安装：

```bash
python scripts/install_assets.py --weights /path/to/final_model.pt \
    --response /path/to/strain_response.npz
```

本地安装器先校验所有输入，再复制文件；下载器在临时文件校验通过后安装。不会把不同版本的文件视为冻结输入。PyTorch checkpoint仅加载团队提供且通过校验的原文件。

## 可选训练与预处理输入

- `b2-inputs.tar.gz`：B2的计数矩阵与蛋白矩阵。安装、上游初始权重和训练命令见 [训练说明](../models/TRAINING.md)。
- `expr_full.tsv.gz`：HIP/HOP已处理条形码表达矩阵。重建响应、原始CEL获取方式及R环境见 [数据准备](PREPROCESSING.md)。

B2 checkpoint附带scFoundation [MODEL_LICENSE](../models/training/vendor/scfoundation/MODEL_LICENSE)，用途受其非商业研究条款约束。数据来源与尚未核实的许可记录见 [DATA_SOURCES.md](DATA_SOURCES.md)。

## 文件格式与历史清单

响应文件含 `strain_orfs`、`array_names`、`compound_inchikeys`、`doses`、`dose_units`、`is_vehicle`、`z_score`、`method`。矩阵形状为5,668菌株×3,850条件。`z_score`为历史键名，实际响应定义见 [数据准备](PREPROCESSING.md)。

`A2_init.tsv` 的原冻结哈希对应CRLF换行，`.gitattributes`保留其字节形式。`MANIFEST.sha256` 是历史运行清单；打包工具另为实际打包文件生成 `MANIFEST.release.sha256`，两者用途不同。
