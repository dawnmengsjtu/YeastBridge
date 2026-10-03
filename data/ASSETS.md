# 模型与数据附件

当前版本见 [submission-v2发布页](https://github.com/dawnmengsjtu/YeastBridge/releases/tag/submission-v2-20261003)。大型checkpoint及响应矩阵通过Release分发，安装器读取 [asset_downloads.json](asset_downloads.json)校验SHA-256。

## 当前完整筛选

```bash
python scripts/install_assets.py --download
python main.py --mode check
```

| 文件 | SHA-256 |
| --- | --- |
| `raw/tier1_models/route_b_v2/final_model.pt` | `281d363fd3ac831f692bff7f3e3a0ae7420b3a221792591c78e5130de09c27ad` |
| `raw/tier1_response/strain_response.npz` | `d174dc0a1cc403b9ab4ce730aa07a11eb51aa290fe04196c2c903d47c486aadf` |

其余当前输入已在仓库中，按 `configs/esm2_joint_tasks_v2.json` 与 `configs/b2_submission_assets.json` 的实际依赖检查。安装器拒绝错误哈希，本地输入先核验再安装。离线安装：

```bash
python scripts/install_assets.py --weights /path/to/b2-submission-v2-final_model.pt \
    --response /path/to/strain_response.npz
```

[完整提交包](https://github.com/dawnmengsjtu/YeastBridge/releases/download/submission-v2-20261003/YeastBridge-full.zip)已包含上述输入和模型，无需重复下载。源码仓库中的数据处理与模型配置均使用相对路径。

## 训练与预处理

`b2-inputs.tar.gz`含B2计数与蛋白矩阵，细胞顺序另存于训练证据中。完整包包含该归档；安装与官方scFoundation初始权重获取见 [TRAINING.md](../models/TRAINING.md)。`expr_full.tsv.gz`用于重新构建HIP/HOP响应，数据获取、原始CEL及R环境见 [PREPROCESSING.md](PREPROCESSING.md)。

scFoundation初始权重通过官方来源获取，不重复打入完整包；提交的最终B2权重保留上游[非商业研究许可](../models/training/vendor/scfoundation/MODEL_LICENSE)。各数据集条款见 [THIRD_PARTY_DATA.md](THIRD_PARTY_DATA.md)。

## 历史版本与文件格式

9月B2仍在[旧资产发布页](https://github.com/dawnmengsjtu/YeastBridge/releases/tag/reproducibility-assets-20261003)，可用 `python scripts/install_assets.py --legacy --download` 安装至旧路径。旧配置仍显式指向旧权重；当前默认入口读取 `models/submission.json`。

HIP/HOP响应含5,668菌株×3,850条件及剂量单位。`z_score`是沿用的历史键名，实际响应定义见数据准备文档。`A2_init.tsv`为基因顺序表，原CRLF字节由 `.gitattributes` 保留。

`MANIFEST.sha256`是历史资产目录；包内 `MANIFEST.release.sha256`核验本次实际分发的所有文件，`BUNDLE.json`记录代码版本与包类型。
