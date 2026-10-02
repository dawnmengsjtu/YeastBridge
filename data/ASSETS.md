# 完整运行所需资产

公开仓库中的预计算嵌入、基因表和基准矩阵已由 `scripts/check_assets.py` 按冻结 SHA-256 校验。两个必需文件尚未随公开仓库提供，也没有已核实的公开下载地址：

| 文件 | SHA-256 | 用途 |
| --- | --- | --- |
| `raw/tier1_models/route_b/final_model.pt` | `112b39afc687cc59fa54c2de22b5e6c3268c082fc94625740971cdd5759f2be7` | B2 投影的权重与偏置 |
| `raw/tier1_response/strain_response.npz` | `d174dc0a1cc403b9ab4ce730aa07a11eb51aa290fe04196c2c903d47c486aadf` | HIP/HOP 处理后的响应矩阵 |

取得原文件后执行：

```bash
python scripts/install_assets.py --weights /path/to/final_model.pt \
    --response /path/to/strain_response.npz
python main.py --mode check
```

安装器先校验所有输入，再复制文件；不会将不同版本的资产视为原冻结输入。它不下载或执行外部代码。PyTorch checkpoint 应仅使用团队提供且校验通过的原文件。

响应文件的字段由执行器读取：`strain_orfs`、`compound_inchikeys`、`doses`、`is_vehicle`、`z_score`，以及可选的 `method`。`z_score` 的形状为菌株数×条件数。化合物和剂量选择规则见执行器注释及原结果的 `execute_summary.json`。

原始数据入口是 [E-MTAB-2391](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-2391)。该 accession 不能直接替代项目处理后的响应矩阵；原始数据到冻结矩阵的完整预处理脚本、剂量单位和运行记录仍需从原实验环境补齐。

`A2_init.tsv` 的原冻结哈希对应 CRLF 换行。仓库通过 `.gitattributes` 保留该文件的字节形式，避免系统换行转换导致误报。

`MANIFEST.sha256` 是历史运行清单，保留用于核对旧结果；它不代表当前仓库的发布清单。打包工具会为实际打包文件单独生成 `MANIFEST.release.sha256`。
