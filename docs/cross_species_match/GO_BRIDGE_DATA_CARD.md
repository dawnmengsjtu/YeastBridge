# GO-BP / GO-slim 跨物种映射数据卡（v1）

状态：`READY_FOR_LOCKED_BENCHMARK`。这里的“就绪”只表示共享功能投影已经补齐，
不表示人类与酵母联合模型已经胜过任一单独模型。

## 做了什么

映射范围在看任何 benchmark 标签之前冻结为：

- 人类：Norman CRISPRa 数据的全部 3,066 个响应基因符号；
- 酵母：SGA 显著互作表中出现过的 query/array ORF 并集，共 5,848 个；
- 公共空间：GO `biological_process` 的 `is_a` / `part_of` 祖先，最终投影到物种无关的
  `goslim_generic`；
- 默认证据配置保留 IEA，并把 IEA 权重预注册为 0.35；同时给出完全独立的
  `curated_no_iea` 覆盖率审计。没有为了提高覆盖率静默删掉或加入某类证据。

构建器不读取扰动标签、训练/测试划分、模型预测或药物响应。GO 边只负责把两种物种的
基因放入同一个功能坐标系，不能被解释为“同一个 KO”、药效、靶点结合或联合模型增益。

## 实测覆盖率

| 范围 | 总基因 | GO-slim（含 IEA） | 去 IEA | 落入双物种共享模块 |
|---|---:|---:|---:|---:|
| Norman 响应基因 | 3,066 | 2,555（83.33%） | 2,448（79.84%） | 2,458（80.17%） |
| SGA query/array ORF | 5,848 | 4,442（75.96%） | 4,336（74.15%） | 4,375（74.81%） |

共得到 70 个有覆盖的通用 GO-slim BP 模块，其中 52 个在人类和酵母两侧都有基因。
`module_counts.tsv` 逐模块给出人类/酵母基因数，不能把这些基因数当作独立生物样本数。

## 冻结来源

- GO basic ontology：release `2026-07-26`，
  [官方说明](https://geneontology.org/docs/download-ontology/)
- Generic GO slim：release `2026-07-26`，
  [官方说明](https://geneontology.org/docs/go-subset-guide/)
- Human UniProt GAF：generated `2026-07-28`，GO release `2026-07-26`
- Yeast MOD GAF：generated `2026-08-04`，GO release `2026-08-02`
- 下载入口和 2026 年新命名：
  [GO annotation downloads](https://geneontology.org/docs/download-go-annotations/downloads/)
- GO 数据产品许可：
  [CC BY 4.0 与引用政策](https://geneontology.org/docs/go-citation-policy/)

精确 URL、下载时间、字节数、SHA-256、版本和许可记录在
`configs/go_bridge.sources.json`。当前 ontology 与酵母 GAF 相差一个周发布；构建器会计数
无法解析的 GO term，不会静默吞掉。此次构建未观察到无法解析的 term。

## 可提交文件

- `assets/go_bridge_v1/projection.shared.tsv`：benchmark 直接消费的五列标准表：
  `species, gene, module_id, evidence, weight`；
- `assets/go_bridge_v1/gene_goslim_projection.tsv.gz`：每个 gene-module 边的证据码、证据等级、
  直接 GO term、距离、含/不含 IEA 权重；
- `assets/go_bridge_v1/gene_bp_ancestors.tsv.gz`：可审计的 BP 祖先支持明细；
- `assets/go_bridge_v1/module_counts.tsv`：模块覆盖；
- `assets/go_bridge_v1/coverage_report.json`：输入/输出哈希和覆盖报告。

服务器隔离目录：
`/public/home/mengxl/dzy/yeastbridge_re_mvp_stage_20260901/go_bridge_assets_v1`。
原始项目没有被修改。

## 复现命令

```bash
python scripts/build_go_bridge_projection.py \
  --go-basic /path/to/go-basic.obo \
  --goslim /path/to/goslim_generic.obo \
  --human-gaf /path/to/HUMAN-uniprot.gaf.gz \
  --yeast-gaf /path/to/YEAST-mod.gaf.gz \
  --norman-npz /path/to/norman_task.npz \
  --sga-tsv /path/to/sga_significant.tsv.gz \
  --evidence-config configs/go_bridge.evidence.json \
  --output-dir /new/locked/output/directory
```

每次改 ontology、GAF、证据权重或作用域都应使用新的输出目录并重新冻结哈希；不得在看到
test 结果后调 GO-slim、证据等级或祖先关系。
