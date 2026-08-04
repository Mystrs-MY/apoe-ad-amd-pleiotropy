# A1 文献优先循环蛋白面板：检索与筛选策略

版本：1.0
检索日期：2026-07-10
目标：识别已发表的 AD、总体 AMD、干性 AMD/地理萎缩和湿性/新生血管 AMD 的 proteome-wide 或大规模循环蛋白遗传因果研究，用于构建 literature-prioritized circulating-protein panel。该检索不代表 A1 自行完成全蛋白质组扫描。

## 数据库与来源层级

| 来源 | 层级 | 当前用途 | 状态 |
|---|---|---|---|
| PubMed/MEDLINE | T1 | 生物医学主检索、PMID、MeSH 和摘要 | 自动检索 |
| CrossRef | T1 | DOI、正式发表状态和跨期刊补充检索 | 自动检索 |
| 出版社/PMC 正文与补充材料 | T1 | 人工核验蛋白名单、效应、阈值、平台和 coloc/HEIDI | 逐篇核验 |
| OpenAlex | T1 聚合索引 | PubMed/CrossRef 结果不足时的发现补充 | 备用 |
| Semantic Scholar | T2 | 相关文献和引用链补充 | 仅在 T1 不足时 |
| bioRxiv/medRxiv | T2 | 记录预印本及其正式发表 lineage | 不作为默认 Tier 1 |
| Scopus/Web of Science/Embase | 机构依赖 | 人工补检与查重 | 当前环境无稳定机构 API，列入人工核验 |

## 完整检索式

机器可读检索式保存于 literature/search_queries.yml。PubMed 包含四个主查询和两个 serum-protein 补充查询：

1. AD_proteome_MR：AD + MR/SMR/HEIDI/colocalization + proteome/circulating protein/pQTL/platform。
2. AMD_proteome_MR：总体及亚型 AMD + MR/SMR/HEIDI/colocalization + proteome/circulating protein/pQTL/platform。
3. AD_causal_protein_platforms：AD + Olink/SomaScan/deCODE/UKB-PPP/INTERVAL/Fenland/SCALLOP + causal/MR/coloc。
4. AMD_causal_protein_platforms：AMD/GA/nAMD + 同一平台词 + causal/MR/coloc。
5. AD_serum_protein_MR_supplement：AD + MR/SMR/HEIDI/coloc + serum protein/proteome。
6. AMD_serum_protein_MR_supplement：AMD/GA/nAMD + MR/SMR/HEIDI/coloc + serum protein/proteome。

CrossRef 使用四条自然语言补充检索：

- Alzheimer disease proteome-wide Mendelian randomization plasma proteins
- age-related macular degeneration proteome-wide Mendelian randomization plasma proteins
- geographic atrophy circulating proteins Mendelian randomization
- neovascular AMD circulating proteins Mendelian randomization

不设置语言过滤。正式同行评议论文优先；预印本保留 publication lineage，但默认不进入最高可信主分析层。

## 去重规则

1. DOI 统一为小写，移除 https://doi.org/ 前缀和末尾标点后精确去重。
2. DOI 缺失时，按规范化题名 token Jaccard 相似度至少 0.90 且首作者姓氏相同去重。
3. 正式论文优先于预印本；元数据完整记录优先。
4. 文献层去重不等于证据独立。全文提取后继续按 pQTL 数据、outcome GWAS、instrument、platform 和 publication lineage 标记 exact duplicate、partially overlapping、same pQTL different outcome、different pQTL same outcome、independent pQTL and outcome 或 replication unclear。

## 标题摘要筛选

纳入候选：

- 结局明确为 AD、总体 AMD 或指定 AMD 亚型；
- 暴露为循环血浆/血清蛋白；
- proteome-wide 或大规模蛋白候选空间；
- 使用 MR、SMR/HEIDI、遗传预测蛋白分析或等价遗传因果框架；
- 至少可识别 pQTL 来源或平台。

排除：

- 纯观察性蛋白关联；
- 仅 biomarker prediction/机器学习；
- 仅脑脊液或组织蛋白且无循环层；
- 仅代谢物、脂质或细胞计数；
- 非 AD/AMD 结局；
- 综述、社论、方案或无原始分析。

## 全文筛选

Layer 1 收录所有满足疾病、循环蛋白和大规模遗传因果分析的阳性蛋白来源研究。

Layer 2 主面板优先要求：

- multiple-testing corrected；
- cis-pQTL 或 cis-dominant；
- 效应方向与单位可提取；
- 蛋白、gene、UniProt/assay/platform 可可靠映射；
- coloc/HEIDI 或独立复制优先；
- outcome GWAS 定义可追溯；
- 复杂位点有 LD、条件或共定位处理。

Tier 2 包括单平台、无 coloc、mixed cis/trans、trans-only、复制不足、仅 FDR/研究内宽阈值、outcome 不完全匹配或数据追溯不完整。

全文或关键补表不可得、只有 nominal P、映射无法确认、复杂位点无处理、或无法区分循环与组织蛋白者不进入 Layer 2，并保留排除理由。

## 数据提取与人工核验

- 摘要只用于初筛，不从摘要猜测未报告蛋白、效应或 PP.H4。
- 每个结构化字段必须能回溯至正文、补表、数据库或官方平台 manifest。
- full-text/supplement unavailable、not reported、mapping unresolved 和 requires manual verification 均作为显式状态。
- 同一蛋白的多篇报告先生成 study-level records，再生成 evidence summary；不简单累计论文数。

## 当前自动化输出

- literature/search_results_deduplicated.tsv：自动检索和 DOI/题名去重后的候选记录。
- literature/study_screening.tsv：人工筛选状态与排除原因。
- literature/data_extraction_template.tsv：逐研究蛋白提取模板。
- data_raw/literature_search：各数据库原始响应。
- logs/literature_search_run.json：检索日期、查询、计数和错误。

## 当前限制

Scopus、Web of Science 和 Embase 在本环境无稳定机构 API。自动结果将通过 PubMed、CrossRef、PMC/出版社正文和必要时 OpenAlex补充；机构数据库补检列入人工核验清单。任何无法获得关键补表的研究均不会被自动升入 Layer 2。
