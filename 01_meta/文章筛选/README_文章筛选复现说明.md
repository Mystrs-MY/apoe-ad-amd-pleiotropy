# 文章筛选复现说明

本目录用于保存 AMD 与 AD/dementia 系统评价检索、去重、题名摘要筛选、全文候选复核和全文排除理由的可追溯材料。

## 目录结构

| 路径 | 内容 | 用途 |
|---|---|---|
| `01_raw_exports/` | PubMed、Embase、Web of Science、Cochrane CENTRAL 的原始导出文件 | 从数据库导出记录重新开始去重和筛选 |
| `02_dedup_and_screening_outputs/` | 去重表、重复组、题名摘要筛选表、全文评估候选表、Cochrane 补充解析表和最终候选 JSON | 复核记录数、重复数量、去重后数量和 88 篇全文评估候选的形成过程 |
| `03_reports_and_audits/` | 检索去重汇总、88 篇候选复核报告、检索式提取和检索变更审计 | 复核检索策略变更、Cochrane n=79 修正、候选文献与既往纳入研究的关系 |
| `Excluded_Studies_reason_summary.csv` | 全文排除理由分类汇总 | 写作 PRISMA、Methods 或 Supplementary Methods 时使用 |

## 推荐复核顺序

1. 从 `01_raw_exports/` 核对四个主要数据库导出：`PubMed.nbib`、`Embase.ris`、`WOS_1.ris`、`WOS_2.ris`，并结合 Cochrane CENTRAL 手工导出记录。
2. 查看 `02_dedup_and_screening_outputs/dedup_screening_summary.json`，确认可用原始记录、重复记录和去重后记录数量。
3. 使用 `duplicate_groups.csv` 复核重复合并逻辑。
4. 使用 `deduplicated_screened_records*.csv` 和 `full_text_assessment_candidates*.csv` 复核题名摘要筛选与全文候选形成过程。
5. 使用 `_final_candidate_selection.json` 和 `03_reports_and_audits/检索去重与全文评估候选汇总.md` 核对最终 88 篇全文评估候选。
6. 使用 `03_reports_and_audits/88篇候选与既往纳入研究复核报告.md` 核对 88 篇候选与既往 10 篇纳入研究之间的关系。
7. 正式全文排除逐篇理由见 `../tables/TableS6_Excluded_Studies.csv`；排除理由分类汇总见本目录 `Excluded_Studies_reason_summary.csv`。

## 投稿相关说明

- 正式补充表应使用 `<project_root>\01_meta\tables\TableS6_Excluded_Studies.csv`。
- `Excluded_Studies_reason_summary.csv` 是写作和核对用汇总表，不作为单独的 Table S6 提交。
- 本目录不保存桌面临时副本；桌面中的 S6 duplicate 文件不作为项目依据。
