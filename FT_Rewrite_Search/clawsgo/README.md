# ClawsGO 的材料（与本地工作分开）

这个文件夹只放 ClawsGO 产出的东西和发给它的东西。**本地自己做、已核对的工作**在项目根目录：`code/`、`results/`、`experiments/LOG.md`、`lessons/`、`notes/`。ClawsGO 的结果在核对之前，不并入这些目录，也不引用为结论。

| 子目录 | 内容 |
|---|---|
| `reports/` | ClawsGO 交回的报告：`Literature_review.pdf`（L4 调研，10-28）、`实验验证.pdf`（实验方案，10-28）、`实验报告.pdf`（试点结果，10-29） |
| `prompts/` | 发给 ClawsGO 的 prompt |
| `handoff/` | 交给 ClawsGO 的材料包：`handoff_clawsgo_2026-10-29/` 及同名 zip（补充实验计划、代码快照、LOG、标准答案、参考资料） |
| `returns/` | ClawsGO 交回的结果包，按日期建子文件夹原样存放，不修改。已有：`2026-10-29_supplementary/`（补充实验，未核对） |

## 核对流程（结果交回后）

1. 把交回的包原样解压到 `returns/<日期>/`。
2. 逐项核对，记录写到 `notes/` 下（如 `ClawsGO_return_verification_<日期>.md`）：
   - 回归测试：在本地用交回的代码重跑 E1.2、E1.5，结果必须和 `results/` 的标准答案一致；
   - 每个实验：对照计划的通过标准，抽查关键数字（本地重跑或独立复算）；
   - 代码改动：和 `handoff/` 里的代码快照做 diff，确认没有另起一套、规则仍按 Jeon 第 III 节；
   - 引用和数字：核对 arXiv ID，标出推断。
3. 核对通过的部分才合并：代码进 `code/`，结果进 `results/`，在 `experiments/LOG.md` 里注明"来自 ClawsGO，<日期> 核对"。

已有的核对记录：`notes/Literature_review_verification_2026-10-28.md`、`notes/Experiment_design_verification_2026-10-28.md`、`notes/Pilot_report_and_supplement_review_2026-10-29.md`。
