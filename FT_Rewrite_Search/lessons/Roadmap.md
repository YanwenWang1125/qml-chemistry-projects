# 阶段 0 课程路线图：FT_Rewrite_Search

目标：11-08 前通过阶段 0 验收（L04 的 gate）。按 [learning-style.md](../../../Plan/learning-style.md) 自上而下学：先看懂要改的电路，再学工具，最后学改写规则。项目计划见 [FT_Rewrite_Search_Plan.md](../FT_Rewrite_Search_Plan.md)。

讲义按 `Quantum Learning/.claude/skills/teach-section` 的流程和模板写成：先写验证代码并实际运行，讲义里只放真实输出。环境：`C:/Users/Administrator/venvs/qml/Scripts/python.exe`，在各课的 `teaching/code/` 目录下运行脚本。共用工具 `ftlib.py` 放在 L01 的 code 目录里，后面几课通过路径引用它。

建议的阅读顺序（每课相同）：讲义 → 原文对应页 → 自测 → gate。

| 课 | 时间 | 主题 | 原文 | 代码与验证 | Gate | 讲义 |
|---|---|---|---|---|---|---|
| L01 | 10/28 – 10/30 | 容错条件与 Shor 型 syndrome 提取 gadget | Jeon 第 II A、II D 节（PDF 第 2–3 页），图 1（第 6 页） | `ftlib.py`、`L01_check_gadgets.py`：三个 gadget 的无噪声正确性、资源数、全部单故障的 $\mathrm{FT}_1$ 检查 | 说明 cat 态和 flag 的作用、$\mathrm{FT}_1$ 条件、无 flag 版本违反在哪 | [L01_teaching.md](L01_FT_and_Shor_Gadget/teaching/L01_teaching.md) |
| L02 | 10/30 – 11/02 | 用 Stim 做带噪声模拟：$p$ 与 $p^2$ 的 scaling | Jeon 第 V A 节（噪声模型）；Stim 文档 | `L02_noisy_gadgets.py`：每点 $10^6$ shots 的蒙特卡洛；`L02_first_order.py`：手算一阶系数 | 解释 $\mathrm{FT}_1 \Rightarrow p^2$，以及模拟和论文逻辑错误率的差别 | [L02_teaching.md](L02_Stim_Noisy_Gadgets/teaching/L02_teaching.md) |
| L03 | 11/02 – 11/05 | Fault equivalence | 2506.17181 第 3.1 节（定义 3.1）、定义 3.4、定理 3.11、命题 3.12 | `L03_fault_equivalence.py`：四组电路的单故障等价检查 | 定义；为什么保持等价就保持容错；为什么错误率可以不同 | [L03_teaching.md](L03_Fault_Equivalence/teaching/L03_teaching.md) |
| L04 | 11/05 – 11/08 | Jeon 的四条改写规则与 Bell reduction | Jeon 第 III、IV 节（PDF 第 3–5 页），附录 B（第 13 页） | `L04_bell_reduction.py`：最小例子的 Bell reduction，以及 Shor gadget 前后对比 | **阶段 0 验收**：手推一次 Bell reduction，并说明两个条件和"等价但错误率不同" | [L04_teaching.md](L04_Jeon_Rewrite_Rules/teaching/L04_teaching.md) |

## 已核对的事实（来自代码输出）

- 论文图 1 的两个电路：基线 5 个辅助比特、9 个 CNOT、深度 5；优化后 3、7、6；与表 I 一致。两者都满足 $\mathrm{FT}_1$。去掉 flag 的版本有 17 个违规单故障。
- 触发 flag 的单故障：我们数到 43 和 34 个，论文表 II 是 42 和 33 个，各差 1 个，原因未核实（可能是 flag 比特 SPAM 错误的建模不同）。
- 蒙特卡洛：无 flag 版本斜率约 1，手算一阶系数 $c_1 \approx 2.07$ 与模拟吻合；基线和优化后电路斜率约 1.86 和 1.96。
- 基线与优化后电路、Bell reduction 最小例子的 A 与 B：都通过单故障等价检查，但错误率都相差约 2.4 倍。

## 尚未覆盖（阶段 1 再学）

- ZX-calculus 本身，以及规则的 ZX 证明（Jeon 附录 A）。
- 空闲噪声、$Z$ 型 gadget、完整的多轮协议和 flag 条件查表解码。
- paritea 的使用。
