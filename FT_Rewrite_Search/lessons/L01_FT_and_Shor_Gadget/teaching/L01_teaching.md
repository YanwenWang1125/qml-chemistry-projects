# L01 容错条件与 Shor 型 syndrome 提取 gadget

材料：Jeon et al., *Automated reduction of fault-tolerant circuits*, arXiv:2610.09749，第 II A、II D 节（PDF 第 2–3 页）和图 1（PDF 第 6 页，渲染图见 [figures/jeon2026_page6_fig1.png](figures/jeon2026_page6_fig1.png)）。本地 PDF：`Plan/literature_sources_2026-10-21/jeon2026.pdf`。

代码：[code/ftlib.py](code/ftlib.py)（后面几课共用的故障传播工具），[code/L01_check_gadgets.py](code/L01_check_gadgets.py)，输出 [code/L01_output.txt](code/L01_output.txt)。

## 1. 这一课讲什么

本项目要改写的对象，是"测量一个稳定子"的小电路，叫 gadget。这一课先把它看懂：它测什么、为什么要用 cat 态、为什么还要加一个 flag 比特，以及"容错"在这里的精确定义（$\mathrm{FT}_1$）。最后用代码对论文图 1 的两个电路（基线电路和 Jeon 优化后的电路）逐个检查所有单故障，确认它们都满足 $\mathrm{FT}_1$，而去掉 flag 的版本不满足。

和已学内容的关系：这一课需要稳定子、Steane 码、syndrome 的基本概念（`Quantum Learning/QEC_Prereq` 里已讲过）。

## 2. 符号表

| 符号 | 含义 | 容易混的地方 |
|---|---|---|
| $[[7,1,3]]$ | Steane 码：7 个数据比特，1 个逻辑比特，码距 3 | 码距 $d=3$，最多纠 $t=1$ 个错误 |
| $q_0,\dots,q_6$ | 数据比特 | |
| $q_7,q_8,\dots$ | 辅助比特（ancilla），编号和论文附录 B 一致 | |
| $S$ | 稳定子群，$2^6=64$ 个元素 | 生成元只有 6 个，群元素有 64 个 |
| $X_{0123}$ | 简写，表示 $X_0X_1X_2X_3$ | 这一课的 gadget 测的就是它 |
| $t$ | 能纠正的错误数，$t=\lfloor (d-1)/2 \rfloor$ | Steane 码 $t=1$ |
| $F$ | 一个故障模式（fault）：电路里哪些位置出了错 | **fault 发生在电路里** |
| $\mathrm{wt}(F)$ | 故障事件的个数 | 一个两比特门后面的 $X\otimes X$ 也只算 1 个故障 |
| $E_G(F)$ | 故障 $F$ 跑完 gadget $G$ 后，留在数据比特上的 Pauli 错误 | **error 留在数据上** |
| $\mathrm{wt}_S(E)$ | $E$ 乘以任意稳定子后能达到的最小权重 | 不是 $E$ 本身的权重 |
| flag | 专门用来报警的辅助比特，测量结果为 1 表示"可能出事了" | 不参与算 syndrome |

## 3. 概念：先白话，再定义

**Steane 码的稳定子。** 6 个生成元：3 个 $X$ 型、3 个 $Z$ 型，作用在相同的三组比特上（论文表 II 的标题给出了这三组）：

$$
\{0,1,2,3\},\quad \{1,2,4,5\},\quad \{2,3,5,6\}
$$

测稳定子不能直接测数据（会破坏量子态），所以要借助辅助比特，把"这个稳定子是 $+1$ 还是 $-1$"读到辅助比特上。

**最朴素的做法，以及它的问题。** 一个辅助比特制备成 $\lvert + \rangle$，依次对 $q_0,q_1,q_2,q_3$ 做 CNOT（辅助比特是控制），最后在 $X$ 基测量。问题在于：如果辅助比特在第一个和第二个 CNOT 之间出了一个 $X$ 错误，这个错误会被后面三个 CNOT 复制到 $q_1,q_2,q_3$ 上。一个故障变成了数据上的多比特错误（这叫 hook error），而码距 3 的码只能纠 1 个。

**Shor 的做法：cat 态。** 用 4 个辅助比特，先做成 cat 态

$$
\lvert \mathrm{cat} \rangle = \frac{1}{\sqrt{2}}\left(\lvert 0000 \rangle + \lvert 1111 \rangle\right)
$$

然后第 $j$ 个辅助比特只对第 $j$ 个数据比特做 CNOT。这样每个辅助比特只碰一个数据比特，辅助比特上的单个 $X$ 错误最多传到一个数据比特。

**新问题：制备 cat 态的电路本身也会传播错误。** 制备 cat 态要用 CNOT 树，树上的一个故障可能让两个 cat 比特同时翻转，接下来就会传到两个数据比特上。所以要在 cat 态碰数据之前，先检查一下它是不是好的：用一个 flag 比特测两个 cat 比特的 $Z$ 奇偶。好的 cat 态奇偶一定是 0；如果是 1，说明制备时出了错。

**正式定义：$\mathrm{FT}_t$**（论文式 (2)，Jeon 引用 Chamberland–Beverland 2018）。一个 syndrome 提取 gadget $G$ 是 $t$-容错的，如果对任意故障模式 $F$：

$$
\left(\mathrm{wt}(F) \le t\right) \wedge \left(F \text{ 没有被 flag 检测到}\right) \;\Rightarrow\; \mathrm{wt}_S\!\left(E_G(F)\right) \le \mathrm{wt}(F)
$$

其中

$$
\mathrm{wt}_S(E) = \min_{s \in S} \mathrm{wt}(E s)
$$

读法：不超过 $t$ 个故障、而且 flag 没报警的时候，数据上留下的错误（按稳定子等价算）不能比故障数多。对 Steane 码 $t=1$，就是"任何一个没被发现的单故障，最多在数据上留下 1 比特的错误"。

为什么要"按稳定子等价算"：差一个稳定子的两个错误，对码空间的作用完全一样（论文第 II A 节）。比如数据上出现 $X_0X_1X_2X_3$，它本身就是稳定子，对逻辑信息没有任何影响，权重按 0 算。

## 4. 一步步推导

### 4.1 为什么 cat 比特 $X$ 测量结果的奇偶等于稳定子的值

第 1 步：设数据处在 $X_{0123}$ 的本征态，$X_{0123}\lvert \psi \rangle = \lambda \lvert \psi \rangle$，$\lambda = \pm 1$。

第 2 步：cat 态的第 $j$ 个比特控制数据 $q_j$ 上的 $X$。cat 态里 $\lvert 0000 \rangle$ 分量不动数据，$\lvert 1111 \rangle$ 分量对数据作用 $X_{0123}$：

$$
\frac{1}{\sqrt{2}}\left(\lvert 0000 \rangle\lvert \psi \rangle + \lvert 1111 \rangle X_{0123}\lvert \psi \rangle\right)
= \frac{1}{\sqrt{2}}\left(\lvert 0000 \rangle + \lambda \lvert 1111 \rangle\right)\lvert \psi \rangle
$$

读法：数据没有变，稳定子的值 $\lambda$ 被"踢"到了 cat 态的相对相位上（phase kickback）。

第 3 步：在 $X$ 基测 4 个 cat 比特，结果 $m_j \in \{0,1\}$。$X^{\otimes 4}$ 作用在 $\lvert 0000 \rangle + \lambda \lvert 1111 \rangle$ 上：

$$
X^{\otimes 4}\left(\lvert 0000 \rangle + \lambda \lvert 1111 \rangle\right) = \lvert 1111 \rangle + \lambda \lvert 0000 \rangle = \lambda\left(\lvert 0000 \rangle + \lambda \lvert 1111 \rangle\right)
$$

所以这个态是 $X^{\otimes 4}$ 的本征态，本征值为 $\lambda$。4 个单比特 $X$ 测量结果的乘积等于 $X^{\otimes 4}$ 的值，即 $(-1)^{m_7+m_8+m_9+m_{10}} = \lambda$。

结论：**cat 比特测量结果的奇偶就是 syndrome 比特。** 单个 $m_j$ 是随机的，只有奇偶有意义。

### 4.2 Pauli 错误怎么穿过 CNOT 和 H

后面所有的故障分析都只用这几条规则（标准结论，教材没写）：

| 门 | 进去的错误 | 出来的错误 |
|---|---|---|
| CNOT，控制 $c$，目标 $t$ | $X_c$ | $X_c X_t$（$X$ 从控制复制到目标） |
| | $Z_t$ | $Z_c Z_t$（$Z$ 从目标复制到控制） |
| | $X_t$ 或 $Z_c$ | 不变 |
| H | $X$ | $Z$ |
| | $Z$ | $X$ |

$Y = iXZ$ 按 $X$ 和 $Z$ 两部分分别传。在 $Z$ 基测量时，带 $X$ 或 $Y$ 的比特结果会翻转，带 $Z$ 的不翻转。

### 4.3 基线电路（图 1a）逐门读

辅助比特 $q_7,\dots,q_{10}$ 做 cat 态，$q_{11}$ 是 flag。

| 步 | 门 | 作用 |
|---|---|---|
| 1 | $H$ on $q_7$ | $q_7$ 变成 $\lvert + \rangle$ |
| 2 | CX $7\to 8$ | |
| 3 | CX $7\to 9$，CX $8\to 10$ | 三步做出 4 比特 cat 态 |
| 4 | CX $7\to 11$ | flag 开始测 $Z_7 Z_8$ |
| 5 | CX $7\to 0$，CX $8\to 11$ | flag 测完；$q_7$ 开始碰数据 |
| 6 | CX $8\to 1$，CX $9\to 2$，CX $10\to 3$ | 每个 cat 比特只碰一个数据比特 |
| 7 | $H$ on $q_7..q_{10}$，测量全部 | cat 比特在 $X$ 基测量，flag 在 $Z$ 基测量 |

flag 为什么在好的 cat 态下一定是 0：CX $7\to 11$ 和 CX $8\to 11$ 把 $q_7$、$q_8$ 的 $Z$ 值加到 $q_{11}$ 上，flag 的结果是 $z_7 \oplus z_8$。cat 态里 $z_7 = z_8$（要么全 0，要么全 1），所以是 0。

### 4.4 一个具体故障：$q_8$ 的制备错误

故障：$q_8$ 制备成 $\lvert 0 \rangle$ 后马上出了一个 $X$，即 $X_8$。

- 经过 CX $7\to 8$：$X_8$ 在目标上，不变。
- 经过 CX $8\to 10$：$X_8$ 在控制上，复制到目标，变成 $X_8 X_{10}$。
- 经过 CX $8\to 11$：$X_8$ 又复制到 flag，变成 $X_8 X_{10} X_{11}$。flag 在 $Z$ 基测量，$X_{11}$ 让它翻转为 1，**被检测到**。
- 经过 CX $8\to 1$、CX $10\to 3$：数据上留下 $X_1 X_3$。

**如果没有 flag**，这个故障就不会被发现，数据上留下 $X_1X_3$。按稳定子等价算：

$$
X_1X_3 \cdot X_0X_1X_2X_3 = X_0X_2
$$

权重还是 2。和其余 6 个 $X$ 型稳定子相乘只会更重。所以 $\mathrm{wt}_S(X_1X_3) = 2 > 1$，违反 $\mathrm{FT}_1$。

## 5. 完整数值例子：手算资源，和代码对照

**基线电路**：辅助比特 4 个 cat 加 1 个 flag，共 5 个。CNOT：cat 态 3 个，flag 2 个，碰数据 4 个，共 9 个。CNOT 深度按 4.3 节的分层：第 2 步 1 层，第 3 步 1 层，第 4 步 1 层，第 5 步 1 层，第 6 步 1 层，共 5 层。

**Jeon 优化后的电路（图 1b）**：辅助比特 $q_7$（$\lvert + \rangle$，$X$ 测量）、$q_8$（$\lvert 0 \rangle$，$X$ 测量）、$q_9$（flag），共 3 个。门的顺序（从图上读出）：

$$
\mathrm{CX}_{7\to 2},\ \mathrm{CX}_{7\to 9},\ \mathrm{CX}_{7\to 0},\ \mathrm{CX}_{7\to 8},\ \mathrm{CX}_{7\to 9},\ \mathrm{CX}_{8\to 3},\ \mathrm{CX}_{7\to 1}
$$

共 7 个 CNOT。前 6 个里除 $\mathrm{CX}_{8\to 3}$ 外都经过 $q_7$，必须一个接一个；$\mathrm{CX}_{8\to 3}$ 可以和最后的 $\mathrm{CX}_{7\to 1}$ 并行。所以深度是 6。

代码输出的资源数（见第 6 节）：基线 5 / 9 / 5，优化后 3 / 7 / 6，和手算一致，也和论文表 I 一致。

注意优化后的电路里，$q_7$ 碰了 3 个数据比特，flag 由 $q_7$ 自己的两个 CNOT 夹住 $\mathrm{CX}_{7\to 0}$ 和 $\mathrm{CX}_{7\to 8}$。为什么这样仍然容错，见自测第 5 题。

## 6. 代码结果

运行：

```text
C:/Users/Administrator/venvs/qml/Scripts/python.exe L01_check_gadgets.py
```

真实输出：

```text
========================================================================
baseline without verification
  resources: {'ancillas': 4, 'cnots': 7, 'cnot_depth': 3}
  noiseless, data |+>^7 (X0X1X2X3 = +1): syndrome parity seen [0], flag seen [0]
  noiseless, data Z0|+>^7 (X0X1X2X3 = -1): syndrome parity seen [1], flag seen [0]
  single faults: 128, raise the flag: 0, undetected: 128
  undetected faults by data weight mod S: {0: 45, 1: 66, 2: 17}
  FT_1: False
    violation: prep X on q8                             data error X1X3         weight mod S = 2
    violation: IX on (q7,q8) after CX7->8               data error X1X3         weight mod S = 2
    violation: IY on (q7,q8) after CX7->8               data error X1X3         weight mod S = 2
    violation: XI on (q7,q8) after CX7->8               data error X0X2         weight mod S = 2
    violation: XZ on (q7,q8) after CX7->8               data error X0X2         weight mod S = 2
    violation: YI on (q7,q8) after CX7->8               data error X0X2         weight mod S = 2
    ... 17 violations in total
========================================================================
baseline (Fig. 1a)
  resources: {'ancillas': 5, 'cnots': 9, 'cnot_depth': 5}
  noiseless, data |+>^7 (X0X1X2X3 = +1): syndrome parity seen [0], flag seen [0]
  noiseless, data Z0|+>^7 (X0X1X2X3 = -1): syndrome parity seen [1], flag seen [0]
  single faults: 160, raise the flag: 43, undetected: 117
  undetected faults by data weight mod S: {0: 51, 1: 66}
  FT_1: True
========================================================================
reduced (Fig. 1b)
  resources: {'ancillas': 3, 'cnots': 7, 'cnot_depth': 6}
  noiseless, data |+>^7 (X0X1X2X3 = +1): syndrome parity seen [0], flag seen [0]
  noiseless, data Z0|+>^7 (X0X1X2X3 = -1): syndrome parity seen [1], flag seen [0]
  single faults: 120, raise the flag: 34, undetected: 86
  undetected faults by data weight mod S: {0: 31, 1: 55}
  FT_1: True
```

每项检查说明了什么：

| 检查 | 说明了什么 |
|---|---|
| noiseless 两行 | 无噪声时三个电路都在正确测 $X_{0123}$：本征值 $+1$ 时奇偶为 0，$-1$ 时为 1；flag 始终为 0 |
| resources | 基线 5 / 9 / 5、优化后 3 / 7 / 6，和论文表 I 一致 |
| 无 flag 版本：17 个违规 | 4.4 节手算的 $q_8$ 制备错误是其中第一条；其余都来自 cat 态制备树上的故障 |
| 基线和优化后：FT_1 True | 所有未被检测的单故障都只留下权重 $\le 1$ 的错误 |

**这些检查的边界：**
- 故障模型是代码里定义的这一套：制备后 $X$、H 后 $X/Y/Z$、CNOT 后 15 种两比特 Pauli、测量前 $X$；没有空闲错误。换一套故障模型，计数会变。
- 只检查了测 $X_{0123}$ 的一个 gadget，没有检查 $Z$ 型 gadget 和一整轮 6 个 gadget。
- **和论文对不上的一处**：触发 flag 的单故障，我们数到基线 43 个、优化后 34 个；论文表 II 是 42 和 33 个，各差 1 个。可能的原因是对 flag 比特 SPAM 错误的建模不同（我们把制备错误和测量错误各算一次）。这一点没有核实，请对照论文第 V A 节。
- 两个电路的接线是从论文图 1 读出来的，请对照 PDF 第 6 页核对一遍。

## 7. 常见误区

| 错误想法 | 为什么容易这样想 | 正确的理解 |
|---|---|---|
| 门越少越好 | 少一个门就少一个出错位置 | 无 flag 版本门更少，但不容错；要在容错的前提下减少资源 |
| 故障权重等于数据错误的权重 | 两个词都叫"权重" | $\mathrm{wt}(F)$ 数电路里出错的位置，$\mathrm{wt}(E)$ 数数据上出错的比特；容错条件比的正是这两者 |
| 数据上出现 $X_0X_1X_2X_3$ 是 4 比特错误 | 看起来动了 4 个比特 | 它是稳定子，$\mathrm{wt}_S = 0$，对逻辑信息没有影响 |
| CNOT 只会把错误从控制传到目标 | 只记住了 $X$ 的规则 | $Z$ 反过来，从目标传到控制 |
| flag 报警就要丢掉这次结果 | 很多态制备协议确实这样做 | Jeon 的协议里 flag 报警是用来切换解码表的，不丢 shot（第 V B 3 节） |
| 单个 cat 比特的测量结果就是 syndrome | 4 个结果看起来各自有意义 | 单个结果是随机的，只有奇偶有意义 |

## 8. 自测

**1.** 对 Steane 码，求 $\mathrm{wt}_S(X_0X_1X_2)$。

<details><summary>答案</summary>

$X_0X_1X_2 \cdot X_0X_1X_2X_3 = X_3$，权重 1。和其他稳定子相乘不会更低（不可能是 0，因为 $X_0X_1X_2$ 不是稳定子：它和 $Z_{2356}$ 只在 $q_2$ 重叠，反对易）。所以 $\mathrm{wt}_S = 1$。
</details>

**2.** 基线电路里，flag $q_{11}$ 在 CX $7\to 11$ 和 CX $8\to 11$ 之间出了一个 $Z$ 错误。它最后造成什么影响？

<details><summary>答案</summary>

$Z_{11}$ 在 CX $8\to 11$ 的目标上，复制到控制：$Z_8 Z_{11}$。$Z_{11}$ 之后没有门了，flag 在 $Z$ 基测量，不翻转。$Z_8$ 经过 CX $8\to 1$（$Z$ 在控制上，不传到目标）不变，再经过 $H_8$ 变成 $X_8$，于是 $q_8$ 的测量结果翻转。syndrome 奇偶翻转，数据上没有错误。这是一个"测错 syndrome"的故障，不违反 $\mathrm{FT}_1$；它靠多轮重复测量来处理。
</details>

**3.**（边界情况）$q_7$ 制备成 $\lvert 0 \rangle$ 后马上出了一个 $X$（在 $H$ 之前）。会怎样？

<details><summary>答案</summary>

$X_7$ 经过 $H_7$ 变成 $Z_7$，即 $\lvert + \rangle$ 变成了 $\lvert - \rangle$。$Z_7$ 在之后所有 CNOT 里都是控制上的 $Z$，不传播。最后 $H_7$ 又把它变回 $X_7$，$q_7$ 的测量翻转，syndrome 翻转，数据和 flag 都不受影响。效果和第 2 题一样，只是测错一位 syndrome。
</details>

**4.**（新情形）要测一个权重为 6 的 $X$ 型稳定子，用 Shor 型 gadget 需要几个 cat 比特、几个碰数据的 CNOT？为什么每个 cat 比特只能碰一个数据比特？

<details><summary>答案</summary>

6 个 cat 比特，6 个碰数据的 CNOT（cat 态本身的制备和验证另算）。每个 cat 比特只碰一个数据比特，是为了让一个 cat 比特上的单个 $X$ 错误最多传到一个数据比特。如果一个辅助比特碰多个数据比特，它中途出错就会传到多个数据比特（hook error），除非另外用 flag 把这种情况抓住。
</details>

**5.** 优化后的电路里，$q_7$ 在 $\mathrm{CX}_{7\to 2}$ 之后、第一个 $\mathrm{CX}_{7\to 9}$ 之前出了一个 $X$。flag 抓不到它（为什么？），但电路仍然满足 $\mathrm{FT}_1$。请说明原因。

<details><summary>答案</summary>

$X_7$ 依次经过：$\mathrm{CX}_{7\to 9}$（flag 翻转一次），$\mathrm{CX}_{7\to 0}$（数据 $X_0$），$\mathrm{CX}_{7\to 8}$（$X_8$），$\mathrm{CX}_{7\to 9}$（flag 再翻转一次，抵消，所以抓不到），$\mathrm{CX}_{8\to 3}$（$X_8$ 传到数据 $X_3$），$\mathrm{CX}_{7\to 1}$（数据 $X_1$）。数据上留下 $X_0X_1X_3$，乘以稳定子 $X_0X_1X_2X_3$ 得 $X_2$，所以 $\mathrm{wt}_S = 1$，满足 $\mathrm{FT}_1$。flag 只需要抓那些会留下"按稳定子算仍然超过 1"的错误。
</details>

## 9. Gate

不看笔记，用 3–5 句话回答：Shor 型 gadget 为什么用 cat 态，为什么还需要 flag？$\mathrm{FT}_1$ 的条件是什么？去掉 flag 的基线电路具体违反在哪里？

<details><summary>参考答案</summary>

用一个辅助比特依次碰 4 个数据比特时，它中途的一个 $X$ 错误会传到多个数据比特；cat 态让每个辅助比特只碰一个数据比特，避免这种传播。但制备 cat 态的 CNOT 树本身会让一个故障翻转两个 cat 比特，所以在 cat 态碰数据之前，用 flag 测两个 cat 比特的 $Z$ 奇偶，把这种坏的 cat 态找出来。$\mathrm{FT}_1$ 要求：任何一个没被 flag 检测到的单故障，在数据上留下的错误按稳定子等价算权重不超过 1。去掉 flag 后，比如 $q_8$ 的制备错误会经 CX $8\to 10$ 扩散，最后留下 $X_1X_3$，乘以 $X_{0123}$ 后仍是权重 2，违反 $\mathrm{FT}_1$；代码一共找到 17 个这样的单故障。
</details>

## 10. 下一课

L02 用 Stim 加上真实比例的噪声做蒙特卡洛模拟：一次只看单故障是枚举（这一课），多个故障同时发生要靠采样。我们会看到，满足 $\mathrm{FT}_1$ 的电路，"flag 没报警但数据坏了"的概率随 $p^2$ 增长，不满足的随 $p$ 增长。
