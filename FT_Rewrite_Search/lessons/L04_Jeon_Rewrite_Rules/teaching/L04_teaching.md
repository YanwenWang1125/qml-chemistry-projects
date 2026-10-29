# L04 Jeon 的四条改写规则，以及 Bell reduction

材料：Jeon et al.（arXiv:2610.09749）第 II C 节、第 III 节（式 (5)–(12)，PDF 第 3–5 页）、第 IV 节（搜索算法）、附录 B（PDF 第 13 页）。渲染图：[figures/jeon2026_page-04.png](figures/jeon2026_page-04.png)（四条规则），[figures/jeon2026_page-13.png](figures/jeon2026_page-13.png)（附录 B）。附录 A 里的 ZX 证明先当黑盒。

代码：[code/L04_bell_reduction.py](code/L04_bell_reduction.py)，输出 [code/L04_output.txt](code/L04_output.txt)。共用工具在 L01 的 [ftlib.py](../../L01_FT_and_Shor_Gadget/teaching/code/ftlib.py)。

## 1. 这一课讲什么

L03 讲了"保持故障等价"这个要求，这一课讲 Jeon 实际使用的四条具体规则。其中三条（受限交换、basis swap、target swap）不减少资源，只调整电路的形状；真正省资源的是 Bell reduction，每用一次删掉一个辅助比特和一个 CNOT。我们先在最小的例子上手推并用代码验证 Bell reduction，再看 Jeon 怎样连用两次，把基线电路（5 个辅助比特、9 个 CNOT）变成优化后的电路（3 个、7 个）。学完这一课，阶段 0 的验收就可以做了。

## 2. 符号表

| 符号 | 含义 | 容易混的地方 |
|---|---|---|
| $\hat{=}$ | 论文里表示"故障等价" | 不只是功能相同 |
| $\lvert \Phi^+ \rangle$ | Bell 态 $\frac{1}{\sqrt{2}}(\lvert 00 \rangle + \lvert 11 \rangle)$ | |
| $q_1, q_2$（规则 (9) 里） | Bell 对的两个比特：$q_1$ 从 $\lvert + \rangle$ 出发做控制，$q_2$ 从 $\lvert 0 \rangle$ 出发做目标 | 和我们 gadget 的 $q_7$ 等编号无关 |
| $U, V$ | Bell 对建立之后，$q_1$ 和 $q_2$ 各自经历的操作（可以涉及其他比特） | |
| $V^{\mathsf{T}}$ | $V$ 的转置；对这一课的例子，就是把 $q_2$ 上的操作顺序倒过来 | |
| $m, q$ | 被删掉的测量结果 $m$，以及其余测量结果的奇偶 $q$ | |
| Bell reduction $(a, b)$ | 附录 B 的写法：去掉 $b$，由 $a$ 接管 | |

## 3. 概念：先白话，再定义

四条规则的图见 [figures/jeon2026_page-04.png](figures/jeon2026_page-04.png)。

**规则 1：受限交换（式 (5)、(6)）。** 两个 CNOT 共用一个比特：
- 式 (5)：共用的比特是两个 CNOT 的**控制**，而且刚制备成 $\lvert + \rangle$；
- 式 (6)：共用的比特是两个 CNOT 的**目标**，而且刚制备成 $\lvert 0 \rangle$。

这时两个 CNOT 可以交换顺序。另外两个端点必须是不同的比特。

无噪声时，共用控制（或共用目标）的两个 CNOT 本来就对易，所以这条规则的意义全在噪声上：交换顺序会改变故障的传播路径，规则保证在这种初始化条件下交换仍然故障等价。论文明确说，共用比特是任意输入态时不能用。

**规则 2：basis swap（式 (7)）。** 一个 CNOT 从 $\lvert + \rangle$ 控制到 $\lvert 0 \rangle$，制备出 Bell 态 $\lvert \Phi^+ \rangle$。Bell 态对交换两个比特是对称的，所以可以把两根线的角色对调：原来的目标线变成从 $\lvert + \rangle$ 出发的控制。这不是插入一个物理的 SWAP 门，而且只对这种初始化成立。

**规则 3：target swap（式 (8)）。** 三个比特分别初始化为 $\lvert + \rangle, \lvert + \rangle, \lvert 0 \rangle$。先做 CNOT $q_2 \to q_3$，把 $q_2, q_3$ 做成 Bell 对；之后从 $q_1$ 出发的 CNOT，打在 $q_2$ 上和打在 $q_3$ 上是故障等价的。白话：Bell 对的两个比特"互为替身"，对其中一个做 $X$ 等于对另一个做 $X$，即 $X_2 \lvert \Phi^+ \rangle = X_3 \lvert \Phi^+ \rangle$。

**规则 4：Bell reduction（式 (9)）。** $q_1$（$\lvert + \rangle$）和 $q_2$（$\lvert 0 \rangle$）用一个 CNOT 做成 Bell 对，之后 $q_1$ 经历 $U$，$q_2$ 经历 $V$ 后被测量。如果满足两个条件：

1. Bell 对建立后，两个比特各干各的，直到测量；没有任何门同时作用在 $q_1$ 和 $q_2$ 上，即联合演化可以写成 $U \otimes V$；
2. $q_2$ 的测量结果 $m$ 不单独使用（不当 flag、不控制后续操作），只参与一个奇偶 $p = m \oplus q$。

那么整个结构可以换成**一根线**：从 $\lvert 0 \rangle$ 出发，先做 $V^{\mathsf{T}}$，再做 $U$。省掉一个比特和那个建 Bell 对的 CNOT。

为什么条件 2 必要：只保留奇偶时，记录 $(m, q)$ 和 $(0, m \oplus q)$ 在操作上没有区别，可以把 $m$ 固定为 0，它的值并进剩下的奇偶里（论文式 (10)）。这只是重新标记经典记录，不是后选择。

**另外三条规则的作用**（论文第 IV 节）：它们本身不降成本，只用来把电路"摆"成 Bell reduction 能用的样子。Jeon 的搜索把"若干条使能规则 + 一次 Bell reduction"打包成一步，每一步恰好删掉一个辅助比特和一个 CNOT。为了控制分支数，只有当使能规则能**立刻**暴露出一个 Bell 对时才使用，不做多步的投机串联。这正是 L4 想用 ML 放开的限制。

## 4. 一步步推导：在最小的例子上做 Bell reduction

目标：测量两个数据比特上的 $X_0X_1$。

**电路 A（两个辅助比特）：**

$$
q_2: \lvert + \rangle,\quad q_3: \lvert 0 \rangle,\quad \mathrm{CX}_{2\to 3},\quad \mathrm{CX}_{2\to 0},\quad \mathrm{CX}_{3\to 1},\quad q_2, q_3 \text{ 在 } X \text{ 基测量，只用 } m_2 \oplus m_3
$$

**第 1 步：检查条件。** $\mathrm{CX}_{2\to 3}$ 之后，$q_2$ 只碰 $d_0$，$q_3$ 只碰 $d_1$，没有同时作用在 $q_2, q_3$ 上的门：条件 1 成立。只用 $m_2 \oplus m_3$：条件 2 成立。

**第 2 步：对号入座。** 规则 (9) 里的 $q_1$ 是我们的 $q_2$，规则里的 $q_2$ 是我们的 $q_3$。

$$
U = \left(\mathrm{CX}_{2\to 0},\ \text{然后 } H,\ \text{测量}\right), \qquad V = \left(\mathrm{CX}_{3\to 1},\ \text{然后 } H\right)
$$

**第 3 步：写出 $V^{\mathsf{T}}$。** 转置把顺序倒过来：$V$ 是"先 CX 后 H"，$V^{\mathsf{T}}$ 就是"先 H 后 CX"。$H$ 是对称矩阵，$H^{\mathsf{T}} = H$；CNOT 对控制比特是对角的，对控制线做转置不改变它。

**第 4 步：拼成一根线。** 从 $\lvert 0 \rangle$ 出发，先 $V^{\mathsf{T}}$ 再 $U$：

$$
\lvert 0 \rangle \xrightarrow{H} \lvert + \rangle \xrightarrow{\mathrm{CX}\to d_1} \xrightarrow{\mathrm{CX}\to d_0} \xrightarrow{H} \text{测量}
$$

**第 5 步：用受限交换整理。** 两个 CNOT 共用的控制比特刚制备成 $\lvert + \rangle$，满足式 (5)，可以交换成先 $d_0$ 后 $d_1$。得到

**电路 B（一个辅助比特）：**

$$
q_2: \lvert + \rangle,\quad \mathrm{CX}_{2\to 0},\quad \mathrm{CX}_{2\to 1},\quad q_2 \text{ 在 } X \text{ 基测量}
$$

少了一个辅助比特（$q_3$）和一个 CNOT（$\mathrm{CX}_{2\to 3}$）。

**独立验证：两个电路测的确实是同一个算符**（标准的 Heisenberg 绘景推导，教材没写）。把"最后测的算符"倒着推过每个门。CNOT 满足 $\mathrm{CX}\, X_c\, \mathrm{CX} = X_c X_t$，所以"门之后测 $X_c X_t$"等于"门之前测 $X_c$"，而"门之后测 $X_t$"等于"门之前测 $X_t$"。

电路 A，最后测 $X_2X_3$：

$$
\begin{aligned}
X_2 X_3 &\xleftarrow{\ \mathrm{CX}_{3\to 1}\ } X_2 X_3 X_1 \\
&\xleftarrow{\ \mathrm{CX}_{2\to 0}\ } X_2 X_0 X_3 X_1 \\
&\xleftarrow{\ \mathrm{CX}_{2\to 3}\ } X_2 X_0 X_1
\end{aligned}
$$

最后一步：$X_2X_3$ 在 $\mathrm{CX}_{2\to 3}$ 之后，等于之前的 $X_2$。$q_2$ 初始为 $\lvert + \rangle$，是 $X_2$ 的 $+1$ 本征态，所以测量结果就是 $X_0X_1$ 的值。

电路 B，最后测 $X_2$：

$$
X_2 \xleftarrow{\ \mathrm{CX}_{2\to 1}\ } X_2 X_1 \xleftarrow{\ \mathrm{CX}_{2\to 0}\ } X_2 X_0 X_1
$$

同样得到 $X_0X_1$。两个电路无噪声时测的是同一个量。

## 5. 完整数值例子：资源和故障位置

| | 电路 A | 电路 B |
|---|---|---|
| 辅助比特 | 2 | 1 |
| CNOT | 3 | 2 |
| CNOT 深度 | 2 | 2 |
| 单故障位置-类型组合（L01 的故障模型） | 58 | 38 |

手数 B 的 38：2 个 CNOT 各 15 种两比特 Pauli，共 30；1 次制备 1 种；2 个 H 各 3 种，共 6；1 次测量 1 种。$30 + 1 + 6 + 1 = 38$。A：3 个 CNOT 共 45，2 次制备 2，3 个 H 共 9（$q_2$ 开头一个、末尾两个），2 次测量 2，合计 58。

代码输出也是 58 和 38（第 6 节）。

**再到 Jeon 的 Shor gadget。** 比较附录 B 的起点和终点（[figures/jeon2026_page-13.png](figures/jeon2026_page-13.png)）：

| 原来 | 碰的数据比特 | 优化后 |
|---|---|---|
| $q_7$ | $d_0$ | $q_7$，碰 $d_2, d_0, d_1$ |
| $q_8$ | $d_1$ | 被删除，由 $q_7$ 接管（reduction $(7,8)$） |
| $q_9$ | $d_2$ | 被删除，由 $q_7$ 接管（reduction $(7,9)$） |
| $q_{10}$ | $d_3$ | 改名为 $q_8$，仍碰 $d_3$ |
| $q_{11}$（flag） | | 改名为 $q_9$ |

原来 flag 由 $q_7$ 和 $q_8$ 各打一个 CNOT；$q_8$ 并进 $q_7$ 以后，变成 $q_7$ 自己打两个 CNOT，这就是 L01 第 5 节看到的"$q_7$ 打两次 flag"。两次 reduction 各省 1 个辅助比特、1 个 CNOT：$5 \to 3$，$9 \to 7$。

这张对照表是比较起点和终点电路得出的；中间每一步用了哪些使能规则，请看附录 B 的中间那张图。

## 6. 代码结果

```text
Part 1: Bell reduction on the smallest example
  A: Bell pair of ancillas   {'ancillas': 2, 'cnots': 3, 'cnot_depth': 2}, single-fault locations: 58, noiseless syndrome {'X0X1=+1': [0], 'X0X1=-1': [1]}
  B: one ancilla             {'ancillas': 1, 'cnots': 2, 'cnot_depth': 2}, single-fault locations: 38, noiseless syndrome {'X0X1=+1': [0], 'X0X1=-1': [1]}
  effects only in A: none
  effects only in B: none
  passes the single-fault check: True

  Monte Carlo at p = 0.01, 2000000 shots, data in |0>^2:
    A: Bell pair of ancillas   P(data X error mod X0X1) = 0.02567
    B: one ancilla             P(data X error mod X0X1) = 0.01065

Part 2: Jeon et al.'s Shor gadget before and after two Bell reductions
  baseline (Fig. 1a)   {'ancillas': 5, 'cnots': 9, 'cnot_depth': 5}, FT_1: True, undetected single faults: 117
  reduced (Fig. 1b)    {'ancillas': 3, 'cnots': 7, 'cnot_depth': 6}, FT_1: True, undetected single faults: 86
```

每项结果说明了什么：

| 结果 | 说明了什么 |
|---|---|
| A、B 的无噪声 syndrome 都正确 | 和第 4 节 Heisenberg 推导一致：两个电路测的都是 $X_0X_1$ |
| 单故障检查通过 | A 和 B 在 L03 的简化检查下故障等价，和 Bell reduction 规则的结论一致 |
| A 的数据错误率 0.0257，B 为 0.0107 | **同一个故障等价类，错误率差了约 2.4 倍**：B 少了一个比特、一个 CNOT，出错位置少 |
| Part 2 | 两次 Bell reduction 后资源从 5 / 9 / 5 变成 3 / 7 / 6，仍满足 $\mathrm{FT}_1$（L01 已逐个检查） |

**边界：**
- 单故障检查只是故障等价的必要条件的数值验证（见 L03 第 4 节），不是证明。规则本身的证明在论文附录 A（ZX 演算），这里当黑盒。
- 数据错误率只统计了数据上的 $X$ 错误（模 $X_0X_1$），没有空闲噪声。
- Shor gadget 的中间推导步骤没有在代码里重现，只验证了起点和终点。

## 7. 常见误区

| 错误想法 | 为什么容易这样想 | 正确的理解 |
|---|---|---|
| 受限交换没有意义，共用控制的 CNOT 本来就对易 | 只考虑了无噪声的情况 | 无噪声时确实对易；规则保证的是有噪声时交换后仍故障等价，所以要求"刚制备"这个条件 |
| basis swap 就是插一个 SWAP 门 | 图上两根线互换了 | 不加任何门，只是利用 Bell 态的对称性重新解释两根线的角色；只对 $\lvert + \rangle$、$\lvert 0 \rangle$ 这种初始化成立 |
| Bell reduction 是后选择 $m = 0$ | 论文说"把 $m$ 固定为 0" | 不丢弃任何 shot；$m$ 的值被并进剩下的奇偶，只是重新标记经典记录 |
| 只要有一个 Bell 对就能用 Bell reduction | 只记住了"Bell 对" | 还需要两个条件：之后各自独立演化；被删的测量结果只参与奇偶 |
| 每用一次规则都会省资源 | 规则统称"改写" | 只有 Bell reduction 省资源；另外三条只调整形状，而且可能暂时让电路"不变好" |
| 规则保持故障等价，所以改写后错误率不变 | 混淆了故障等价和错误率（L03） | A 和 B 故障等价，错误率差约 2.4 倍 |

## 8. 自测

**1.**（数值）基线 Shor gadget 一整轮测 6 个稳定子，每个 gadget 都做两次 Bell reduction。一整轮的辅助比特制备次数和 CNOT 数各是多少？和论文对上了吗？

<details><summary>答案</summary>

每个 gadget 从 5 个辅助比特、9 个 CNOT 降到 3 个、7 个。6 个 gadget：辅助比特制备 $6 \times 3 = 18$ 次（原来 30），CNOT $6 \times 7 = 42$ 个（原来 54）。和论文摘要里的 $30 \to 18$、$54 \to 42$ 一致。
</details>

**2.** 在电路 A 里，如果最后的协议要单独用 $m_3$（比如把它当 flag），还能做 Bell reduction 吗？为什么？

<details><summary>答案</summary>

不能。条件 2 要求被删除的测量结果只参与奇偶。如果单独使用 $m_3$，记录 $(m_3, q)$ 和 $(0, m_3 \oplus q)$ 就不再等价，删掉 $q_3$ 会丢失信息。
</details>

**3.**（边界情况）把电路 A 改成：建好 Bell 对后，再加一个 $\mathrm{CX}_{2\to 3}$，然后才去碰数据。能直接对原来那个 Bell 对做 Bell reduction 吗？

<details><summary>答案</summary>

不能直接做。新加的 $\mathrm{CX}_{2\to 3}$ 同时作用在 $q_2$ 和 $q_3$ 上，违反条件 1（联合演化不能写成 $U \otimes V$）。（两个连续的 $\mathrm{CX}_{2\to 3}$ 会抵消，先把它们消掉，再做 Bell reduction 是另一回事；那需要先确认"消掉"这一步本身是故障等价的。）
</details>

**4.**（新情形）Jeon 的搜索只在使能规则能"立刻"暴露 Bell 对时才用它们。举一个理由说明，为什么允许多步串联可能找到更好的电路，以及为什么这会让搜索变难。

<details><summary>答案</summary>

有些 Bell 对要连续用两三条使能规则才能暴露出来；只看一步的搜索永远到不了这些电路，可能错过资源更少或错误率更低的结果（论文自己承认这个限制）。但使能规则本身不降成本，允许串联后，每一步可选的规则数会乘起来，搜索树分支爆炸；而且中间步骤看起来"没变好"，需要往前看几步才能判断值不值得走。这正是价值网络和搜索策略要解决的问题。
</details>

**5.** 规则 (9) 推出的单线电路是"先碰 $d_1$ 再碰 $d_0$"，我们的电路 B 是"先碰 $d_0$ 再碰 $d_1$"。为什么可以这样调换？

<details><summary>答案</summary>

两个 CNOT 共用的控制比特刚制备成 $\lvert + \rangle$，另外两个端点 $d_0, d_1$ 不同，满足受限交换式 (5)，交换后仍故障等价。
</details>

## 9. Gate（阶段 0 验收）

不看笔记完成：

1. 在纸上对电路 A 做一次 Bell reduction，写出每一步（检查两个条件、写出 $U$ 和 $V$、得到单线电路、用受限交换整理），最后得到电路 B。
2. 用 3–5 句话说明：Bell reduction 的两个条件是什么、各自为什么必要；为什么 A 和 B 故障等价，错误率却不同。

<details><summary>参考答案</summary>

1. 见第 4 节第 1–5 步。
2. 条件一：Bell 对建立后两个比特各自独立演化，没有门同时作用在两者上，这样整个过程才能写成 $U \otimes V$，利用 $(U \otimes V)\lvert \Phi^+ \rangle = (U V^{\mathsf{T}} \otimes I)\lvert \Phi^+ \rangle$ 把 $V$ 挪到另一根线上。条件二：被删掉的测量结果只参与奇偶，这样把它固定为 0、把它的值并进剩下的奇偶，在操作上没有区别，不丢信息，也不是后选择。A 和 B 的每种"悄悄出错"的方式都能一一对应，所需故障数不变，所以故障等价；但 B 少了一个比特和一个 CNOT，出错的位置更少，所以错误率更低（模拟中约 0.011 对 0.026）。
</details>

## 10. 下一步

阶段 0 到这里结束。阶段 1（11-09 起）要把这四条规则写成代码，在 [[7,1,3]] 上复现 Jeon 的两个案例：Shor gadget 每个只有 7 个可达电路，Steane gadget 有 391 个候选。届时会用 paritea 做独立的故障等价检查，用这一课的 Bell reduction 例子做单元测试。
