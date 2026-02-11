## Flow index：圈空间中的 nowhere-zero 流

本文档说明 **flow index** 脚本 `compute_flow_index.py`，用于计算连通多重图在圈空间中的 nowhere-zero 流的最小最大边范数。

**文件：**

- `compute_flow_index.py`：为每条边分配圈空间系数，通过带约束优化求解 nowhere-zero 流的最小最大范数，并输出 flow value $`r = M + 1`$。

**依赖：** Python 3.8+、`networkx`、`numpy`、`scipy`、`matplotlib`。

---

## 0. 建议的文件放置

将以下文件放在同一目录：

- `compute_flow_index.py`
- `README.md`（英文）
- `README_CN.md`（本文件，中文）

---

## 1. 数学背景

设 $`G`$ 为连通无向多重图，$`n`$ 个顶点、$`m`$ 条边。**圈空间**的维数为 $`k = m - n + 1`$。在选定的一组基下，每条边 $`e`$ 对应一个系数向量 $`c_e \in \mathbb{Z}^k`$。

我们寻找映射 $`X : \mathbb{R}^k \to \mathbb{R}^d`$，使得对每条边 $`e`$ 有

$$
f_e = \sum_{i=1}^{k} c_{e,i} \, X_i \in \mathbb{R}^d
$$

满足：

- **无处为零**：对所有边 $`e`$ 有 $`\|f_e\|_p \ge 1`$；
- **流守恒**：由圈空间构造，每个顶点的净流（Out $`-`$ In）恒为零。

**目标**是在所有可行 $`X`$ 中最小化 $`M = \max_e \|f_e\|_p`$。**Flow value** 定义为 $`r = M + 1`$。

---

## 2. 快速使用

### 2.1 运行 demo

```bash
python compute_flow_index.py
```

会运行内置示例（Petersen 图，$`d=2`$、$`p=2`$、L2 范数），打印图结构、可视化图、执行优化，并输出边向量与流守恒检查。

### 2.2 在代码中使用

```python
import networkx as nx
from compute_flow_index import (
    assign_edge_vectors_multigraph,
    solve_with_constraints,
    visualize_multigraph,
)

# 构建 MultiGraph，每条边带 key=idx, idx=idx
MG = nx.MultiGraph()
MG.add_nodes_from(range(10))
edges = [(i, (i+1)%5) for i in range(5)] + [(i,i+5) for i in range(5)] + [(i+5,(i+2)%5+5) for i in range(5)]
for idx, (u, v) in enumerate(edges):
    MG.add_edge(u, v, key=idx, idx=idx)

# 分配圈空间系数并求解
edge_vectors, tree_edge_keys = assign_edge_vectors_multigraph(MG)
best_M, final_edge_vectors = solve_with_constraints(MG, edge_vectors, d=2, p=2, repeats=20, maxiter=1000)

print(f"Flow value r = {best_M + 1:.6f}")
```

---

## 3. 参数说明

在 `main_example()` 中可调整：

- **`d`**：流向量的维数（默认 2）。
- **`p`**：范数 $`\|\cdot\|_p`$ 的阶（默认 2，即欧氏范数）。
- **`repeats`**：优化器随机重启次数（默认 40）。
- **`maxiter`**：每次 SLSQP 的最大迭代次数（默认 5000）。

也可通过取消注释示例图定义或自定义 `edges` 来更换图（轮图、圈图、偶极子图等）。

---

## 4. API

### 4.1 assign_edge_vectors_multigraph

**`assign_edge_vectors_multigraph(MG: nx.MultiGraph) -> (edge_vectors, tree_edge_keys)`**

- **`MG`**：连通 `networkx.MultiGraph`，无自环。每条边需有 `key`，可选 `idx`。
- **返回**：`edge_vectors` 为字典，`(u, v, key) ->` 长度为 $`k = m - n + 1`$ 的整数向量（圈空间系数）；`tree_edge_keys` 为生成树边列表。

### 4.2 solve_with_constraints

**`solve_with_constraints(MG, edge_vectors, d=2, p=2, repeats=20, maxiter=1000)`**

- **`MG`**、**`edge_vectors`**：来自 `assign_edge_vectors_multigraph`。
- **`d`**：流向量维数。
- **`p`**：范数阶（如 2 表示 L2）。
- **`repeats`**：随机重启次数。
- **`maxiter`**：每次 SLSQP 的迭代上限。
- **返回**：`(best_M, final_edge_vectors)`，其中 `best_M` 为最小最大边范数，`final_edge_vectors` 为边到 $`\mathbb{R}^d`$ 流向量的映射。

### 4.3 visualize_multigraph

**`visualize_multigraph(G, tree_edges, edge_vectors)`**

绘制图结构（Petersen 图使用自定义布局，其他图使用 Kamada–Kawai）。

### 4.4 辅助函数

- **`print_graph_structure(G)`** — 打印顶点、边及圈空间维数。
- **`print_edge_vectors(MG, final_edge_vectors)`** — 按 idx 打印每条边的流向量及范数。
- **`print_flow_conservation(MG, final_edge_vectors, p=2)`** — 检查各顶点的流守恒（Out $`-`$ In $`= 0`$）。

---

## 5. 示例图

脚本中提供了注释掉的示例：

- **Petersen 图**（默认）：10 点、15 边，$`k=6`$。
- **轮图** $`W_k`$：中心点 + 外圈。
- **圈图** $`C_k`$：单圈。
- **偶极子图** $`Mu_k`$：两点间 $`k`$ 条重边。

---

## 6. 性能说明

- 优化使用 SLSQP，约束为 $`\|f_e\|_p \ge 1`$，目标为最小化 $`\max_e \|f_e\|_p`$。多次随机重启有助于找到更优解。
- 主要耗时在优化器；圈空间分配与可视化为 $`O(m)`$、$`O(m^2)`$ 量级。
- 图较大时可适当增大 `repeats` 和 `maxiter`。

