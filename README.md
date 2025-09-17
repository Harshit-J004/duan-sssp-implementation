# Breaking the Sorting Barrier: Fast SSSP Implementation

A C++ implementation of the single-source shortest path algorithm from the paper **"Breaking the Sorting Barrier for Directed Single-Source Shortest Paths"** by Ran Duan, Jiayi Mao, Xiao Mao, Xinkai Shu, Longhui Yin (July 31, 2025), achieving **O(m log^(2/3) n)** time complexity.

## Overview

This implementation breaks the traditional **O((n + m) log n)** sorting barrier of Dijkstra's algorithm using advanced techniques:
- **Recursive separators** for graph decomposition
- **Batched data structures** with range-bounded blocks  
- **Degree transformation** to handle high-degree vertices
- **RadixHeap** for efficient priority queue operations

## Algorithm Features

### Core Components

- **DataStructureD**: Custom data structure with small/large block management for batched vertex processing
- **DegreeTransform**: Converts high-degree vertices into cycles to maintain sparsity
- **RadixHeap**: Specialized heap using bit manipulation for sub-logarithmic operations
- **Recursive Solver**: Multi-level approach with separator-based graph decomposition

### Key Innovations

1. **Separator-based decomposition** using randomized algorithms with early termination
2. **Block-structured batching** with distance-range constraints  
3. **Adaptive parameter selection**: k = ⌊log^(1/3) n⌋, t = ⌊log^(2/3) n⌋
4. **Multi-level recursion** with exponentially growing batch capacities
5. **Degree transformation** creating auxiliary cycles for high-degree vertices

## Performance

- **Time Complexity**: O(m log^(2/3) n) - breaking the sorting barrier
- **Space Complexity**: O(n + m) plus auxiliary structures
- **Practical Performance**: High constant factors, best for large sparse graphs

## Usage

Simply compile and run the provided code:

```bash
g++ -O2 -std=c++17 SSSP_1.cpp -o sssp
./sssp < input.txt
```

## Input Format

```
n m
a1 b1 c1
a2 b2 c2
...
am bm cm
```

Where:
- `n` = number of vertices (1-indexed)
- `m` = number of edges  
- `ai bi ci` = directed edge from `ai` to `bi` with weight `ci`

The algorithm computes shortest distances from vertex 1 to all vertices and outputs them space-separated.

## Implementation Details

### Algorithm Parameters
- `k = max(1, ⌊log^(1/3) n⌋)`: Block size and separator threshold
- `t = max(1, ⌊log^(2/3) n⌋)`: Recursion branching factor
- `lmax = max(1, min(15, ⌈log n / t⌉))`: Maximum recursion depth

### Data Structure Components
- **Small blocks**: Dynamic containers with distance-range bounds
- **Large blocks**: Sorted fixed-capacity overflow containers
- **Batch operations**: Efficient bulk insertions with preprocessing
- **Vertex tracking**: Hash-based duplicate detection and removal

### Graph Preprocessing
- **Degree transformation**: Vertices with degree > 2 replaced by auxiliary cycles
- **Edge preservation**: Original connectivity maintained through zero-weight paths
- **Graph expansion**: New vertex count can be up to O(n + m)

## Algorithmic Flow

1. **Initialization**: Apply degree transform, set parameters
2. **Recursive decomposition**: Find separators using randomized BFS
3. **Batched processing**: Use DataStructureD for efficient vertex management  
4. **Base case**: Standard Dijkstra for small subproblems
5. **Cleanup**: RadixHeap-based final distance computation

## Theoretical Background

The algorithm achieves sub-quadratic performance through:

1. **Separator decomposition**: Limits recursive subproblem sizes
2. **Batched data structures**: Amortizes expensive operations
3. **Degree reduction**: Maintains graph sparsity invariants
4. **Multi-level approach**: Balances recursion depth vs batch size

## Limitations

- **High constant factors**: May be slower than Dijkstra on small/medium graphs
- **Memory overhead**: Complex auxiliary data structures
- **Implementation complexity**: Many interacting components
- **Parameter sensitivity**: Performance depends on graph structure

## References

- Ran Duan, Jiayi Mao, Xiao Mao, Xinkai Shu, Longhui Yin. "Breaking the Sorting Barrier for Directed Single-Source Shortest Paths." July 31, 2025.
- Theoretical analysis of sorting lower bounds for shortest paths
- Advanced techniques in graph algorithms and data structures

## License

MIT License - Feel free to use and modify for research or educational purposes.

## Contributing

This is a research implementation. Contributions for optimizations, bug fixes, or additional features are welcome.
