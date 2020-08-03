# DESCRIPTION
Belief-propagation guided decimation (BPD) as a solver for the maximum feedback vertex set problem. This program is applicable on a single graph instance. The imput graph file has the following form:

> N   M       % first row specifies N (vertex number) and M (edge number)

> i_1  j_1                                            % undirected edge (i_1, j_1)

> i_2  j_2                                            % undirected edge (i_2, j_2)

> .    .

> .    .

> .    .

> i_M  j_M                                          % undirected edge (i_M, j_M)


The program reads only the first EdgeNumber edges from the imput graph, with EdgeNumber being explicitly specified (EdgeNumber <= M, of course).

Each vertex has two states (unoccupied, 0; occupied, 1). An unoccupied vertex belongs to the constructed feedback vertex set S, while an occupied vertex does not belong to S.

Initially all the vertices are occupied and active. At each round of decimation, some of the remaining active vertices are fixed to be un-occupied sequentially. After a vertex is fixed to be un-occupied (and it becomes inactive), the remaining active subgraph is simplified (during which more vertices become inactive). During the simplication process, all the leaf vertices in the active subgraph are removed recursively (and they are fixed to be occupied).

# REFERENCES
For details about the replica-symmetric mean field theory and the spin glass model on which this algorithm was based, please consult the following references:

> 1. Hai-Jun Zhou, "Spin Glasses and Message Passing" (Science Press, Beijing,2015), chapter 7, pages 218--238.
> 2. Hai-Jun Zhou, "Spin glass approach to the feedback vertex set problem",European Physical Journal B 86: 455 (2013).

# COMPILE
To generate the executive file, you can simply compile as

`c++ -O3 -o bpdfvs.exe FVSbpdV01.cpp`

!!!please make sure some of the key parameters, such as input graph name, number
of edges, number of vertices, output file names, as appropriately specified in 
the FVSbpdV01.cpp file.

# USE
After you successfully compiled the code, you can then run the algorithm as

`bpdfvs.exe`

Good luck!

# HISTORY

- 06.09.2015--08.09.2015: FVSbpdV01.cpp revision. Make the algorithm works even if some of the vertices have extremely many attached edges.
- 06.09.2015: FVSbpdV01.cpp (copied from FVSbpdv00.cpp)
- 29.06.2013--05.07.2013: FVSbpdv00.cpp

# PROGRAMMER

- Hai-Jun Zhou
- Institute of Theoretical Physics, Chinese Academy of Sciences
- Zhong-Guan-Cun East Road 55, Beijing 100190
- email: zhouhj@itp.ac.cn
- webpage: http://home.itp.ac.cn/~zhouhj
