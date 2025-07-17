# Shortest Path Problem on Affine Subspaces

This directory contains the source code for the Shortest Path Problem (SPP) on Affine Subspaces (AS) described [here](../../report/SPP-affine_subspaces.pdf).

In particular, the [`Graph.jl`](./Graph.jl) (resp. [`AffineSubspace.jl`](./AffineSubspace.jl)) file contains the structure of a directed graph (resp. of an AS). 

The [`Problem.jl`](./Problem.jl) file contains the structure of our particular problem. The main component is a set of AS, a source and a target.

An implementation of the classical SPP is available in the file [`SPP.jl`](./SPP.jl).

Finally, the files [`SPP_intersections_vertices.jl`](./SPP_intersections_vertices.jl) and [`SPP_intersections_edges.jl`](./SPP_intersections_edges.jl) contains respectively the implementations of the SPP on AS, considering the intersections either as the vertices or as the edges of the graph, as described in the section 2 of the [report](../../report/SPP-affine_subspaces.pdf)