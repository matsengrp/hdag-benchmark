## A Lark-based newick and BEAST nexus parser

This parser is designed to properly load a variety of newick formats, as well
as formats for extended newick (`[&NHX...]` and `[&attr]`) format node
annotations, including annotations that contain commas and nested lists, which
can confuse packages like ete3 and Dendropy.
The parsing target is lists of `ete3.Tree` objects.

## Requirements:
* Lark
* ete3

## Usage:
Functions for parsing files can be found in `tree_transformer.py`, for example,
* `parse_newick`
* `parse_newick_file`
* `parse_nexus_file`
* `iter_nexus_trees`


All node attributes are added as `ete3.TreeNode` attributes, preserved by tree
copy, with the same attribute name that appears in the newick annotation.
