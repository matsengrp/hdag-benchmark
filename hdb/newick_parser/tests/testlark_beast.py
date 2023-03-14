from lark import Lark, Transformer
from tree_transformer import TreelistTransformer, NexusTransformer


with open('newick.lark', 'r') as fh:
    beast_translate_parser = Lark(fh, parser='lalr', start='beast_file', transformer=NexusTransformer())

# with open('newick.lark.ebnf', 'r') as fh:
#     newick_parser_inline = Lark(fh, parser='lalr', transformer=TreelistTransformer())
#     # newick_parser = Lark(fh, debug=True)



with open("clade_13.GTR.trees.nexus", 'r') as fh:
    for line in fh:
        if "Begin trees;" in line:
            data = fh.read()

trees = beast_translate_parser.parse(data)
print(trees)
