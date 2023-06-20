import sys

args = sys.argv[2:]
template_file = sys.argv[1]

with open(template_file, 'r') as fh:
    template = fh.read()

print(template.format(args))
