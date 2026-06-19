import sys

with open('sauron.py', 'r') as f:
    content = f.read()

top_import = """import argparse
try:
    import rinpy
except ImportError:
    rinpy = None
"""

content = content.replace("import argparse", top_import)

with open('sauron.py', 'w') as f:
    f.write(content)
