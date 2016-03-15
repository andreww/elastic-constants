#!/usr/bin/env python
# encoding: utf-8
"""
concat_cij_latex.py

Widen a LaTeX table for an additional Cij calculation

Copyright (c) 2016 Andrew Walker. All rights reserved.

"""

import re

end_row_RE = re.compile(r'\\\\\s*$')
def latex_table_cocat(left_file, right_file, out_file=None):

    # We work with the file in memory, so we can overwrite
    with open(left_file, 'r') as f:
        left_lines = f.read().splitlines()

    with open(right_file, 'r') as f:
        right_lines = f.read().splitlines()

    if out_file is None:
        out_file = left_file

    with open(out_file, 'w') as f:
        for left, right in zip(left_lines, right_lines):
            new_line = end_row_RE.sub(right, left)  
            new_line = new_line + r' \\' + ' \n'
            f.write(new_line)

if __name__ == "__main__":
    import sys
    latex_table_cocat(sys.argv[1], sys.argv[2])
    
