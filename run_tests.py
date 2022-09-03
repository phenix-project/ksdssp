from __future__ import print_function

from libtbx import easy_run
from libtbx.test_utils import show_diff
import os
import libtbx.load_env

def exercise () :
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf.pdb",
    test=os.path.isfile)
  if pdb_file is None :
    print("skipping")
    return False
  ksdssp_out = easy_run.fully_buffered(command='phenix.ksdssp "%s"' % pdb_file)
  assert not show_diff("\n".join(ksdssp_out.stdout_lines),
"""\
HELIX    1   1 VAL A   15  ASP A   17  5                                   3
HELIX    2   2 ASP A   38  LEU A   47  1                                  10
HELIX    3   3 SER A   58  ARG A   64  1                                   7
HELIX    4   4 ASP A  120  THR A  135  1                                  16
HELIX    5   5 ALA A  139  ALA A  151  1                                  13
HELIX    6   6 ARG A  166  ALA A  177  1                                  12
HELIX    7   7 ARG A  182  GLN A  208  1                                  27
HELIX    8   8 PRO A  217  LYS A  224  1                                   8
HELIX    9   9 ASP A  229  LEU A  232  1                                   4
HELIX   10  10 ALA A  236  TYR A  250  1                                  15
HELIX   11  11 LEU A  253  ASP A  259  1                                   7
HELIX   12  12 GLN A  264  LEU A  274  1                                  11
SHEET    1   A 4 LEU A  27  SER A  30  0
SHEET    2   A 4 VAL A 156  HIS A 159  1  N  VAL A 156   O  PHE A  28
SHEET    3   A 4 ASP A  51  ASP A  54  1  N  ASP A  51   O  LEU A 157
SHEET    4   A 4 ASP A  74  LEU A  77  1  N  ASP A  74   O  VAL A  52""")
  print("OK")

if __name__ == "__main__" :
  exercise()
