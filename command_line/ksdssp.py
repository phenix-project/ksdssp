import libtbx.load_env
from libtbx.utils import Usage
import os, sys, subprocess

def run(args):
  exe = libtbx.env.under_build("ksdssp/exe/ksdssp")
  if (os.name == "nt"):
    exe += ".exe"
  if (not os.path.isfile(exe)):
    from libtbx.str_utils import show_string
    from libtbx.utils import Sorry
    raise Sorry("Missing phenix.ksdssp executable: %s" % show_string(exe))
  if len(args) == 0 or not os.path.isfile(args[0]) :
    raise Usage("phenix.ksdssp model.pdb")
  subprocess.call([exe] + args)

if (__name__ == "__main__"):
  run(sys.argv[1:])
