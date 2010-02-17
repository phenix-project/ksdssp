# LIBTBX_SET_DISPATCHER_NAME phenix.ksdssp

import libtbx.load_env
import os, sys, subprocess

def run(args):
  exe = libtbx.env.under_build("ksdssp/exe/ksdssp")
  if (os.name == "nt"):
    exe += ".exe"
  if (not os.path.isfile(exe)):
    from libtbx.str_utils import show_string
    from libtbx.utils import Sorry
    raise Sorry("Missing phenix.ksdssp executable: %s" % show_string(exe))
  subprocess.call([exe] + args)

if (__name__ == "__main__"):
  run(sys.argv[1:])
