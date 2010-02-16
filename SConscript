import libtbx
Import("env_base", "env_etc")

ksdssp_scons_env = env_base.Clone(
  LIBS=env_etc.libm)
Export("ksdssp_scons_env")

SConscript("libpdb++/SConscript")
SConscript("ksdssp_src/SConscript")
