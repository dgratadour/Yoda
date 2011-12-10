#include "optimpack-test.i"
require, "utils.i";

func do_test(fatol=, frtol=, sftol=, sgtol=, sxtol=, output=)
{
  keys = ["fatol", "frtol", "sftol", "sgtol", "sxtol"];
  code = "for (i=1;i<=18;++i) optim_test_um, prob=i, method=0, verb=1";
  for (i=1;i<=numberof(keys);++i) {
    key = keys(i);
    val = symbol_def(key);
    if (! is_void(val)) code += ", "+key+"="+val;
  }
  if (! is_void(output)) {
    code += ", output=\""+output+"\";";
    write, open(output, "w"), linesize=1000, format="# %s\n", code;
  } else {
    code += ";";
  }
  eval, code, debug=0;
}

do_test, output="test1",
  fatol="0.0", frtol="1e-12", sftol="1e-3", sgtol="0.9", sxtol="0.1";
do_test, output="test2",
  fatol="0.0", frtol="1e-12", sftol="1e-3", sgtol="0.1", sxtol="0.1";
do_test, output="test3",
  fatol="0.0", frtol="1e-12", sftol="1e-2", sgtol="0.1", sxtol="0.1";
do_test, output="test4",
  fatol="0.0", frtol="1e-12", sftol="5e-2", sgtol="0.1", sxtol="0.1";
