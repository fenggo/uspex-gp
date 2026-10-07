#!/usr/bin/env bash
# PoC launcher for USPEX-GP RL bandit (mirrors the official octave invocation)
OCTROOT=/home/feng/mathlib/octave-10.3.0
export LD_LIBRARY_PATH=$OCTROOT/libmex/.libs:$OCTROOT/libinterp/.libs:$OCTROOT/liboctave/.libs:$OCTROOT/src/.libs

P=$OCTROOT/scripts
PATHS=$P
for d in help ode time strings general specfun plot plot/appearance plot/util plot/draw statistics io signal miscellaneous geometry polynomials prefs elfun testfun deprecated image legacy gui audio linear-algebra java java/org java/org/octave path set optimization special-matrix startup
do
  PATHS="$PATHS:$P/$d"
done
PATHS="$PATHS:$OCTROOT/libinterp:$OCTROOT/libinterp/octave-value:$OCTROOT/libinterp/operators:$OCTROOT/libinterp/template-inst:$OCTROOT/libinterp/dldfcn:$OCTROOT/libinterp/corefcn:$OCTROOT/libinterp/parse-tree:$OCTROOT/examples/data:$OCTROOT/libgui/graphics"

cd /home/feng/uspex-gp
exec $OCTROOT/src/.libs/octave-cli --no-init-path --path="$PATHS" \
  --image-path=.:$P/image \
  --doc-cache-file=$OCTROOT/doc/interpreter/doc-cache \
  --built-in-docstrings-file=$OCTROOT/libinterp/DOCSTRINGS \
  --texi-macros-file=$OCTROOT/doc/interpreter/macros.texi \
  --info-file=$OCTROOT/doc/interpreter/octave.info \
  --quiet USPEX.m
