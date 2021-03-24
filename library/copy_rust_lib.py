#!/usr/bin/env python
import shutil
import os
import sys

src_path = os.path.join(
    sys.path[0], "calculate_distances", "target", "release")
dst_path = sys.path[0]
if os.name == "nt":
    src_path = os.path.join(src_path, "libcalculate_distances.dll")
    dst_path = os.path.join(dst_path, "calculate_distances.pyd")
elif os.name == "posix":
    src_path = os.path.join(src_path, "libcalculate_distances.so")
    dst_path = os.path.join(dst_path, "calculate_distances.so")
shutil.copy(src_path, dst_path)
