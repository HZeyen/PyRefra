# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 13:34:18 2024

@author: Hermann
"""

from pathlib import Path
from runpy import run_path

pkg_dir = Path(__file__).resolve().parent


def execute_script():
    script_pth = pkg_dir / "pyrefra_script.py"
    run_path(str(script_pth), run_name="__main__")
