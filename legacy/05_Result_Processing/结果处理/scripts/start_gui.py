#!/usr/bin/env python3
"""
启动qPCR结果处理器

Usage:
    python start_gui.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.complete_gui import main

if __name__ == '__main__':
    main()

