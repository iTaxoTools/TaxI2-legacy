#!/usr/bin/env python3

import os
import sys
import numpy as np
import tempfile
import tkinter as tk
import multiprocessing

from library.gui import TaxiGUI

resource_path = getattr(sys, '_MEIPASS', sys.path[0])

def gui_main() -> None:
    root = tk.Tk()

    def close_window():
        root.destroy()
        root.quit()

    root.title("TaxI2")
    if os.name == "nt":
        root.wm_iconbitmap(os.path.join(
            resource_path, 'data', 'TaxI2.ico'))

    root.protocol("WM_DELETE_WINDOW", close_window)
    root.rowconfigure(0, weight=1)
    root.columnconfigure(0, weight=1)

    preview_dir = tempfile.mkdtemp()
    TaxiGUI(root, preview_dir=preview_dir)

    root.mainloop()
    root.quit()


def main() -> None:
    np.seterr(all='ignore')
    gui_main()


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
