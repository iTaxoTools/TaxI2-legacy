#!/usr/bin/env python3

import sys
from library.gui_utils import *
from library.programstate import *
import tkinter.messagebox as tkmessagebox
import numpy as np

from library.plot_taxi import Plot


def gui_main(debug: bool) -> None:
    root = tk.Tk()
    root.rowconfigure(0, weight=1)
    root.columnconfigure(0, weight=1)
    root.columnconfigure(2, weight=1)

    programstate = ProgramState(root)

    input_file_chooser = FileChooser(root, label="Input File", mode="open")
    output_dir_chooser = FileChooser(
        root, label="Output Directory", mode="dir")

    aligned_chk = ttk.Checkbutton(
        root, variable=programstate.already_aligned, text="Already aligned")

    distances_frm = ttk.Frame(root, padding=5, relief='sunken')

    for kind in range(NDISTANCES):
        checkbutton = ttk.Checkbutton(
            distances_frm, variable=programstate.distance_options[kind], text=distances_names[kind])
        checkbutton.grid(row=kind, column=0, sticky="w")

    format_chooser = LabeledCombobox(root, label="Input file format", values=list(
        ProgramState.formats.keys()), readonly=True, var=programstate.input_format_name)

    print_alignments_chk = ttk.Checkbutton(
        root, variable=programstate.print_alignments, text="Print alignments")

    def process() -> None:
        with display_errors_and_warnings(debug):
            input_file = input_file_chooser.file_var.get()
            output_dir = output_dir_chooser.file_var.get()
            programstate.process(
                input_file, output_dir)
            plot_input = os.path.join(
                output_dir_chooser.file_var.get(), ProgramState.SUMMARY_STATISTICS_NAME)
            distance_name = [distance for distance, is_chosen in zip(
                distances_short_names, programstate.distance_options) if is_chosen.get()]
            print(plot_input, output_dir, distance_name)
            Plot(plot_input, output_dir, distance_name)
            tkmessagebox.showinfo("Done", "Calculation complete.")

    process_btn = ttk.Button(root, text="Calculate distances", command=process)

    input_file_chooser.grid(row=0, column=0)
    output_dir_chooser.grid(row=0, column=2)

    format_chooser.grid(row=1, column=0)
    distances_frm.grid(row=2, column=1)

    process_btn.grid(row=3, column=1)
    aligned_chk.grid(row=3, column=2)
    print_alignments_chk.grid(row=4, column=2)

    root.mainloop()


def main() -> None:
    if "debug" in sys.argv:
        np.seterr(all='ignore')
        gui_main(debug=True)
    else:
        np.seterr(all='ignore')
        gui_main(debug=False)


if __name__ == "__main__":
    main()
