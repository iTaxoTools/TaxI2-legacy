import os
import sys
import io
import shutil
import logging
import warnings
from typing import Any, Callable, Iterator

import tkinter as tk
import tkinter.ttk as ttk
import tkinter.font as tkfont
import tkinter.messagebox as tkmessagebox
import tkinter.filedialog as tkfiledialog

from library.programstate import *
from library.gui_utils import LabeledCombobox


class TaxiGUI(ttk.Frame):

    def __init__(self, *args: Any, preview_dir, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)

        self.images = {}
        self.images["txt_icon"] = tk.PhotoImage(
            file=os.path.join(sys.path[0], "data/file-text.png"))
        self.images["graph_icon"] = tk.PhotoImage(
            file=os.path.join(sys.path[0], "data/file-graph.png"))
        self.images["log_icon"] = tk.PhotoImage(
            file=os.path.join(sys.path[0], "data/file-log.png"))
        self.preview_dir = preview_dir
        self.programstate = ProgramState(self)

        self.create_top_frame()
        self.create_parameters_frame()
        self.create_filelist_frame()
        self.create_preview_frame()

        ttk.Separator(self, orient="horizontal").grid(
            row=1, column=0, columnspan=3, sticky="we")

        self.input_file = tk.StringVar()
        ttk.Entry(self, textvariable=self.input_file).grid(
            row=2, column=0, columnspan=3, sticky="we")

        self.rowconfigure(3, weight=1)
        self.columnconfigure(2, weight=1)
        self.grid(row=0, column=0, sticky="nsew")

    def create_top_frame(self) -> None:
        top_frame = ttk.Frame(self, relief="sunken", padding=4)
        top_frame.columnconfigure(5, weight=1)
        top_frame.rowconfigure(0, weight=1)
        top_frame.grid(row=0, column=0, columnspan=3, sticky="nsew")

        ttk.Label(top_frame, text="Taxi2",
                  font=tkfont.Font(size=20)).grid(row=0, column=0)
        ttk.Separator(top_frame, orient="vertical").grid(
            row=0, column=1, sticky="nsew")

        for image_key, image_file, text, column, command in (
                ("open_button", "open.png", "open", 2, self.open_command),
                ("save_button", "save.png", "save",
                 3, self.save_command("selected")),
                ("save_all_button", "save_all.png",
                 "save_all", 4, self.save_command("all")),
                ("run_button", "run.png", "run", 5, self.run_command)):
            self.images[image_key] = tk.PhotoImage(
                file=os.path.join(sys.path[0], "data", image_file))
            ttk.Button(top_frame, text=text,
                       image=self.images[image_key], compound="top", style="Toolbutton", padding=(10, 0), command=command).grid(row=0, column=column, sticky="w")

        ttk.Separator(top_frame, orient="vertical").grid(
            row=0, column=6, sticky="nsew")
        self.images["logo"] = tk.PhotoImage(file=os.path.join(
            sys.path[0], "data", "iTaxoTools Digital linneaeus MICROLOGO.png"))
        ttk.Label(top_frame, image=self.images["logo"]).grid(
            row=0, column=7, sticky="nse")

    def open_command(self) -> None:
        path = tkfiledialog.askopenfilename()
        if not path:
            return
        self.input_file.set(os.path.abspath(path))

    def save_command(self, which: str) -> Callable[[], None]:
        """
        which should be either "all" or "selected"
        """
        def command():
            save_folder = tkfiledialog.askdirectory()
            if not save_folder:
                return
            for file in self.outfilenames(which):
                full_filename = os.path.join(self.preview_dir, file)
                shutil.copy(full_filename, save_folder)
        return command

    def run_command(self) -> None:
        pass

    def outfilenames(self, which: str) -> Iterator[str]:
        if which == "all":
            index_list = self.filelist.get_children()
        elif which == "selected":
            index_list = self.filelist.selection()
        else:
            raise ValueError(f"Don't know how to save {which}")
        for index in index_list:
            yield self.filelist.item(index, option="text")

    def create_parameters_frame(self) -> None:
        parameters_frame = ttk.LabelFrame(self, text="Parameters")
        parameters_frame.grid(row=3, column=0, sticky="nsew")
        parameters_frame.rowconfigure(5, weight=1)
        parameters_frame.columnconfigure(0, weight=1)

        ttk.Label(parameters_frame, text="Input file format").grid(
            row=0, column=0, sticky='w')

        format_combobox = ttk.Combobox(parameters_frame, textvariable=self.programstate.input_format_name,
                                       state='readonly', values=list(ProgramState.formats.keys()))
        format_combobox.current(0)
        format_combobox.grid(row=1, column=0, sticky='w')

        distances_frm = ttk.LabelFrame(
            parameters_frame, text="Distances to calculate", padding=5, relief='sunken')
        distances_frm.grid(row=2, column=0, sticky='we')

        for kind in range(NDISTANCES):
            ttk.Checkbutton(
                distances_frm, variable=self.programstate.distance_options[kind], text=distances_names[kind]).grid(
                row=kind, column=0, sticky="w")

        ttk.Checkbutton(
            parameters_frame, variable=self.programstate.already_aligned, text="Already aligned").grid(row=3, column=0, sticky='w')
        ttk.Checkbutton(
            parameters_frame, variable=self.programstate.print_alignments, text="Print alignments").grid(row=4, column=0, sticky='w')

    def create_filelist_frame(self) -> None:
        filelist_frame = ttk.Labelframe(self, text="Files")
        filelist_frame.rowconfigure(0, weight=1)
        filelist_frame.columnconfigure(0, weight=1)
        filelist_frame.grid(row=3, column=1, sticky="nsew")

        self.filelist = ttk.Treeview(filelist_frame,
                                     height=15, selectmode="extended", show="tree")
        self.filelist.grid(row=0, column=0, sticky="nsew")

        filelist_scroll = ttk.Scrollbar(filelist_frame,
                                        orient='vertical', command=self.filelist.yview)
        self.filelist.configure(yscrollcommand=filelist_scroll.set)
        filelist_scroll.grid(row=0, column=1, sticky="nsew")

        filelist_scroll_x = ttk.Scrollbar(filelist_frame,
                                          orient='horizontal', command=self.filelist.xview)
        self.filelist.configure(xscrollcommand=filelist_scroll_x.set)
        filelist_scroll_x.grid(row=1, column=0, sticky="nsew")

        self.filelist.bind("<<TreeviewSelect>>", self.preview_selected)

    def icon_for_file(self, filename) -> tk.PhotoImage:
        TXT_EXTS = {".txt", ".tab", ".tsv", ".csv"}
        _, ext = os.path.splitext(filename)
        if ext in TXT_EXTS:
            return self.images["txt_icon"]
        elif ext == ".log":
            return self.images["log_icon"]
        else:
            return self.images["graph_icon"]

    def fill_file_list(self) -> None:
        def by_ext(name):
            name, ext = os.path.splitext(name)
            return (ext, name)

        for filename in sorted(os.listdir(self.preview_dir), key=by_ext):
            name = os.path.basename(filename)
            img = self.icon_for_file(name)
            self.filelist.insert(parent="", index="end", text=name, image=img)

    def create_preview_frame(self) -> None:
        self.preview_frame = ttk.LabelFrame(self, text="Preview")
        self.preview_frame.rowconfigure(0, weight=1)
        self.preview_frame.columnconfigure(0, weight=1)
        self.preview_frame.grid(row=3, column=2, sticky="nsew")

        self.preview = tk.Text(
            self.preview_frame, height=15, width=30, wrap="none")
        self.preview.grid(row=0, column=0, sticky="nsew")

        yscroll = ttk.Scrollbar(
            self.preview_frame, orient='vertical', command=self.preview.yview)
        self.preview.config(yscrollcommand=yscroll.set)
        yscroll.grid(row=0, column=1, sticky="nsew")

        xscroll = ttk.Scrollbar(
            self.preview_frame, orient='horizontal', command=self.preview.xview)
        self.preview.config(xscrollcommand=xscroll.set)
        xscroll.grid(row=1, column=0, sticky="nsew")

    def preview_selected(self, _) -> None:
        self.preview.delete("1.0", "end")
        if not self.filelist.selection():
            return
        selected_index = self.filelist.selection()[-1]
        self.preview_frame.configure(
            text=f'Preview - {self.filelist.item(selected_index, option="text")}')
        file_to_preview = os.path.join(
            self.preview_dir, self.filelist.item(selected_index, option="text"))
        TXT_EXTS = {".txt", ".tab", ".tsv", ".csv", ".log"}
        IMG_EXTS = {".gif", ".png", ".pbm", ".pgm", ".ppm", ".pnm"}
        _, ext = os.path.splitext(file_to_preview)
        if ext in TXT_EXTS:
            self.preview_txt(file_to_preview)
        elif ext in IMG_EXTS:
            self.preview_img(file_to_preview)
        else:
            self.no_preview(file_to_preview)

    def preview_txt(self, filename) -> None:
        with open(filename) as file:
            self.preview.insert("1.0", file.read())

    def preview_img(self, filename) -> None:
        self.images["current"] = tk.PhotoImage(file=filename)
        self.preview.image_create("1.0", image=self.images["current"])

    def no_preview(self, _) -> None:
        self.preview.insert("1.0", "Preview is not possible")


def test_look() -> None:
    root = tk.Tk()
    root.rowconfigure(0, weight=1)
    root.columnconfigure(0, weight=1)
    gui = TaxiGUI(root, preview_dir="/tmp/out_dir")
    gui.fill_file_list()
    root.mainloop()
