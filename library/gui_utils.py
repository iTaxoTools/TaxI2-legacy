import os.path
import tkinter as tk
import tkinter.filedialog
import tkinter.ttk as ttk
import tkinter.messagebox
from typing import Literal, Any, Dict, Tuple, List, Optional
from contextlib import contextmanager
import warnings


def try_relpath(path: str) -> str:
    """
    Return relative path to the file, if possible
    """
    try:
        return os.path.relpath(path)
    except:
        return os.path.abspath(path)


class FileChooser():
    """
    Creates a frame with a label, entry and browse button for choosing files
    """

    def __init__(self, parent: Any, *, label: str, mode: Literal["open", "save"]):
        self.frame = ttk.Frame(parent)
        self.frame.columnconfigure([0, 1], weight=1)
        self.label = ttk.Label(self.frame, text=label)
        self.file_var = tk.StringVar()
        self.entry = ttk.Entry(self.frame, textvariable=self.file_var)
        if mode == "open":
            self._dialog = tk.filedialog.askopenfilename
        elif mode == "save":
            self._dialog = tk.filedialog.asksaveasfilename

        def browse() -> None:
            if (newpath := self._dialog()):
                try:
                    newpath = os.path.relpath(newpath)
                except:
                    newpath = os.path.abspath(newpath)
                self.file_var.set(newpath)

        self.button = ttk.Button(self.frame, text="Browse", command=browse)

        self.label.grid(row=0, column=0, sticky='nws')
        self.entry.grid(row=1, column=0, sticky='nwse')
        self.button.grid(row=1, column=1)
        self.grid = self.frame.grid


class LabeledEntry():
    """
    Group of a label, entry and a string variable
    """

    def __init__(self, parent: tk.Misc, *, label: str):
        self.frame = ttk.Frame(parent)
        self.label = ttk.Label(self.frame, text=label)
        self.var = tk.StringVar()
        self.entry = ttk.Entry(self.frame, textvariable=self.var)
        self.frame.columnconfigure(1, weight=1)
        self.label.grid(column=0, row=0)
        self.entry.grid(column=1, row=0)
        self.grid = self.frame.grid


class LabeledCombobox():
    """
    Group of a label, Combobox and a string variable
    """

    def __init__(self, parent: tk.Misc, *, label: str, values: List[str], readonly: bool):
        self.frame = ttk.Frame(parent)
        self.label = ttk.Label(self.frame, text=label)
        self.var = tk.StringVar()
        self.combobox = ttk.Combobox(
            self.frame, textvariable=self.var, values=values)
        if readonly:
            self.combobox.configure(state='readonly')
            if values:
                self.combobox.current(0)
        self.frame.columnconfigure(1, weight=1)
        self.label.grid(column=0, row=0)
        self.combobox.grid(column=1, row=0)
        self.grid = self.frame.grid


class Listbox():
    """
    Wrapper for a read-only tk.Listbox with a method that returns the selection
    """

    def __init__(self, parent: tk.Misc, *, height: int, selectmode: Literal["browse", "extended"], values: List[str]) -> None:
        self.list = values
        self.listbox = tk.Listbox(
            parent, height=height, selectmode=selectmode, listvariable=tk.StringVar(value=values))
        self.grid = self.listbox.grid

    def selection(self) -> List[str]:
        """
        Returns the list of items that are currently selected
        """
        return [self.list[i] for i in self.listbox.curselection()]


class ColumnSelector():

    def __init__(self, parent: tk.Misc) -> None:
        self.notebook = ttk.Notebook(parent)
        self.lists: List[Listbox] = []
        self.grid = self.notebook.grid

    def set_columns(self, columns: Dict[str, List[str]]) -> None:
        for i in range(self.notebook.index("end")):
            self.notebook.forget(i)
        self.lists = []
        for column_name in columns:
            self.lists.append(Listbox(self.notebook, height=10,
                                      values=columns[column_name], selectmode="extended"))
            self.notebook.add(self.lists[-1].listbox, text=column_name)

    def selection(self) -> Optional[Tuple[str, List[str]]]:
        if self.notebook.index("end") == 0:
            return None
        i = self.notebook.index("current")
        return (self.notebook.tab(i)["text"], self.lists[i].selection())


class FileListChooser():
    """
    Contains a listbox with file names and "Add Files", "Add Directory" and "Remove Files" buttons
    """

    def __init__(self, parent: tk.Misc, *, label: str) -> None:
        self.frame = ttk.Frame(parent)
        self.frame.rowconfigure(4, weight=1)
        self.frame.columnconfigure(0, weight=1)
        self.label = ttk.Label(self.frame, text=label)
        self.listbox = tk.Listbox(self.frame, height=10, selectmode='extended')
        self.yscroll = ttk.Scrollbar(
            self.frame, orient='vertical', command=self.listbox.yview)
        self.listbox.configure(yscrollcommand=self.yscroll.set)
        self.add_files_btn = ttk.Button(
            self.frame, text="Add Files", command=self.add_files)
        self.add_directory_btn = ttk.Button(
            self.frame, text="Add Directory", command=self.add_directory)
        self.remove_files_btn = ttk.Button(
            self.frame, text="Remove Files", command=self.remove_files)

        self.label.grid(row=0, column=0)
        self.listbox.grid(row=1, column=0, rowspan=4, sticky="nsew")
        self.yscroll.grid(row=1, column=1, rowspan=4, sticky="ns")
        self.add_files_btn.grid(row=1, column=2, sticky="w")
        self.add_directory_btn.grid(row=2, column=2, sticky="w")
        self.remove_files_btn.grid(row=3, column=2, sticky="w")

        self.grid = self.frame.grid

    def add_files(self) -> None:
        self.listbox.insert('end', *map(
            try_relpath, tk.filedialog.askopenfilenames()))

    def add_directory(self) -> None:
        if not (dirname := tk.filedialog.askdirectory()):
            return
        self.listbox.insert('end', *map(try_relpath, os.listdir(dirname)))

    def remove_files(self) -> None:
        for i in reversed(self.listbox.curselection()):
            self.listbox.delete(i)

    def file_list(self) -> List[str]:
        return self.listbox.get(0, 'end')


@contextmanager
def display_errors_and_warnings() -> Any:
    try:
        with warnings.catch_warnings(record=True) as warns:
            yield
            for w in warns:
                tk.messagebox.showwarning("Warning", str(w.message))
    except FileNotFoundError as ex:
        tk.messagebox.showerror("Error", ex.strerror)
    # except Exception as ex:
    #     tk.messagebox.showerror("Error", str(ex))
