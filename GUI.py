import tkinter as tk
from tkinter import filedialog, messagebox
from pathlib import Path
import graph as gp
import KEGG_0_3 as kg
import translate as tr

# Translation helper 
def translate_compound(text):
    text = tr.ukr_to_eng(text)
    text = tr.amino_acids(text)
    return text


# Main window
root = tk.Tk()
root.title("KEGG Pathway Drawer")
root.geometry("540x480")
root.resizable(True, True)

# Directory variable  
script_dir = Path(__file__).parent.resolve()
dir_path_var = tk.StringVar(value=str(script_dir))

# Outer frame with padding 
outer = tk.Frame(root, padx=16, pady=16)
outer.pack(expand=True, fill="both")

# Section: Start compound 
tk.Label(outer, text="Start compound:", anchor="w").grid(
    row=0, column=0, sticky="w", pady=4)
entry_start = tk.Entry(outer, width=34)
entry_start.insert(0, "glucose")
entry_start.grid(row=0, column=1, columnspan=2, sticky="ew", pady=4, padx=4)

# Section: Organism
tk.Label(outer, text="Organism:", anchor="w").grid(
    row=1, column=0, sticky="w", pady=4)
entry_organism = tk.Entry(outer, width=34)
entry_organism.insert(0, "Escherichia coli")
entry_organism.grid(row=1, column=1, columnspan=2, sticky="ew", pady=4, padx=4)

# Mode selector
mode_var = tk.StringVar(value="single")

mode_frame = tk.LabelFrame(outer, text="Mode", padx=8, pady=4)
mode_frame.grid(row=2, column=0, columnspan=3, sticky="ew", pady=6)

tk.Radiobutton(mode_frame, text="Single target",
               variable=mode_var, value="single",
               command=lambda: toggle_mode()).pack(side="left", padx=8)
tk.Radiobutton(mode_frame, text="Multiple targets",
               variable=mode_var, value="multi",
               command=lambda: toggle_mode()).pack(side="left", padx=8)

# Single-target panel  
single_frame = tk.Frame(outer)
single_frame.grid(row=3, column=0, columnspan=3, sticky="ew")

tk.Label(single_frame, text="Target compound:", anchor="w").grid(
    row=0, column=0, sticky="w", pady=4)
entry_single_target = tk.Entry(single_frame, width=34)
entry_single_target.insert(0, "піруват")
entry_single_target.grid(row=0, column=1, sticky="ew", pady=4, padx=4)
single_frame.columnconfigure(1, weight=1)

# Multi-target panel
multi_frame = tk.Frame(outer)
# (not gridded yet — shown only in multi mode)

tk.Label(multi_frame, text="Target compounds:", anchor="w").grid(
    row=0, column=0, sticky="nw", pady=4)

# Scrollable list of target entries
target_list_frame = tk.Frame(multi_frame, bd=1, relief="sunken")
target_list_frame.grid(row=0, column=1, sticky="ew", padx=4)
multi_frame.columnconfigure(1, weight=1)

target_entries = []   # holds tk.Entry widgets


def add_target_row(default_text=""):
    """Append a new target-compound row to the list."""
    row_frame = tk.Frame(target_list_frame)
    row_frame.pack(fill="x", pady=2)

    entry = tk.Entry(row_frame, width=28)
    entry.insert(0, default_text)
    entry.pack(side="left", padx=(4, 0))
    target_entries.append(entry)

    def remove_this():
        target_entries.remove(entry)
        row_frame.destroy()

    btn_minus = tk.Button(row_frame, text="−", width=2,
                          command=remove_this,
                          bg="#EF9A9A", relief="flat")
    btn_minus.pack(side="left", padx=4)


def add_target_via_button():
    add_target_row()


btn_add = tk.Button(multi_frame, text="＋  Add target",
                    command=add_target_via_button,
                    bg="#A5D6A7", relief="flat", padx=4)
btn_add.grid(row=1, column=1, sticky="w", padx=4, pady=4)

# Seed with two default rows
add_target_row("glycerin")
add_target_row("Palmitate")


def toggle_mode():
    if mode_var.get() == "single":
        multi_frame.grid_remove()
        single_frame.grid(row=3, column=0, columnspan=3, sticky="ew")
    else:
        single_frame.grid_remove()
        multi_frame.grid(row=3, column=0, columnspan=3, sticky="ew")


# Start in single mode
toggle_mode()

# Directory selector
dir_frame = tk.Frame(outer)
dir_frame.grid(row=4, column=0, columnspan=3, sticky="ew", pady=6)
dir_frame.columnconfigure(1, weight=1)

tk.Label(dir_frame, text="Output folder:", anchor="w").grid(
    row=0, column=0, sticky="w", pady=4)
dir_entry = tk.Entry(dir_frame, textvariable=dir_path_var, state="readonly", width=28)
dir_entry.grid(row=0, column=1, sticky="ew", pady=4, padx=4)


def select_directory():
    folder = filedialog.askdirectory(title="Select output folder")
    if folder:
        dir_path_var.set(folder)


tk.Button(dir_frame, text="Browse…", command=select_directory).grid(
    row=0, column=2, padx=4)

# Nodes-per-row control
layout_frame = tk.Frame(outer)
layout_frame.grid(row=5, column=0, columnspan=3, sticky="w", pady=2)

tk.Label(layout_frame, text="Nodes per row:").pack(side="left")
nodes_per_row_var = tk.IntVar(value=4)
tk.Spinbox(layout_frame, from_=2, to=10,
           textvariable=nodes_per_row_var, width=4).pack(side="left", padx=6)

# Cache management
def clear_cache():
    if messagebox.askyesno(
            "Clear graph cache",
            "This deletes compound_graph_with_reactions.pkl and enzyme_cache.pkl "
            "so they are rebuilt from KEGG on the next run.\n\n"
            "Do this after updating the app to fix reaction-ID bugs, "
            "or when you want a fresh graph for a different organism.\n\n"
            "Continue?"):
        kg.clear_graph_cache()
        messagebox.showinfo("Done", "Cache cleared. The next search will rebuild it from KEGG.")

tk.Button(outer, text="🗑  Clear graph cache", command=clear_cache,
          fg="#B71C1C", relief="flat").grid(row=6, column=2, sticky="e", pady=4)

# Execute button
def run_process():
    start    = translate_compound(entry_start.get().strip())
    organism = entry_organism.get().strip()
    out_dir  = dir_path_var.get()
    n_per_row = nodes_per_row_var.get()

    if not start:
        messagebox.showerror("Input error", "Please enter a start compound.")
        return
    if not organism:
        messagebox.showerror("Input error", "Please enter an organism.")
        return

    try:
        if mode_var.get() == "single":
            # Single-target mode
            target = translate_compound(entry_single_target.get().strip())
            if not target:
                messagebox.showerror("Input error", "Please enter a target compound.")
                return

            out_path = out_dir + "/gui_output"
            result   = kg.find_path(start, target, organism)
            print(result)
            gp.create_pathway_pdf(result, n_per_row, out_path)
            messagebox.showinfo("Done", f"PDF saved to:\n{out_path}.pdf")

        else:
            # Multi-target mode
            raw_targets = [e.get().strip() for e in target_entries if e.get().strip()]
            if not raw_targets:
                messagebox.showerror("Input error", "Please add at least one target compound.")
                return

            translated_targets = [translate_compound(t) for t in raw_targets]

            out_path = out_dir + "/gui_multi_output"
            all_paths = kg.find_paths(start, translated_targets, organism)

            # Filter out empty paths (compounds unreachable) and warn
            valid_paths  = [p for p in all_paths if p]
            valid_names  = [t for t, p in zip(translated_targets, all_paths) if p]
            failed_names = [t for t, p in zip(translated_targets, all_paths) if not p]

            if failed_names:
                messagebox.showwarning(
                    "Some paths not found",
                    f"No path found for:\n  {chr(10).join(failed_names)}\n\n"
                    "These will be skipped."
                )

            if not valid_paths:
                messagebox.showerror("Error", "No valid paths found for any target.")
                return

            # Infer start name from first step of first valid path
            start_display = valid_paths[0][0][0] if valid_paths[0] else start

            gp.create_multi_target_pathway(
                pathway_lists=valid_paths,
                target_names=valid_names,
                start_name=start_display,
                nodes_per_row=n_per_row,
                filename=out_path
            )
            messagebox.showinfo("Done", f"PDF saved to:\n{out_path}.pdf")

    except Exception as exc:
        messagebox.showerror("Error", str(exc))
        raise   # also print traceback to console


btn_run = tk.Button(outer, text="▶  Execute",
                    command=run_process,
                    bg="#2196F3", fg="white",
                    font=("Helvetica", 11, "bold"),
                    padx=12, pady=6, relief="flat")
btn_run.grid(row=6, column=0, columnspan=3, pady=16)

outer.columnconfigure(1, weight=1)

if __name__ == "__main__":
    root.mainloop()
