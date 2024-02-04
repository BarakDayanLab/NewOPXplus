import tkinter as tk


class BDDialog:

    def __init__(self):
        self.result = None
        pass

    def handle_button(self, button):
        self.res_text = self.entry.get("1.0", tk.END)
        self.res_button = button

        # Destroy the dialog window
        self.root.destroy()
        pass

    def prompt(self, title=None, message=None, button1_text='', button2_text=''):
        # Create a root window
        self.root = tk.Tk()
        if title is not None:
            self.root.title("Experiment completed")
        self.root.geometry("500x270")
        self.root.wm_resizable(False, False)

        # Create a dialog window
        #self.dialog = tk.Toplevel(root)
        self.dialog = tk.Frame(self.root)

        self.dialog.pack(fill="both", expand=True)

        # Create a label with a message
        label = tk.Label(self.dialog, text='Enter the comment to save with experiment', font=("Bold", 12))  # "Bold", 12)
        label.pack(pady=10)

        # Create an entry with the StringVar
        self.entry = tk.Text(self.dialog, height=10)
        self.entry.pack(side='top', padx=10, pady=10)

        # Create three buttons with custom labels and commands
        button1 = tk.Button(self.dialog, text=button1_text, font=("Bold", 12), command=lambda: self.handle_button(button1_text))
        button1.pack(side='right', padx=5, pady=5)

        button2 = tk.Button(self.dialog, text=button2_text, font=("Bold", 12), command=lambda: self.handle_button(button2_text))
        button2.pack(side='right', padx=5, pady=5)

        # button3 = tk.Button(self.dialog, text='Something', command=self.handle_button)
        # button3.pack(side='right', padx=5, pady=5)

        # Start the main loop
        self.root.mainloop()

        return self.res_text, self.res_button

if __name__ == "__main__":

    bdd = BDDialog()

    res_text, res_button = bdd.prompt(title='Experiment Completed', message='Enter the comment to save with experiment', button1_text='Ignore', button2_text='Shager!')

    pass
