#!/usr/bin/python
# -*- coding: utf-8 -*-

from tkinter import *
from tkinter import filedialog as tkFileDialog
from PIL import Image, ImageTk
import os
import script

class ProteinModellerApp(Frame):
  
    def __init__(self, parent):
        Frame.__init__(self, parent)   
        self.parent = parent        
        self.init_ui()
        
    def init_ui(self):
        # Set up the main UI elements and configure the window
        self.parent.title("CIIT Protein Modeller: FragMod")
        self.pack(fill=BOTH, expand=True)
        
        # Create the menu bar and other components
        self.create_menu_bar()
        self.load_default_image()
        self.create_labels()
        self.create_text_input()
        self.create_submit_button()

    def create_menu_bar(self):
        # Create the menu bar with File menu options
        menubar = Menu(self.parent)
        self.parent.config(menu=menubar)
        
        file_menu = Menu(menubar)
        file_menu.add_command(label="Output Directory", underline=0, command=self.output_path)
        file_menu.add_command(label="Exit", underline=0, command=self.on_exit)
        menubar.add_cascade(label="File", menu=file_menu)

    def load_default_image(self):
        # Load and display the default image on the GUI
        self.Path = '1.gif'
        self.img = ImageTk.PhotoImage(Image.open(self.Path))
        label_image = Label(self, image=self.img)
        label_image.grid(row=1, pady=10)

    def create_labels(self):
        # Create and display descriptive labels
        label_text = (
            "We propose a hybrid or integrated methodology (FragMod) to solve\n"
            "the riddle of protein structure prediction, our methodology\n"
            "involves the combinatorial approach by integrating\n"
            "multiple modeling techniques.\n\n"
            "Please Enter or Paste the primary sequence of the query protein:"
        )
        label_description = Label(self, text=label_text, justify=LEFT)
        label_description.grid(row=2, pady=10)

    def create_text_input(self):
        # Create a text input box with a border and increased font size
        text_input = Text(self, bg="white", fg="blue", width=60, height=10, relief="solid", bd=2, font=("Helvetica", 12))
        text_input.grid(row=3, pady=10)
        self.text = text_input  # Store the reference to access text content later

    def create_submit_button(self):
        # Create a submit button
        submit_button = Button(self, text="Submit", command=self.run_modelling)
        submit_button.grid(row=4, pady=10)

    def output_path(self):
        # Open a dialog to select the output directory and save it to a file
        ftypes = [('Text Files', '*.txt'), ('All files', '*')]
        out_path = tkFileDialog.askdirectory(initialdir='.')
        path_file = open('config.txt', 'w')
        path_file.write(str(out_path))
        path_file.close()

    def run_modelling(self):
        # Run the modelling script with the input data from the text box
        input_data = self.get_text_input_data()
        if input_data:
            # Assuming that script.modelit is a function that needs to be called
            # You might need to import the script module at the beginning of your code
            # import script
            # script.modelit(input_data)
            print(f"Running modelling with input data:\n{input_data}")

    def get_text_input_data(self):
        # Retrieve the content from the text input box
        return self.text.get("1.0", END).strip()

    def on_exit(self):
        # Close the application
        self.quit()

def main():
    # Create and run the main Tkinter application
    root = Tk()
    app = ProteinModellerApp(root)
    root.geometry("450x600+100+100")
    root.mainloop()  

if __name__ == '__main__':
    main()
