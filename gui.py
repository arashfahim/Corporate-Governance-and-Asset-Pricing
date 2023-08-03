import tkinter as tk
import colorsys
import random
from tkinter import messagebox
import ast
import equations as eqn





class ParamEntry(object):
    """"Base class for the user interface"""
    def __init__(self):
        self.j = 0
        self.i = 0
        self.color = ['red', 'green', 'blue', 'yellow', 'purple', 'orange']
        self.params_dict = {
            "\u03BC:": '1.2', #mu
            "\u03B3:": '4.0', #gamma
            "\u0072:": '3.0', #r
            "\u03BB": '0.5', #lambda
            "\u03C3:": r'[0.4, 0.2, 0.1]', #sigma
            "\u03C1:": r'[0.0, 0.5, 1.0]' #rho
        }
        self.entry_width = int(len(self.params_dict["\u03C3:"])/5*3.5)
        # print(self.entry_width)
        self.num_dict  = {
            "range": "1.0",
            "Number_of_points": "2000",
            "df0": "5",
            "Stop_criterion_for_F_ODE": "0.001",
            "Upper_dy_lim": "1000",
            "Lower_dy_lim": "0.1",
            "num_iterations": "200"
        }
        self.sim_dict = {
            "num_paths": "100",
            "max_length": "500"
        }
        self.dict ={
            "Parameters" : self.params_dict,
            "Numerical_settings" : self.num_dict,
            # "Simulation_settings" : self.sim_dict,#add if simulatioins added to the
            "color": "red"
        }
        self.frame_list = []
        self.dict_list = []
        self.j = 0
        self.lbl_y = 5
        # self.range = 1.7
        # self.df0 = 4
        # self.Upper_dy_lim = 80
        # self.Lower_dy_lim = 0.1
        # self.num_iterations = 1000
        # self.eps = 0.008
        # self.N = 500
        # self.num_paths = 100
        # self.max_length = 500

    def center(self):
        w = self.window.winfo_reqwidth()
        h = self.window.winfo_reqheight()
        ws = self.window.winfo_screenwidth()
        hs = self.window.winfo_screenheight()
        x = (ws/5) - (w/5)
        y = (hs/5) - (h/5)
        self.window.geometry('+%d+%d' % (x, y))


    '''creates table of parameters'''
    def dict_2_table(self,dict,string):
        lbl = tk.Label(
            master = self.tmp_box
        )
        lbl.config(
                relief=tk.RIDGE,
                text=string,
                bg="black",
                fg = "white",
                width = 30,
                font=("Arial", 14)
        )
        lbl.grid(row = self.i, columnspan=2, padx=1, sticky="N")
        self.i += 1
        for name, val in dict.items():
            lbl = tk.Label(
                    master=self.tmp_box
            )
            lbl.config(
                    relief=tk.RIDGE,
                    text=name,
                    bg="gray",
                    anchor="w",
                    width=22,
                    font=("Arial", 14)
                )
            lbl.grid(row = self.i, column = 0, padx=1, sticky="W")
            entry = tk.Entry(
                    master=self.tmp_box
            )
            entry.config(
                    relief=tk.GROOVE,
                    bg="black",
                    fg = "white",
                    insertbackground='white',
                    width = self.entry_width,
                    font=("Arial", 14)
                )
            entry.insert("0",val)
            entry.grid(row =self.i , column = 1, padx=1)
            self.i += 1
        self.i += 1
        self.j += 1

    def initiate(self):
        self.window = tk.Tk()
        self.window.title("Principal Agent parameters")
        self.window.geometry("1000x600")
        self.center()
        self.base_frm = tk.Frame(
            master = self.window,
            width=200,
            bg="white"
        )
        self.base_frm.grid(row = 0, columnspan=6, padx=1, sticky="W")
        self.param_box()
        btn_frm = tk.Frame(
            master=self.window)
        btn_frm.config(
            bg = "white"
        )
        btn_frm.grid(row = 1, column = 0, sticky= "W", padx=1, pady=1)

        add_btn = tk.Button(master=btn_frm,
            text="New Frame (up to 6)",
            width=20,
            height=5,
            bg="blue",
            fg="red",
            font=("Arial", 20),
            command = self.param_box
        )
        add_btn.grid(row = 0, column = 0, sticky= "W", padx=1, pady=1)

        save_btn = tk.Button(master=btn_frm,
            text="Save Parameters",
            width=20,
            height=5,
            bg="red",
            fg="blue",
            font=("Arial", 20),
            command = lambda:
                self.get_all_entry_widgets_text_content(self.frame_list)
        )
        save_btn.grid(row = 0, column = 1, sticky= "W")


    def param_box(self):
        '''Creates a window to show and modify values of parameters and other quantities'''
        try:
            color = next(iter(self.color))
            del self.color[0]
        except Exception as e:
            raise r"You cannot add more than 6 set of parameters. Kick on 'Save Parameters' to continue."

        self.tmp_box = tk.Frame(
            master = self.base_frm,
            bg=color
        )
        self.tmp_box.grid(row = self.i, column = self.j, sticky="W")
        self.frame_list.append(self.tmp_box)
        for j in range(0,len(self.dict)-1):
            string = list(self.dict.keys())[j]
            dict = self.dict[string]
            # print(dict)
            self.dict_2_table(dict,string)
        self.i = 0

    def get_all_entry_widgets_text_content(self,parent_frames):
        if self.dict_list:
            self.dict_list = []
        """Stopped at 07/21/2021, 5:38pm to do Vijay corrections to the math paper."""
        tmp_dict = dict.fromkeys([
            "Parameters",
            "Numerical_settings",
            "Simulation_settings"
        ])
        for frame in parent_frames:
            children_widgets = frame.winfo_children()
            tmp_list=[]
            for child_widget in children_widgets:
                if child_widget.winfo_class() == 'Entry':
                    tmp_list.append(child_widget.get())
            params_dict = {
                list(self.params_dict.keys())[i]: ast.literal_eval(tmp_list[i]) for i in list(range(0,6))
            }
            num_dict  = {
                list(self.num_dict.keys())[i]: ast.literal_eval(tmp_list[i+6]) for i in list(range(0,7))
            }
            # sim_dict = {
            #     list(self.sim_dict.keys())[i]: ast.literal_eval(tmp_list[i+13]) for i in list(range(0,2))
            # }#no simulation exists in the code we add when a simulation is necessary , sim_dict
            tmp_list = [params_dict, num_dict]#no simulation exists in the code we add when a simulation is necessary , sim_dict
            tmp_dict = {
                list(tmp_dict.keys())[i]: tmp_list[i] for i in range(0,len(tmp_list))
            }
            tmp_dict["color"] = frame["background"]
            self.dict_list.append(tmp_dict)
        print("Parameters saved.")

    def on_closing(self):
        if messagebox.askokcancel("Quit",
            "Do you want to quit and run the rest of the program?"):
            self.window.destroy()

    def quit(self):
        self.window.destroy()
