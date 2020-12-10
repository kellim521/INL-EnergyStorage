# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 22:03:05 2020

@author: Riley
"""
import PowerBlock_Riley
import tkinter as tk

#From PowerBlock_Riley.py
#This line runs the PowerBlock_Riley.py file to get the values of every temperature.
M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11,M12,M13,M14,A,B,C,D,E,F,A2,B2,D2,E2,F2= PowerBlock_Riley.Cycle(1.0,155)

#These lines convert the temperature values to integers with two decimal places.
M1 = str(round(M1, 2))
M2 = str(round(M2, 2))
M3 = str(round(M3, 2))
M4 = str(round(M4, 2))
M5 = str(round(M5, 2))
M6 = str(round(M6, 2))
M7 = str(round(M7, 2))
M8 = str(round(M8, 2))
M9 = str(round(M9, 2))
M10 = str(round(M10, 2))
M11 = str(round(M11, 2))
M12 = str(round(M12, 2))
M13 = str(round(M13, 2))
M14 = str(round(M14, 2))
A = str(round(A, 2))
B = str(round(B, 2))
C = str(round(C, 2))
D = str(round(D, 2))
E = str(round(E, 2))
F = str(round(F, 2))
A2 = str(round(A2, 2))
B2 = str(round(B2, 2))
D2 = str(round(D2, 2))
E2 = str(round(E2, 2))
F2 = str(round(F2, 2))

#These lines are the code for the GUI. 
#The text2.insert on line 54 is where you can change what is being displayed on the right of the GUI.
root = tk.Tk()
text1 = tk.Text(root, height=50, width=120)
photo = tk.PhotoImage(file='./Power Block.png')
text1.insert(tk.END, '\n')
text1.image_create(tk.END, image=photo)
text1.pack(side=tk.LEFT)
text2 = tk.Text(root, height=50, width=35)
scroll = tk.Scrollbar(root, command=text2.yview)
text2.configure(yscrollcommand=scroll.set)
text2.tag_configure('big', font=('Arial', 15))
text2.insert(tk.END,"\nM1 = " + str(M1) + "\nM2 = " + str(M2) + "\nM3 = " + str(M3) + "\nM4 = " + str(M4) + "\nM5 = " + str(M5) + "\nM6 = " + str(M6) + "\nM7 = " + str(M7) + "\nM8 = " + str(M8) + "\nM9 = " + str(M9) + "\nM10 = " + str(M10) + "\nM11 = " + str(M11) + "\nM12 = " + str(M12) + "\nM13 = " + str(M13) + "\nM14 = " + str(M14) + "\nA = " + str(A) + "\nB = " + str(B) + "\nC = " + str(C) + "\nD = " + str(D) + "\nE = " + str(E) + "\nF = " + str(F) + "\nA2 = " + str(A2) + "\nB2 = " + str(B2) + "\nD2 = " + str(D2) + "\nE2 = " + str(E2) + "\nF2 = " + str(F2) + "\n", 'big')
text2.pack(side=tk.LEFT)
scroll.pack(side=tk.RIGHT, fill=tk.Y)
root.mainloop()