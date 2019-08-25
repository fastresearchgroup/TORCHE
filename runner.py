"""
Author(s): Dakota Brown*
*Virginia Commonwealth University
*Fluids in Advanced Systems and Technology Research Group

Script for automating calculations using the TORCHE package based on
values read from a file, whose results are to be either printed to another file,
exported to another Python script, or both.

Date last edited: August 24, 2019
"""


"""-----------------------------------------------------------------------------
Imports
-----------------------------------------------------------------------------"""
import io
import TORCHE as torche
import sys
from TORCHE import Torche


"""-----------------------------------------------------------------------------
Helper method that simplifies the task of running one or more functions on a
list of operands (operands MUST be TORCHE objects).

Method iterates through the list of operands and applies the function, or list
of functions, provided to each operand.

Currently, there is no syntatic sugar added to the results, i.e. the results
are just numbers with no descriptions, units, or rounding. This can be added
later as necessary.
-----------------------------------------------------------------------------"""
def apply(func, operand):
    result = []

    for i in range(len(operand)):
        result.append([f(operand[i]) for f in func])
    return result

"""-----------------------------------------------------------------------------
If script is passed an argument, that argument should be a readable file
containing physical values pertinent to the TORCHE package functions. If an
argument was not passed to the script, the script then prompts the user for such
a file.
-----------------------------------------------------------------------------"""
if (len(sys.argv) > 1):
    filename_in = sys.argv[1]
else:
    filename_in = input("Provide a file to read from: ")

while True:
    try:
        input_file = open(filename_in, 'r')
        break
    except IOError:
        print("The specified file could not be found, or could not be opened.\n")
        filename_in = input("Provide a file to read from: ")

"""-----------------------------------------------------------------------------
Parse file data and create object(s) to represent data
-----------------------------------------------------------------------------"""
# TODO: read values from file (don't yet know the expected format, i.e. XML, JSON)
torche_objs = []

while True: # TODO: change while condition to represent an EOF being reached
    # assignment of values to variables depends on format of text file
    torche_objs.append(Torche(rho, a, b, d, geom, Pr, Pr_w, N_rows, vel, Re))

"""-----------------------------------------------------------------------------
Collect results of running TORCHE functions
-----------------------------------------------------------------------------"""
# TODO: simplify function parameters; pass them as object instead
# TODO: if multiple data points are provided, need to iterate through them
function_list = [torche.dP_Zu, torche.dP_GG, torche.HT_Zu, torche.HT_GG]

results = apply(function_list, torche_objs)
"""-----------------------------------------------------------------------------
Prompt user for next action, i.e. print, export, or both.

***Note:
    it is possible to program this script to accept a secondary argument, one
    that specifies which action to take at this point. This will increase
    effeciency when being automated on lots of data points.
***
-----------------------------------------------------------------------------"""
print("How would you like to proceed?\n")
print("1) Print to a file\n")
print("2) Export to another script\n")
print("3) Both\n")

option = input("Indicate selection: ")

while True: # Loop until a valid selection is made
    if option == 1:
        # Ask for filename to print to
        filename_out = input("Provide a file to write to: ")

        while True:
            try:
                output_file = open(filename_out, 'w')
                break
            except IOError:
                print("The specified file could not be created, or could not be opened.\n")
                filename_out = input("Provide a file to write to: ")

        # Write data to file (need specific format before programming)
        break
    elif option == 2:
        # Ask for script name
        print("placeholder")
        break
    elif option == 3:
        # Do both!
        print("placeholder")
        break
    else:
        print("Invalid selection, please try again.\n")
        option = input("Indicate selection: ")
