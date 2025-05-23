#coding:utf-8
from .util import * # This was already a relative import
import PySimpleGUIWx as sg
import sys, datetime

from .FreeTSE import FreeTSE # Changed from import * to specific import

def run_gui_application():
    # Original logic from the if __name__ == "__main__": block
    
    # Part 1: ini選択 (INI selection)
    sg.theme("Default1")
    
    menu_def = [] # Empty menu
    
    layout = [
        [sg.Menu(menu_def, tearoff=True)],
        [sg.Text("FreeTSE", font=("Helvetica", 18))],
        [sg.Frame(title="Estimation scenario", layout=[
                [sg.Text("Choose .ini file")], 
                [sg.InputText("", key="ini", size=(33, 1)), sg.FileBrowse()]
            ], relief=sg.RELIEF_SUNKEN)
        ],
        [sg.Frame(title="Visualization: Estimation results", layout=[
                [sg.CBox("Time-space diagram of estimated density", key="tsd_esti", default=True)],
                [sg.CBox("Time-space diagram of probe vehicle speed", key="tsd_speed", default=True)],
                [sg.CBox("Time-space diagram of input data", key="tsd_input", default=True)],
            ], relief=sg.RELIEF_SUNKEN)
        ],
        [sg.Frame(title="Visualization: Validation resutls (if ground truth is available)", layout=[
                [sg.CBox("Time-space diagram of true density", key="tsd_true", default=True)],
                [sg.CBox("Scatter diagrams of true and estimated states", key="scat", default=True)],
                [sg.CBox("Time series of true and estimated states", key="timeseries", default=True)],
            ], relief=sg.RELIEF_SUNKEN)
        ],
        [sg.Frame(title="Save options", layout=[
                [sg.CBox("Save the results", key="save", default=True)],
                [sg.Text("Save file name prefix:"), sg.InputText("res_"+datetime.datetime.now().strftime("%Y%m%d_%H%M%S"), key="fname", size=(20, 1))],
            ], relief=sg.RELIEF_SUNKEN)
        ],
        [sg.Submit(button_text="Estimate"), sg.Cancel()]
    ]
    
    window = sg.Window("FreeTSE", layout, default_element_size=(40, 1), grab_anywhere=False)
    event, values = window.read()
    
    if event in [None, "Cancel"]:
        window.close()
        return # Exit the function if Cancel or window closed

    window.close() # Close the window after getting initial values, or keep it open?
                   # Original code closed it here. Let's stick to that for now.

    # Part 2: 実行 (Execution)
    ini = values["ini"]
    if not ini: # Check if ini is empty
        sg.popup("Error: Please choose your .ini file.")
        # raise Exception("Please choose your .ini file") # Consider if raising an exception is best for GUI
        return # Exit if no ini file selected

    tsd_true = values["tsd_true"]
    tsd_mean = values["tsd_esti"] # tsd_esti was the key for "estimated density"
    tsd_speed = values["tsd_speed"]
    tsd_data = values["tsd_input"] # tsd_input was for "input data"
    timeseries_val = values["timeseries"] # key was "timeseries"
    scat_val = values["scat"] # key was "scat"
    export = values["save"]
    fname_prefix = values["fname"] # key was "fname"
    
    # Original code printed the INI file content, this might be for debugging, can be kept or removed.
    # For a library function, printing to stdout might not be ideal.
    # print(f"Selected INI file: {ini}") 
    # f = open(ini, "r")
    # for l in f:
    #     print(l[:-1])
    # f.close()
    
    # print("\nEstimating, please wait...\n") # Similar to above, console output.
                                          # GUI should provide feedback if possible.

    try:
        tse = FreeTSE(ini=ini, gui_mode=True) # Pass ini path
        tse.estimation()
        tse.accuracy_evaluation() # This method has a print_mode, default is 1 (prints to console)
        if export:
            tse.save_results(fname_prefix + ".csv")
        
        # Call visualize. Ensure visualize can handle sg.popup if gui_mode is True and errors occur.
        # The visualize method in FreeTSE.py already calls show() which displays the plots.
        tse.visualize(
            smooth=tsd_mean, 
            true=tsd_true, 
            speed=tsd_speed, 
            timeseries=timeseries_val, 
            inputdata=tsd_data, 
            scatter=scat_val, 
            save=export, 
            fname=fname_prefix
        )
    except Exception as e:
        sg.popup(f"An error occurred during processing:\n{e}")
        # Optionally re-raise or log the exception here
        # For a GUI, just showing the popup might be the intended behavior.

# Example of how it might be called (for testing, not part of the library code itself)
# if __name__ == '__main__':
# run_gui_application()
