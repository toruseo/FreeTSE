#coding:utf-8
from util import *
import PySimpleGUIWx as sg
import sys, datetime

from FreeTSE import *

if __name__ == "__main__":
	if 1:
		#ini選択
		sg.theme("Default1")
		
		menu_def = []
		
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
		window.close()
		
		#sg.popup(f"event: {event}\n values: {values}")
		
		if event in [None, "Cancel"]:
			sys.exit()

	if 1:
		#実行
		ini = values["ini"]
		if ini == "":
			sg.popup(f"Error: Please choose your .ini file.")
			raise Exception(f"Please choose your .ini file")
		tsd_true = values["tsd_true"]
		tsd_mean = values["tsd_esti"]
		tsd_speed = values["tsd_speed"]
		tsd_data = values["tsd_input"]
		timeseries = values["timeseries"]
		scat = values["scat"]
		export = values["save"]
		fname = values["fname"]
		
		f = open(ini, "r")
		for l in f:
			print(l[:-1])
		f.close()
		
		print("\nEstimating, please wait...\n")
		
		tse = FreeTSE(ini, gui_mode=True)
		tse.estimation()
		tse.accuracy_evaluation()
		if export:
			tse.save_results(fname+".csv")
		tse.visualize(smooth=tsd_mean, true=tsd_true, speed=tsd_speed, timeseries=timeseries, inputdata=tsd_data, scatter=scat, save=export, fname=fname)
