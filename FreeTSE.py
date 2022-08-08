#coding:utf-8
from util import *

import random as pyrand
from pylab import *
import pandas as pd

import sys, os
import copy
import configparser
import warnings

class FreeTSE:
	def __init__(self, ini=None, name="untitled", gui_mode=False):
		"""Initialization
		
		Parameters
		----------
		ini: str
			Path to the input ini file.
		name: str
			Name of the project.
		"""
		self.name = name
		self.gui_mode = gui_mode
		if ini != None:
			self.read_ini(ini)
		print(self.name)
	
	def read_ini(self, ini):
		"""Read .ini to define default scenario
		
		Parameters
		----------
		ini: str
			Path to the input ini file.
		"""
		self.ini = ini
		try:
			cfg = configparser.ConfigParser()
			cfg.read(ini)
			self.name = cfg["Data"]["name"]
			dt = float(cfg["Resolution"]["dt"])
			dx = float(cfg["Resolution"]["dx"])
			
			mint = float(cfg["Resolution"]["time_min"])
			maxt = float(cfg["Resolution"]["time_max"])
			minx = float(cfg["Resolution"]["space_min"])
			maxx = float(cfg["Resolution"]["space_max"])
			
			number_of_lanes = float(cfg["Resolution"]["number_of_lanes"])
			
			self.speed_data_name = cfg["Data"]["speed"]
			self.speed_label_t = cfg["Data"]["speed_label_t"]
			self.speed_label_x = cfg["Data"]["speed_label_x"]
			self.speed_label_v = cfg["Data"]["speed_label_v"]
			if "density" in cfg["Data"]:
				self.density_data_name = cfg["Data"]["density"]
				self.density_label_t = cfg["Data"]["density_label_t"]
				self.density_label_x = cfg["Data"]["density_label_x"]
				self.density_label_k = cfg["Data"]["density_label_k"]
			else:
				self.density_data_name = None
			if "flow" in cfg["Data"]:
				self.flow_data_name = cfg["Data"]["flow"]
				self.flow_label_t = cfg["Data"]["flow_label_t"]
				self.flow_label_x = cfg["Data"]["flow_label_x"]
				self.flow_label_q = cfg["Data"]["flow_label_q"]
			else:
				self.flow_data_name = None
				
			if "GroundTruth" in cfg and "true_density" in cfg["GroundTruth"] and cfg["GroundTruth"]["true_density"] != "None":
				self.groundtruth = True
				self.density_dat_true_name = cfg["GroundTruth"]["true_density"]
				self.true_density_label_t = cfg["GroundTruth"]["density_label_t"]
				self.true_density_label_x = cfg["GroundTruth"]["density_label_x"]
				self.true_density_label_k = cfg["GroundTruth"]["density_label_k"]
				self.flow_dat_true_name = None
			elif "GroundTruth" in cfg and "true_flow" in cfg["GroundTruth"] and cfg["GroundTruth"]["true_flow"] != "None":
				self.groundtruth = True
				self.density_dat_true_name = None
				self.flow_dat_true_name = cfg["GroundTruth"]["true_flow"]
				self.true_flow_label_t = cfg["GroundTruth"]["flow_label_t"]
				self.true_flow_label_x = cfg["GroundTruth"]["flow_label_x"]
				self.true_flow_label_q = cfg["GroundTruth"]["flow_label_q"]
			else:
				self.groundtruth = False

		except:
			if self.gui_mode == True:
				sg.popup(f"Error: Something is wrong with the .ini file `{ini}`")
			raise Exception(f"Something is wrong with the .ini file `{ini}`")
		
		self.set_data(mint, maxt, dt, minx, maxx, dx, number_of_lanes)
	
	def set_scenario(self, name, dt, dx, mint, maxt, minx, maxx, number_of_lanes, speed_data_name, speed_label_t, speed_label_x, speed_label_v, density_data_name=None, density_label_t=None, density_label_x=None, density_label_k=None, flow_data_name=None, flow_label_t=None, flow_label_x=None, flow_label_q=None, density_dat_true_name=None, true_density_label_t=None, true_density_label_x=None, true_density_label_k=None, flow_dat_true_name=None, true_flow_label_t=None, true_flow_label_x=None, true_flow_label_q=None):
		"""Set estimation scenario.
		
		Parameters
		----------
		name: str
			scenario name
		dt: number
			temporal resolution of estimation (s)
		dx: number
			spatial resolution of estimation (m)
		mint: number
			initial time of the time-space region to be estimated
		maxt: number
			last time of the time-space region to be estimated
		minx: number
			upstream-end position of the time-space region to be estimated
		maxx: number
			downstream-end position of the time-space region to be estimated
		number_of_lanes: number
			the number of lanes of the target section
		speed_data_name: str
			path to speed data (probe vehicle data)
		speed_label_t: str
			column name of the time in the speed data
		speed_label_x: str
			column name of the position in the speed data
		speed_label_v: str
			column name of the speed in the speed data
		density_data_name: str
			path to density data (detector data)
			if it does not exist, set `density_data_name=None`
		density_label_t: str
			column name of the time in the density data
		density_label_x: str
			column name of the position in the density data
		density_label_k: str
			column name of the density in the density data
		flow_data_name: str
			path to flow data (detector data)
			if it does not exist, set `flow_data_name=None`
		flow_label_t: str
		flow_label_x: str
		flow_label_q: str
		density_dat_true_name: str
			path to ground truth density data for accuracy validation
			if it does not exist, set `density_dat_true_name=None`
		true_density_label_t: str
		true_density_label_x: str
		true_density_label_k: str
		flow_dat_true_name: str
			path to ground truth flow data for accuracy validation
			if it does not exist, set `flow_dat_true_name=None`
		true_flow_label_t: str
		true_flow_label_x: str
		true_flow_label_q: str
		"""
		self.name = name
		dt = dt
		dx = dx
		mint = mint
		maxt = maxt
		minx = minx
		maxx = maxx
		number_of_lanes = number_of_lanes
		self.speed_data_name = speed_data_name
		self.speed_label_t = speed_label_t
		self.speed_label_x = speed_label_x
		self.speed_label_v = speed_label_v
		if density_data_name != None:
			self.density_data_name = density_data_name
			self.density_label_t = density_label_t
			self.density_label_x = density_label_x
			self.density_label_k = density_label_k
		else:
			self.density_data_name = "None"
		if flow_data_name != None:
			self.flow_data_name = flow_data_name
			self.flow_label_t = flow_label_t
			self.flow_label_x = flow_label_x
			self.flow_label_q = flow_label_q
		else:
			self.flow_data_name = "None"
		if density_dat_true_name != None:
			self.groundtruth = True
			self.flow_dat_true_name = None
			self.density_dat_true_name = density_dat_true_name
			self.true_density_label_t = true_density_label_t
			self.true_density_label_x = true_density_label_x
			self.true_density_label_k = true_density_label_k
		elif flow_dat_true_name != None:
			self.groundtruth = True
			self.density_dat_true_name = None
			self.flow_dat_true_name = flow_dat_true_name
			self.true_flow_label_t = true_flow_label_t
			self.true_flow_label_x = true_flow_label_x
			self.true_flow_label_q = true_flow_label_q
		else:
			self.groundtruth = False
		
		self.set_data(mint, maxt, dt, minx, maxx, dx, number_of_lanes)
	
	def set_data(self, mint=None, maxt=None, dt=None, minx=None, maxx=None, dx=None, number_of_lanes=None):
		"""Read and structure data
		
		"""
		
		#fundamental parameters
		if mint != None:
			self.mint = mint
		if maxt != None:
			self.maxt = maxt
		if dt != None:
			self.dt = dt
		if minx != None:
			self.minx = minx
		if maxx != None:
			self.MAXX = maxx
		if dx != None:
			self.dx = dx
		if number_of_lanes != None:
			self.number_of_lanes = number_of_lanes
		self.tsize = int((self.maxt-self.mint)/self.dt)
		self.xsize = int((self.MAXX-self.minx)/self.dx)
		
		#precision parameters
		self.detector_precision = 0.01**2
		self.initial_precision = (self.number_of_lanes*0.2)**2
		self.concervation_precision = 0.01**2
		self.cv_speed_precision = 3**2
		self.initial_density = 0.2**2
		
		#traffic data loading
		vv_raw = pd.read_csv(self.speed_data_name)
		vv_mod = vstack([vv_raw[self.speed_label_t], vv_raw[self.speed_label_x], vv_raw[self.speed_label_v]]).T
		self.vv = self.griddata_generation(vv_mod, min_value=0, max_value=self.dx/self.dt, interpolation="both")
		self.probe_record = vv_mod
		
		if self.density_data_name not in [None, "None"]:
			kk_raw = pd.read_csv(self.density_data_name)
			kk_mod = vstack([kk_raw[self.density_label_t], kk_raw[self.density_label_x], kk_raw[self.density_label_k]]).T
			self.kk = self.griddata_generation(kk_raw, missing_value=-1, interpolation="time")
			self.detector_record = kk_mod
		elif self.flow_data_name not in [None, "None"]:
			qq_raw = pd.read_csv(self.flow_data_name)
			qq_mod = vstack([qq_raw[self.flow_label_t], qq_raw[self.flow_label_x], qq_raw[self.flow_label_q]]).T
			self.qq = self.griddata_generation(qq_mod, missing_value=-1, interpolation="time")
			self.kk = self.qq/self.vv
			self.kk[self.qq == -1] = -1
			self.detector_record = qq_mod
		else:
			if self.gui_mode == True:
				sg.popup(f"Error: Either density or flow data is required")
			raise Exception("Either density or flow data is required")
		
		#ground truth data loading
		if self.groundtruth:
			if self.density_dat_true_name != None:
				kk_true_raw = pd.read_csv(self.density_dat_true_name)
				kk_true_mod = vstack([kk_true_raw[self.true_density_label_t], kk_true_raw[self.true_density_label_x], kk_true_raw[self.true_density_label_k]]).T
				self.kk_true = self.griddata_generation(kk_true_mod, missing_value=-1)
			else:
				qq_true_raw = pd.read_csv(self.flow_dat_true_name)
				qq_true_mod = vstack([qq_true_raw[self.true_flow_label_t], qq_true_raw[self.true_flow_label_x], qq_true_raw[self.true_flow_label_q]]).T
				self.qq_true = self.griddata_generation(qq_true_mod, missing_value=-1, interpolation="time")
				self.kk_true = self.qq_true/self.vv
				self.kk_true[self.qq_true == -1] = -1
	
	def griddata_generation(self, data_raw, min_value=None, max_value=None, missing_value=0, interpolation=0):
		"""generate grid data from input table
		
		Parameters
		----------
		data_raw: list
			Raw input traffic data.
		min_value: float or None
			Lowerbound of an output value.
		max_value: float or None
			Upperbound of an output value.
		missing_value: float
			Value used to interpolate missing value.
		interpolation: 0 or "both" or "time" or "space"
			Direction of missing value interpolation.
			
		Returns
		-------
		data_array: array
			Traffic data converted to grid data.
		"""
		data_array = zeros([self.tsize, self.xsize]) + missing_value
		
		data_dic = {}
		for l in data_raw:
			t = l[0]
			x = l[1]
			d = l[2]
			n = int((t-self.mint)//self.dt)
			i = int((x-self.minx)//self.dx)
			if (n,i) in data_dic:
				data_dic[n,i].append(d)
			else:
				data_dic[n,i] = [d]
		for (n,i) in data_dic.keys():
			if 0 <= n < self.tsize and 0 <= i < self.xsize:
				data_array[n,i] = average(data_dic[n,i]) #harmonic mean is suitable depending on the probe vehicle data specification, but we cannot know
				if min_value != None:
					if data_array[n,i] < 0:
						data_array[n,i] = 0
				if max_value != None:
					if data_array[n,i] > self.dx/self.dt:
						speed_exceeds_flag = 1
						print("%.1f>%.1f"%(data_array[n,i],self.dx/self.dt), end=" ")
						data_array[n,i] = self.dx/self.dt
		print()
		
		#missing data interpolation
		num_inter = 0
		if interpolation in ("both", "space"):
			for n in range(self.tsize):
				for i in range(1, self.xsize):
					if data_array[n,i] == missing_value and data_array[n][i-1] != missing_value:
						data_array[n,i] = data_array[n][i-1]
						num_inter += 1
		if interpolation in ("both", "time"):
			for n in range(1, self.tsize):
				for i in range(self.xsize):
					if data_array[n,i] == missing_value and data_array[n-1][i] != missing_value:
						data_array[n,i] = data_array[n-1][i]
						num_inter += 1
			for n in range(self.tsize-2, -1, -1):
				for i in range(self.xsize):
					if data_array[n,i] == missing_value and data_array[n+1][i] != missing_value:
						data_array[n,i] = data_array[n+1][i]
						num_inter += 1
		
		return data_array
	
	def smooth_speeddata(self, tagg, xagg):
		vv_new = zeros([self.tsize, self.xsize])
		for t in range(self.tsize):
			for x in range(self.xsize):
				vlist = []
				for tt in range(tagg):
					for xx in range(xagg):
						if t+tt < self.tsize and x+xx < self.xsize:
							vlist.append(self.vv[t+tt,x+xx])
				vv_new[t,x] = average(vlist)
		self.vv = vv_new

	def filtering(self):
		"""traffic state estimation by Kalman filtering
		"""
		self.yy = self.kk
		
		self.k_prio = zeros([self.tsize, self.xsize])
		self.k_post = zeros([self.tsize, self.xsize])
		self.k_smoo = zeros([self.tsize, self.xsize])
		
		xx_prio = ones(self.xsize)*self.initial_density
		xx_post = ones(self.xsize)*self.initial_density
		xx_smoo = ones(self.xsize)*self.initial_density
		
		V_prio = zeros([self.xsize, self.xsize])
		V_post = ones([self.xsize, self.xsize])*self.initial_precision
		V_smoo = zeros([self.xsize, self.xsize])
		#Q = identiry(self.xsize])*self.concervation_precision
		Q =  zeros([self.xsize, self.xsize])
		
		self.k_prio[0] = xx_prio
		self.k_post[0] = xx_post
		
		self.V_prio_dic = {}
		self.V_post_dic = {}
		self.V_smoo_dic = {}
		self.xx_prio_dic = {}
		self.xx_post_dic = {}
		self.xx_smoo_dic = {}
		self.F_dic = {}
		
		for n in range(self.tsize):
			
			F = zeros([self.xsize, self.xsize])
			
			for i in range(self.xsize):
				if i > 0:
					F[i,i] = (1-self.dt/self.dx*self.vv[n,i])
					F[i,i-1] = self.dt/self.dx*self.vv[n,i-1]
					Q[i,i] = self.cv_speed_precision
					Q[i,i-1] = self.cv_speed_precision
				else:
					F[0,0] = 1
					Q[0,0] = self.cv_speed_precision
			
			H = zeros([self.xsize, self.xsize])
			R = zeros([self.xsize, self.xsize])
			for i in range(self.xsize):
				if self.kk[n,i] != -1:
					H[i,i] = 1
					R[i,i] = self.detector_precision
			
			xx_prio = F @ xx_post
			
			V_prio = F @ V_post @ F.T + Q
			
			K = V_prio @ H.T @ pinv(H @ V_prio @ H.T + R)
			
			xx_post = xx_prio + K @ (self.yy[n] - H @ xx_prio)
			
			V_post = V_prio - K @ H @ V_prio
			
			#ad-hoc cleansing
			for i in lange(xx_prio):
				if xx_prio[i] < 0:
					xx_prio[i] = 0
				#if xx_prio[i] > 0.2*self.number_of_lanes:
				#	xx_prio[i] = 0.2*self.number_of_lanes
				if xx_post[i] < 0:
					xx_post[i] = 0
				#if xx_post[i] > 0.2*self.number_of_lanes:
				#	xx_post[i] = 0.2*self.number_of_lanes

			self.k_prio[n] = xx_prio
			self.k_post[n] = xx_post
			
			self.V_prio_dic[n] = copy.copy(V_prio)
			self.V_post_dic[n] = copy.copy(V_post)
			self.xx_prio_dic[n] = copy.copy(xx_prio)
			self.xx_post_dic[n] = copy.copy(xx_post)
			self.F_dic[n] = copy.copy(F)
		
		self.V_smoo_dic[self.tsize-1] = copy.copy(V_post)
		self.xx_smoo_dic[self.tsize-1] = copy.copy(xx_post)
		self.k_smoo[self.tsize-1] = xx_post
	
	def smoothing(self):
		"""traffic state estimation by RST smoothing
		"""
		for n in range(self.tsize-2, -1, -1):
			A = self.V_post_dic[n] @ self.F_dic[n].T @ pinv(self.V_prio_dic[n+1])
			xx_smoo = self.xx_post_dic[n] + A @ (self.xx_smoo_dic[n+1] - self.xx_prio_dic[n+1])
			V_smoo = self.V_post_dic[n] + A @ (self.V_smoo_dic[n+1] - self.V_prio_dic[n+1]) @ A.T

			#without smoothing
			#xx_smoo = self.xx_post_dic[n] 
			#V_smoo = self.V_post_dic[n] 

			#ad-hoc cleansing
			for i in lange(xx_smoo):
				if xx_smoo[i] < 0:
					xx_smoo[i] = 0
				#if xx_smoo[i] > 0.2*self.number_of_lanes:
				#	xx_smoo[i] = 0.2*self.number_of_lanes
			
			
			self.k_smoo[n] = xx_smoo
			
			self.xx_smoo_dic[n] = copy.copy(xx_smoo)
			self.V_smoo_dic[n] = copy.copy(V_smoo)
		
		if sum(xx_smoo.flatten()) == sum(self.xx_post_dic[n].flatten()):
			print("no smoothing!")
			
	
	def estimation(self):
		"""call filtering and smoothing
		"""
		self.filtering()
		self.smoothing()
		self.compute_cum_curves()
	
	def compute_cum_curves(self):
		"""compute cumulative arrival/depature curves
		"""
		self.N = zeros([self.tsize, self.xsize])
		for t in range(self.tsize):
			if t > 0:
				self.N[t,0] = self.N[t-1,0] + self.k_smoo[t,0]*self.vv[t,0]*self.dt
			for x in range(1, self.xsize):
				self.N[t,x] = self.N[t,x-1] - self.k_smoo[t,x]*self.dx
	
	def accuracy_evaluation(self, print_mode=1):
		"""evaluate accuracy against ground truth data
		"""
		if self.groundtruth:
			if self.density_dat_true_name != None:
				if print_mode:
					print("density error")
				self.rmse_prio = rmse(self.kk_true[self.kk_true>0].flatten(), self.k_prio[self.kk_true>0].flatten())
				self.rmse_post = rmse(self.kk_true[self.kk_true>0].flatten(), self.k_post[self.kk_true>0].flatten())
				self.rmse_smoo = rmse(self.kk_true[self.kk_true>0].flatten(), self.k_smoo[self.kk_true>0].flatten())
				
				self.mape_prio = mae(self.kk_true[self.kk_true>0].flatten(), self.k_prio[self.kk_true>0].flatten(), 1)*100
				self.mape_post = mae(self.kk_true[self.kk_true>0].flatten(), self.k_post[self.kk_true>0].flatten(), 1)*100
				self.mape_smoo = mae(self.kk_true[self.kk_true>0].flatten(), self.k_smoo[self.kk_true>0].flatten(), 1)*100
				
			else:
				if print_mode:
					print("flow error")
				self.rmse_prio = rmse(self.qq_true[self.qq_true>0].flatten(), (self.vv[self.qq_true>0]*self.k_prio[self.qq_true>0]).flatten())
				self.rmse_post = rmse(self.qq_true[self.qq_true>0].flatten(), (self.vv[self.qq_true>0]*self.k_post[self.qq_true>0]).flatten())
				self.rmse_smoo = rmse(self.qq_true[self.qq_true>0].flatten(), (self.vv[self.qq_true>0]*self.k_smoo[self.qq_true>0]).flatten())
				
				self.mape_prio = mae(self.qq_true[self.qq_true>0].flatten(), (self.vv[self.qq_true>0]*self.k_prio[self.qq_true>0]).flatten(), 1)*100
				self.mape_post = mae(self.qq_true[self.qq_true>0].flatten(), (self.vv[self.qq_true>0]*self.k_post[self.qq_true>0]).flatten(), 1)*100
				self.mape_smoo = mae(self.qq_true[self.qq_true>0].flatten(), (self.vv[self.qq_true>0]*self.k_smoo[self.qq_true>0]).flatten(), 1)*100
				
			if print_mode:
				#print("RMSE (prio): %.4f"%self.rmse_prio)
				#print("RMSE (post): %.4f"%self.rmse_post)
				print("RMSE: %.4f"%self.rmse_smoo)

				#print("MAPE (prio): %.1f%%"%self.mape_prio)
				#print("MAPE (post): %.1f%%"%self.mape_post)
				print("MAPE: %.1f%%"%self.mape_smoo)
	
	def save_results(self, name):
		"""output estimation results as csv
		"""
		out = [["t (s)", "x (m)", "estimated k (veh/m)", "v (m/s)", "estimated q (veh/s)", "true k (veh/m)"]]
		if self.groundtruth and self.density_dat_true_name == None:
			out[0][-1] = "true q (veh/s)"
		if not self.groundtruth:
			out[0] = out[0][:-1]
		for n in range(self.tsize):
			for i in range(self.xsize):
				l = [
					n*self.dt,
					i*self.dx,
					self.k_smoo[n,i],
					self.vv[n,i],
					self.k_smoo[n,i]*self.vv[n,i],
				]
				if self.groundtruth:
					l.append(self.kk_true[n,i] if self.density_dat_true_name != None else self.qq_true[n,i])
				else:
					pass
				out.append(l)
		writecsv(name, out)
	
	def get_results(self):
		return self.k_smoo*self.vv, self.k_smoo, self.vv
	
	def visualize(self, speed=0, prior=0, posterior=0, smooth=0, posterior_stddev=0, smooth_stddev=0, true=0, observation=0, qk=0, qmax=0.2, cumcurves=0, qmax_cum=0, scatter=0, cum_true=0, timeseries=0, save=0, inputdata=0, fname="img"):
		"""visualize the results
		"""
		if prior:
			figure(figsize=(12,3))
			title("a prior")
			subplots_adjust(top=.9, bottom=.2, left=.125, right=.9, wspace=.2, hspace=.2)
			imshow(self.k_prio.T, origin="lower", interpolation="nearest", aspect="auto", extent=(self.mint, self.maxt, self.minx, self.MAXX), vmin=0, vmax=0.1*self.number_of_lanes, cmap=cm_kusakabe_pb2)
			colorbar().set_label("density (veh/m)")
			xlabel("t (s)")
			ylabel("x (m)")
			grid()
			if save:
				savefig("%s_prior.png"%fname)
		
		if posterior:
			figure(figsize=(12,3))
			title("estimated density")
			subplots_adjust(top=.9, bottom=.2, left=.125, right=.9, wspace=.2, hspace=.2)
			imshow(self.k_post.T, origin="lower", interpolation="nearest", aspect="auto", extent=(self.mint, self.maxt, self.minx, self.MAXX), vmin=0, vmax=0.1*self.number_of_lanes, cmap=cm_kusakabe_pb2)
			colorbar().set_label("estimated density (veh/m)")
			xlabel("t (s)")
			ylabel("x (m)")
			grid()
			if save:
				savefig("%s_posterior.png"%fname)
		
		if smooth:
			figure(figsize=(12,3))
			title("estimated density")
			subplots_adjust(top=.9, bottom=.2, left=.125, right=.9, wspace=.2, hspace=.2)
			imshow(self.k_smoo.T, origin="lower", interpolation="nearest", aspect="auto", extent=(self.mint, self.maxt, self.minx, self.MAXX), vmin=0, vmax=0.1*self.number_of_lanes, cmap=cm_kusakabe_pb2)
			colorbar().set_label("estimated density (veh/m)")
			xlabel("t (s)")
			ylabel("x (m)")
			grid()
			if save:
				savefig("%s_estimated_density.png"%fname)
		
		if posterior_stddev:
			vars = []
			for n in range(1,self.tsize):
				vars.append([self.V_post_dic[n][i,i] for i in range(self.xsize)])
			figure(figsize=(12,3))
			title("stddev of estimated density")
			subplots_adjust(top=.9, bottom=.2, left=.125, right=.9, wspace=.2, hspace=.2)
			imshow(sqrt(array(vars)).T, origin="lower", interpolation="nearest", aspect="auto", extent=(self.mint, self.maxt, self.minx, self.MAXX), vmin=0, vmax=0.05*self.number_of_lanes, cmap=cm_kusakabe_pb2)
			colorbar().set_label("posterior density stddev (veh/m)")
			xlabel("t (s)")
			ylabel("x (m)")
			grid()
			if save:
				savefig("%s_posterior_stddev.png"%fname)
		
		if smooth_stddev:
			vars = []
			for n in range(1,self.tsize):
				vars.append([self.V_smoo_dic[n][i,i] for i in range(self.xsize)])
			figure(figsize=(12,3))
			title("stddev of estimated density")
			subplots_adjust(top=.9, bottom=.2, left=.125, right=.9, wspace=.2, hspace=.2)
			imshow(sqrt(array(vars)).T, origin="lower", interpolation="nearest", aspect="auto", extent=(self.mint, self.maxt, self.minx, self.MAXX), vmin=0, vmax=0.05*self.number_of_lanes, cmap=cm_kusakabe_pb2)
			colorbar().set_label("estimated density stddev (veh/m)")
			xlabel("t (s)")
			ylabel("x (m)")
			grid()
			if save:
				savefig("%s_smooth_stddev.png"%fname)
		
		if true:
			if self.groundtruth:
				figure(figsize=(12,3))
				subplots_adjust(top=.9, bottom=.2, left=.125, right=.9, wspace=.2, hspace=.2)
				title("true density")
				imshow(self.kk_true.T, origin="lower", interpolation="nearest", aspect="auto", extent=(self.mint, self.maxt, self.minx, self.MAXX), vmin=0, vmax=0.1*self.number_of_lanes, cmap=cm_kusakabe_pb2)
				colorbar().set_label("density (veh/m)")
				xlabel("t (s)")
				ylabel("x (m)")
				grid()
				if save:
					savefig("%s_true_density.png"%fname)
		
		if observation:
			figure(figsize=(12,3))
			subplots_adjust(top=.9, bottom=.2, left=.125, right=.9, wspace=.2, hspace=.2)
			title("observed density")
			imshow(self.kk.T, origin="lower", interpolation="nearest", aspect="auto", extent=(self.mint, self.maxt, self.minx, self.MAXX), vmin=0, vmax=0.1*self.number_of_lanes, cmap=cm_kusakabe_pb2)
			colorbar().set_label("density (veh/m)")
			xlabel("t (s)")
			ylabel("x (m)")
			grid()
			if save:
				savefig("%s_observation.png"%fname)
		
		if speed:
			figure(figsize=(12,3))
			subplots_adjust(top=.9, bottom=.2, left=.125, right=.9, wspace=.2, hspace=.2)
			title("observed speed")
			imshow(self.vv.T, origin="lower", interpolation="nearest", aspect="auto", extent=(self.mint, self.maxt, self.minx, self.MAXX), vmin=0, vmax=self.dx/self.dt, cmap=cm_kusakabe_pb2)
			colorbar().set_label("speed (m/s)")
			xlabel("t (s)")
			ylabel("x (m)")
			grid()
			if save:
				savefig("%s_speed.png"%fname)
			
		if inputdata:
			figure(figsize=(12,3))
			subplots_adjust(top=.9, bottom=.2, left=.125, right=.9, wspace=.2, hspace=.2)
			title("input data")
			plot(self.probe_record[:,0], self.probe_record[:,1], "b.", label="speed")
			plot(self.detector_record[:,0], self.detector_record[:,1], "rx", label="detector")
			xlim([self.mint, self.maxt])	
			ylim([self.minx, self.MAXX])
			xlabel("t (s)")
			ylabel("x (m)")
			legend()
			grid()
			if save:
				savefig("%s_inputdata.png"%fname)
		
		if qk:
			figure()
			k_max = 0.1*self.number_of_lanes
			q_max = qmax
			
			q_max = percentile(self.qq_true.flatten(), 95)*1.5
			k_max = percentile(self.kk_true.flatten(), 95)*1.5
			
			hist2d(self.k_smoo[self.k_smoo>0].flatten(), (self.k_smoo[self.k_smoo>0]*self.vv[self.k_smoo>0]).flatten(), bins=20, range=[[0,k_max],[0,q_max]], cmap=cm_kusakabe_pb2)
			colorbar().set_label("flequency")
			xlabel("estimated k (veh/m)")
			ylabel("estimated q (veh/s)")
			grid()
			if save:
				savefig("%s_qk.png"%fname)
		
		if cumcurves:
			tt = arange(self.mint, self.maxt, self.dt)
			
			figure(figsize=(12,3))
			plot(tt, self.N[:,0], "r", label="arrival")
			plot(tt, self.N[:,-1], "b", label="departure")
			grid()
			xlabel("t (s)")
			ylabel("cumulative N (veh)")
			legend()
			
			figure(figsize=(12,3))
			plot(tt, self.N[:,0]-self.N[:,-1], "r", label="# of vehicles")
			ylim(ymin=0)
			grid()
			xlabel("t (s)")
			ylabel("existing N (veh)")
			
			if qmax_cum == 0:
				qmax = (self.N[-1,-1] - self.N[0,-1])/(self.maxt-self.mint)
			figure(figsize=(12,3))
			plot(tt, self.N[:,0]-(tt-self.mint)*qmax, "r", label="arrival")
			plot(tt, self.N[:,-1]-(tt-self.mint)*qmax, "b", label="departure")
			plot([self.mint, self.mint+self.dt*self.xsize], [0, -self.dt*self.xsize*qmax], "--", c="gray", label="ref. slope %.3f"%(-qmax))
			for t in linspace(self.mint, self.maxt, 10):
				plot([t, t+self.dt*self.xsize], [0, -self.dt*self.xsize*qmax], "--", c="gray")
			grid()
			xlabel("t (s)")
			ylabel("oblique N (veh)")
			legend()
		
		#if cum_true and self.groundtruth:


		if scatter and self.groundtruth:
			k_max = 0.1*self.number_of_lanes
#			figure(figsize=(4,4))
#			subplot(111, aspect="equal")
#			title("a prior")
#			hist2d(self.kk_true[self.kk_true>0].flatten(), self.k_prio[self.kk_true>0].flatten(), bins=20, range=[[0,k_max],[0,k_max]], cmap=cm_kusakabe_pb2)
#			colorbar().set_label("flequency")
#			plot([0,k_max],[0,k_max],"r--")
#			xlabel("true density")
#			ylabel("smoothed density")
#			grid()
#			
#			figure(figsize=(4,4))
#			subplot(111, aspect="equal")
#			title("a posterior")
#			hist2d(self.kk_true[self.kk_true>0].flatten(), self.k_post[self.kk_true>0].flatten(), bins=20, range=[[0,k_max],[0,k_max]], cmap=cm_kusakabe_pb2)
#			colorbar().set_label("flequency")
#			plot([0,k_max],[0,k_max],"r--")
#			xlabel("true density")
#			ylabel("smoothed density")
#			grid()
			
			k_max = percentile(self.kk_true.flatten(), 95)*1.5
			
			figure(figsize=(4,4))
			subplot(111, aspect="equal")
			subplots_adjust(top=.9, bottom=.1, left=.2, right=.9, wspace=.2, hspace=.2)
			#title("estimates (MAP/smoothed)")
			hist2d(self.kk_true[self.kk_true>0].flatten(), self.k_smoo[self.kk_true>0].flatten(), bins=20, range=[[0,k_max],[0,k_max]], cmap=cm_kusakabe_pb2)
			colorbar().set_label("flequency")
			plot([0,k_max],[0,k_max],"r--")
			xlabel("true density (veh/m)")
			ylabel("estimated density (veh/m)")
			grid()
			if save:
				savefig("%s_k_scatter.png"%fname)
		
			if self.density_dat_true_name == None:
				q_max = qmax
				
				q_max = percentile(self.qq_true.flatten(), 95)*1.5
				
				figure(figsize=(4,4))
				subplot(111, aspect="equal")
				subplots_adjust(top=.9, bottom=.1, left=.2, right=.9, wspace=.2, hspace=.2)
				#title("estimates (MAP/smoothed)")
				hist2d(self.qq_true[self.qq_true>0].flatten(), (self.vv[self.qq_true>0]*self.k_smoo[self.qq_true>0]).flatten(), bins=20, range=[[0,q_max],[0,q_max]], cmap=cm_kusakabe_pb2)
				colorbar().set_label("flequency")
				plot([0,q_max],[0,q_max],"r--")
				xlabel("true flow (veh/s)")
				ylabel("estimated flow (veh/s)")
				grid()
				if save:
					savefig("%s_q_scatter.png"%fname)
		
		if timeseries and self.groundtruth:
			for x in lange(self.kk_true[0]):
				if self.kk_true[0,x] != -1:
					if self.density_dat_true_name != None:
						figure(figsize=(12,4))
						title("x=%.0f"%(self.dx*x))
						plot(arange(self.mint, self.maxt, self.dt), self.kk_true[:,x], "r", label="ground truth")
						plot(arange(self.mint, self.maxt, self.dt), self.k_smoo[:,x], "b--", label="estimated")
						legend(loc="best")
						xlabel("t")
						ylabel("flow")
						ylim(ymin=0)
					else:
						figure(figsize=(12,4))
						title("x=%.0f"%(self.dx*x))
						plot(arange(self.mint, self.maxt, self.dt), self.qq_true[:,x], "r", label="ground truth")
						plot(arange(self.mint, self.maxt, self.dt), self.vv[:,x]*self.k_smoo[:,x], "b--", label="estimated")
						legend(loc="best")
						xlabel("t")
						ylabel("flow")
						ylim(ymin=0)
					if save:
						savefig("%s_timeseries_%d.png"%(fname, x))
		
		show()

if __name__ == "__main__":
	tse = FreeTSE()
	tse.set_scenario(
		name = "ngsim_trajectory",
		dt = 4,
		dx = 100,
		mint = 0,
		maxt = 800,
		minx = 0,
		maxx = 500,
		number_of_lanes = 5,
		speed_data_name = "./dat/ngsim_sampled_trajectories.csv",
		speed_label_t = "t",
		speed_label_x = "x",
		speed_label_v = "v",
		density_data_name = None,
		density_label_t = "t",
		density_label_x = "x",
		density_label_k = "k",
		flow_data_name = "./dat/ngsim_grid_flow_200m.csv",
		flow_label_t = "t",
		flow_label_x = "x",
		flow_label_q = "q",
		density_dat_true_name = None,
		true_density_label_t = "t",
		true_density_label_x = "x",
		true_density_label_k = "k",
		flow_dat_true_name = "./dat/ngsim_grid_flow_400m.csv",
		true_flow_label_t = "t",
		true_flow_label_x = "x",
		true_flow_label_q = "q"
	)
	tse.estimation()
	tse.accuracy_evaluation()
	fname = "res_test"
	tse.save_results(fname+".csv")
	tse.visualize(smooth=1, true=1, speed=1, timeseries=1, inputdata=1, save=1, fname=fname)
	
	q, k, v = tse.get_results()
	print("flow", q)
	print("density", k)
	print("speed", v)