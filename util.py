#coding:utf-8

def lange(l):
	"""range(len(l))を略すためだけの関数
	
	Parameters
	----
	l : list
	"""
	return range(len(l))

def readcsv(name, mode="str", header_row=0):
	"""
	CSVファイルを読み込む関数
	name: 読み込み対象ファイル
	mode: 読み込みモード
		"str" or 指定なし: string
		"int": int
		"float": float
		その他： 自動判別
	header_row: 1行目のヘッダを飛ばす場合1
	"""
	
	out = []
	in_file = open(name, "r")
	
	for line in in_file:
		if header_row == 1:
			#ヘッダ行を飛ばす
			header_row = 0
			continue
		
		if mode == "str":
			out.append(line[:-1].split(","))
		elif mode == "int":
			out.append([int(i) for i in line[:-1].split(",")])
		elif mode == "float":
			out.append([float(i) for i in line[:-1].split(",")])
		else:
			out_tmp2 = []
			if line[-1] in [str(i) for i in range(10)]: 
				line_tmp = line
			else:
				line_tmp = line[:-1]
			for i in line_tmp.split(","):
				try:
					out_tmp2.append(int(i))
				except ValueError:
					try:
						out_tmp2.append(float(i))
					except ValueError:
						out_tmp2.append(i)
			out.append(out_tmp2)
	in_file.close()
	return out

def writecsv(name,data):
	"""
	CSVファイルを書き込む関数
	name: 書き込み対象ファイル
	data: 書き込む2次元配列
	"""
	out_file=open(name,"w")
	for line in data:
		for j in range(len(line)):
			# 書き込む1行の値が最後の場合には改行コードを
			# 挿入し，それ以外はタブコードを挿入する
			if j==len(line)-1:
				out_file.write(str(line[j])+"\n")
			else:
				out_file.write(str(line[j])+",")
	out_file.close()

def ave(data, zero=None):
	if len(data) == 0:
		return zero
	
	return sum(data)/len(data)


def mae(data1, data2, p=0, nodata=0):
	if p == 0:
		return ave([abs(data1[i]-data2[i]) for i in lange(data1)], nodata)
	else:
		return ave([abs((data1[i]-data2[i])/data1[i]) for i in lange(data1) if data1[i] != 0], nodata)

from math import sqrt

def rmse(data1, data2, p=0, nodata=0):
	if p == 0:
		return sqrt(ave([(data1[i]-data2[i])**2 for i in lange(data1)], nodata))
	else:
		return sqrt(ave([((data1[i]-data2[i])/data1[i])**2 for i in lange(data1) if data1[i] != 0], nodata))

def define_colormap(red, green, blue, gradationsteps=1024, name=None):
	"""
	Colorbar from
	- Kusakabe, T., Iryo, T. and Asakura, Y.: Barcelo, J. and Kuwahara, M. (Eds.) Data mining for traffic flow analysis: Visualization approach, Traffic Data Collection and its Standardization, Springer, 57-72, 2010
	"""
	from matplotlib.colors import LinearSegmentedColormap
	cdict = {"red"  : red, "green": green, "blue" : blue}
	if name == None:
		name = "cmap"
	return LinearSegmentedColormap(name, cdict, gradationsteps)

cm_kusakabe_pb2 = define_colormap(( (0,1,1),(0.01,1,220/255), (1/6,198/255,198/255), (2/6,175/255,175/255), (3/6,149/255,149/255), (4/6,121/255,121/255), (5/6, 87/255, 87/255), (1,  2/255,  2/255) ),
								  ( (0,1,1),(0.01,1,220/255), (1/6,189/255,189/255), (2/6,159/255,159/255), (3/6,129/255,129/255), (4/6,100/255,100/255), (5/6, 71/255, 71/255), (1, 37/255, 37/255) ),
								  ( (0,1,1),(0.01,1,220/255), (1/6,218/255,218/255), (2/6,215/255,215/255), (3/6,212/255,212/255), (4/6,208/255,208/255), (5/6,203/255,203/255), (1,197/255,197/255) ))
