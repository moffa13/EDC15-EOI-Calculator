import os
import struct


class DoubleLinkedList():
	def __init__(self, initList = []):
		self.list = []
		self.head = None
		self.tail = None
	
		for i in initList:
			self.add(i)

	def add(self, data):
		if self.head is None:
			node = DLLNode(self.list, data)
			self.list.append(node)
			self.head = node
			self.tail = node
		else:
			node = DLLNode(self.list, data, self.tail, None)
			self.list.append(node)
			self.tail.next = node
			self.tail = node


class DLLNode():
	def __init__(self, listT, data = 0, prev = None, next = None):
		self.data = data
		self.prev = prev
		self.list = listT
		self.i = len(self.list)
		self.next = next
	def skip(self, n):
		n = int(n)
		return self.list[self.i + n]


class Map:
	def __init__(self, edcList = None, address = 0, options = {}):
		self.edcList = edcList
		self.address = address
		self.options = options
		if edcList is not None:
			self.load()

	def writeToFile(self, file):
		if self.isMapWithoutAxis():
			# ONLY DATA LIKE SOI
			file = open(os.path.join(os.path.dirname(__file__), file), 'rb+')
			file.seek(self.address) # go to data because file is displacement is in bytes

			for y in range(self.yAxisSize):
				for x in range(self.xAxisSize):
					val = int(round((self.get(x, y) - self.options["map"]["add"]) / self.options["map"]["multiplier"]))
					v = struct.pack('<h', val)
					file.write(v)
			file.close()
		else:
			print("NOT YET IMPLEMENTED !")
	def __str__(self):
		ret = "    "
		for i in range(self.xAxisSize):
			axisValue = self.xAxis[i]
			ret += "{0:10}".format(axisValue)
		ret += "\n\n"
		for y in range(self.yAxisSize):
			yAxisValue = self.yAxis[y]
			ret += "{0:>4} ".format(yAxisValue)
			for x in range(self.xAxisSize):
				realValue = self.get(x, y)
				ret += "{:10.2f}".format(realValue)
			ret += "\n"

		return ret

	'''
	Takes two axis values with corresponding values and returns the interpolated between the two
	yMin, yMax are the axis values
	min, max are the values
	y is the value we want to interpolate between yMin and yMax
	'''
	@staticmethod 
	def rawInterpolate(min, max, yMin, yMax, y):
		diffMin = max - min
		diffY = yMax - yMin
		if diffY == 0 or diffMin == 0:
			return min # or max ?
		realDiff = y - yMin
		res = (diffMin / diffY) * realDiff
		return min + res

	'''
	Returns from the axis values from which ones we should interpolate from
	'''
	@staticmethod 
	def findUpperLowerBoundaries(headerList, x):
		if not headerList:
			return None

		if x <= headerList[0]:
			return (headerList[0], headerList[0])

		if x >= headerList[-1]:
			return (headerList[-1], headerList[-1])

		for i in range(1, len(headerList)):
			if x == headerList[i]:
				return (headerList[i], headerList[i])
			elif x < headerList[i]:
				return (headerList[i - 1], headerList[i])

		return (headerList[-1], headerList[-1])
	
	'''
	Returns value from a 3d map, None if no exact match
	'''
	def getValue(self, x, y):
		for i in range(self.xAxisSize):
			xAxisValue = self.xAxis[i]
			if xAxisValue == x :
				for j in range(self.yAxisSize):
					yAxisValue = self.yAxis[j]
					if yAxisValue == y:
						return self.get(i, j)
		return None
	'''
	Returns a value from interpolating in a map
	'''
	def interpolate(self, reqX, reqY):

		xRange = Map.findUpperLowerBoundaries(self.xAxis, reqX)
		yRange = Map.findUpperLowerBoundaries(self.yAxis, reqY)

		if xRange[0] == xRange[1] and yRange[0] == yRange[1]:
			return self.getValue(xRange[0], yRange[0])
		elif xRange[0] == xRange[1]:
			v1 = self.getValue(xRange[0], yRange[0])
			v2 = self.getValue(xRange[0], yRange[1])
			return Map.rawInterpolate(v1, v2, yRange[0], yRange[1], reqY)
		elif yRange[0] == yRange[1]:
			v3 = self.getValue(xRange[0], yRange[1])
			v4 = self.getValue(xRange[1], yRange[1])
			return Map.rawInterpolate(v3, v4, xRange[0], xRange[1], reqX)
		else:
			v1 = self.getValue(xRange[0], yRange[0])
			v2 = self.getValue(xRange[0], yRange[1])
			v3 = self.getValue(xRange[1], yRange[0])
			v4 = self.getValue(xRange[1], yRange[1])
			xInterpolate = Map.rawInterpolate(v1,v3, xRange[0], xRange[1], reqX)
			xInterpolate2 = Map.rawInterpolate(v2,v4, xRange[0], xRange[1], reqX)
			return Map.rawInterpolate(xInterpolate, xInterpolate2, yRange[0], yRange[1], reqY)
	def setInterpolate(self, reqX, reqY, target):
		x1, x2 = Map.findUpperLowerBoundaries(self.xAxis, reqX)
		y1, y2 = Map.findUpperLowerBoundaries(self.yAxis, reqY)

		q11 = self.interpolate(x1, y1)
		q21 = self.interpolate(x2, y1)
		q12 = self.interpolate(x1, y2)
		q22 = self.interpolate(x2, y2)
		
		return Map.rawSetInterpolate(reqX, reqY, x1, x2, y1, y2, q11, q21, q12, q22, target)
	@staticmethod
	def rawSetInterpolate(reqX, reqY, x1, x2, y1, y2, q11, q21, q12, q22, target):

		dx = x2 - x1
		dy = y2 - y1

		if dx == 0 and dy == 0:
			q11 = target
		elif dx == 0:
			wy = (reqY - y1) / dy
			interpolated = q11 * (1 - wy) + q12 * wy
			delta = target - interpolated
			w1 = 1 - wy
			w2 = wy
			norm = w1**2 + w2**2
			q11 += delta * w1 / norm
			q12 += delta * w2 / norm

			q21 = None
			q22 = None
			
		elif dy == 0:
			wx = (reqX - x1) / dx
			interpolated = q11 * (1 - wx) + q21 * wx
			delta = target - interpolated
			w1 = 1 - wx
			w2 = wx
			norm = w1**2 + w2**2
			q11 += delta * w1 / norm
			q21 += delta * w2 / norm

			q12 = None
			q22 = None
			
		else:
			wx = (reqX - x1) / dx
			wy = (reqY - y1) / dy

			interpolated = (
				q11 * (1 - wx) * (1 - wy) +
				q21 * wx * (1 - wy) +
				q12 * (1 - wx) * wy +
				q22 * wx * wy
			)

			delta = target - interpolated

			w11 = (1 - wx) * (1 - wy)
			w21 = wx * (1 - wy)
			w12 = (1 - wx) * wy
			w22 = wx * wy

			norm = w11**2 + w21**2 + w12**2 + w22**2

			q11 += delta * w11 / norm
			q21 += delta * w21 / norm
			q12 += delta * w12 / norm
			q22 += delta * w22 / norm

			

		return q11, q21, q12, q22, x1, x2, y1, y2

	def load(self):
		listStart = self.edcList.head
		mapStart = listStart.skip(self.address / 2)
		mapStart = mapStart.next # skip number descriptor

		if self.isNormalMap():

			mapStart = self.loadY(mapStart)
			mapStart = mapStart.next
			mapStart = self.loadX(mapStart)
			self.loadData(mapStart)
		elif self.isMapWithoutAxis():
			mapStart = listStart.skip(self.options["x"]["address"] / 2)
			mapStart = mapStart.next
			self.loadX(mapStart)

			mapStart = listStart.skip(self.options["y"]["address"] / 2)
			mapStart = mapStart.next
			self.loadY(mapStart)

			mapStart = listStart.skip(self.address / 2)
			self.loadData(mapStart)
		elif self.is1DMap():
			mapStart = self.loadX(mapStart)
			self.yAxisSize = 1
			self.yAxis = [0]
			self.loadData(mapStart)
		else:
			raise Exception("Unknown error")


	def is1DMap(self):
		return "onlyX" in self.options and self.options["onlyX"]

	def isMapWithoutAxis(self):
		return "onlyData" in self.options and self.options["onlyData"]

	def isNormalMap(self):
		return not self.is1DMap() and not self.isMapWithoutAxis()

	def loadX(self, node):
		self.xAxisSize = node.data
		self.xAxis = []
		node = node.next
		for i in range(self.xAxisSize):
			v = node.data * self.options["x"]["multiplier"] + self.options["x"]["add"]
			if "func" in self.options["x"]:
				v = self.options["x"]["func"](v)
			self.xAxis.append(v)
			node = node.next
		return node

	def loadY(self, node):
		self.yAxisSize = node.data
		self.yAxis = []
		node = node.next
		for i in range(self.yAxisSize):
			v = node.data * self.options["y"]["multiplier"] + self.options["y"]["add"]
			if "func" in self.options["y"]:
				v = self.options["y"]["func"](v)
			self.yAxis.append(v)
			node = node.next
		return node

	def loadData(self, node=None):
		self.data = []
		
		for i in range(self.yAxisSize):
			self.data.append([])
			for j in range(self.xAxisSize):
				if node is not None:
					v = node.data * self.options["map"]["multiplier"] + self.options["map"]["add"]
					node = node.next
				else:
					v = None
				if "map" in self.options and "func" in self.options["map"]:
					v = self.options["map"]["func"](v)
				self.data[i].append(v)
				
		return node
	def copyAxis(self, otherMap):
		self.xAxis = otherMap.xAxis
		self.yAxis = otherMap.yAxis
		self.xAxisSize = otherMap.xAxisSize
		self.yAxisSize = otherMap.yAxisSize
		self.loadData()
	def setV(self, x, y, v):
		self.data[y][x] = v
	def get(self, x, y):
		return self.data[y][x]

def findInjectionfromSOI(durationsArray, selectorMap, iq, rpm, soi):
	durationsReversed = selectorMap.xAxis.copy()
	durationsReversed.reverse() # get selector axis in low to high order (0,5,10,15,....)

	bound = Map.findUpperLowerBoundaries(durationsReversed, soi) # find 2 duration maps associated with that SOI

	fromDur = selectorMap.getValue(bound[0], 0)
	toDur = selectorMap.getValue(bound[1], 0) 

	fromDurValue = durationsArray[fromDur].interpolate(rpm, iq) # in durations, mg is Y, rpm is X, opposite of SOI
	toDurValue = durationsArray[toDur].interpolate(rpm, iq) # in durations, mg is Y, rpm is X, opposite of SOI

	realInjectionVal = Map.rawInterpolate(fromDurValue, toDurValue, bound[0], bound[1], soi)
	return realInjectionVal

def displayNeededChanges(q11, q21, q12, q22, x1, x2, y1, y2):

	print("Do the changes as follows:")
	if x1 == x2 and y1 == y2:
		print(f"{'':>6}{x1:>8}")
		print("-" * 16)
		print(f"{y1:>5} | {q11:.2f}")
	elif x1 == x2:
		print(f"{'':>6}{x1:>8}")
		print("-" * 16)
		print(f"{y1:>5} | {q11:.2f}")
		print(f"{y2:>5} | {q12:.2f}")
	elif y1 == y2:
		print(f"{'':>6}{x1:>8}{x2:>8}")
		print("-" * 26)
		print(f"{y1:>5} | {q11:.2f}   {q21:.2f}")
	else:
		print(f"{'':>6}{x1:>8}{x2:>8}")
		print("-" * 26)
		print(f"{y1:>5} | {q11:.2f}   {q21:.2f}")
		print(f"{y2:>5} | {q12:.2f}   {q22:.2f}")

def setInjectionfromSOI(durationsArray, selectorMap, iq, rpm, soi, setv):
	durationsReversed = selectorMap.xAxis.copy()
	durationsReversed.reverse() # get selector axis in low to high order (0,5,10,15,....)

	bound = Map.findUpperLowerBoundaries(durationsReversed, soi) # find 2 duration maps associated with that SOI

	fromDur = selectorMap.getValue(bound[0], 0)
	toDur = selectorMap.getValue(bound[1], 0) 

	fromDurValue = durationsArray[fromDur].interpolate(rpm, iq) # in durations, mg is Y, rpm is X, opposite of SOI
	toDurValue = durationsArray[toDur].interpolate(rpm, iq) # in durations, mg is Y, rpm is X, opposite of SOI

	newValues = Map.rawSetInterpolate(soi, 0, bound[0], bound[1], 0, 0, fromDurValue, toDurValue, None, None, setv)

	print(newValues)

	print("Changes needed on map " + str(fromDur) + ':')
	displayNeededChanges(*durationsArray[fromDur].setInterpolate(rpm, iq, newValues[0]))
	print("Changes needed on map " + str(toDur) + ':')
	displayNeededChanges(*durationsArray[toDur].setInterpolate(rpm, iq, newValues[1]))

def findBestSOIToMatchEOI(rpm, iq, soi, maxSOIDeviation, maxSOI, increment, targetEOI):

	currentEOI = soi - findInjectionfromSOI(durations, selector, iq, rpm, soi)

	if currentEOI > targetEOI: # If current eoi is better than targetEOI, do not touch
		return soi

	bestSOI = soi
	bestEOI = 100

	#first = soi - maxSOIDeviation # we start at lowest possible soi
	first = soi # maybe we should not try with lower soi as it could set very low soi because lower soi = lower durations
	while first < soi + maxSOIDeviation and first <= maxSOI:
		val = first - findInjectionfromSOI(durations, selector, iq, rpm, first) # find different EOI for same iq, rpm but different soi
		diff = abs(val - targetEOI)
		if min(diff, bestEOI) == diff:
			bestEOI = diff
			bestSOI = first

		first = first + increment
	
	return bestSOI

def readToFile(file, bytes):
	file = open(os.path.join(os.path.dirname(__file__), file), 'rb')
	i=0
	while i < 524288:
		byte = file.read(2)
		byte_le = struct.unpack('<h', byte)
		bytes.add(byte_le[0])
		i+=2

def dry_air_density(temp_celsius, pressure_pa=101325):
    R = 287.058  # J/(kg·K), specific gas constant for dry air
    temp_kelvin = temp_celsius + 273.15
    density = pressure_pa / (R * temp_kelvin)
    return density

durationsOptions = {
	'map': {
		'multiplier': 0.023437,
		'add': 0
	},
	'x': {
		'multiplier': 1,
		'add': 0
	},
	'y': {
		'multiplier': 0.01,
		'add': 0
	}
}

IQByMapOptions = {
	'map': {
		'multiplier': 0.01,
		'add': 0
	},
	'x': {
		'multiplier': 1,
		'add': 0
	},
	'y': {
		'multiplier': 1,
		'add': 0
	}
}

BoostOptions = {
	'map': {
		'multiplier': 1,
		'add': 0
	},
	'x': {
		'multiplier': 0.01,
		'add': 0
	},
	'y': {
		'multiplier': 1,
		'add': 0
	}
}

SOIOptions = {
	"onlyData": True,
	'map': {
		'multiplier': -0.023437,
		'add': 78
	},
	'x': {
		'address': 0x00000,
		'multiplier': 0.01,
		'add': 0
	},
	'y': {
		'address': 0x00000,
		'multiplier': 1,
		'add': 0
	}
}

selectorOptions = {
	"onlyX": True,
	'x': {
		'func': lambda x : int(round(x)),
		'multiplier': -0.023437,
		'add': 78
	},
	'map':{
		'func': lambda x : int(round(x)),
		'multiplier': 0.003906,
		'add': 0
	}
}

TLOptions = {
	'map': {
		'multiplier': 0.010000,
		'add': 0
	},
	'x': {
		'multiplier': 1,
		'add': 0
	},
	'y': {
		'multiplier': 1,
		'add': 0
	}
}

edc15File = "polo-blt-moffa.bin"
bytes = DoubleLinkedList()
readToFile(edc15File, bytes)

# Codeblock
# Usually 0x40000, 0x50000, 0x60000
baseAddr = 0x60000

selector = Map(bytes, baseAddr + 0x1544c, selectorOptions)
durations = []
durations.append(Map(bytes, baseAddr + 0x146a0, durationsOptions))
durations.append(Map(bytes, baseAddr + 0x14798, durationsOptions))
durations.append(Map(bytes, baseAddr + 0x14a1e, durationsOptions))
durations.append(Map(bytes, baseAddr + 0x14ca4, durationsOptions))
durations.append(Map(bytes, baseAddr + 0x14f2a, durationsOptions))
durations.append(Map(bytes, baseAddr + 0x151b0, durationsOptions))
SOIOptions["x"]["address"] = baseAddr + 0x1887A
SOIOptions["y"]["address"] = baseAddr + 0x18856

IQByMap = Map(bytes, baseAddr + 0x0E3C8, IQByMapOptions)

SOI = Map(bytes, baseAddr + 0x1a15a, SOIOptions)

Boost = Map(bytes, baseAddr + 0x16C64, BoostOptions)

TL = Map(bytes, baseAddr + 0x0D91E, TLOptions)

TLINJ = Map(None, 0)
TLINJ.copyAxis(TL)

TLSOI = Map(None, 0)
TLSOI.copyAxis(TL)

TLEOI = Map(None, 0)
TLEOI.copyAxis(TL)

TLAFR = Map(None, 0)
TLAFR.copyAxis(TL)

TLBoost = Map(None, 0)
TLBoost.copyAxis(TL)

TLTQ = Map(None, 0)
TLTQ.copyAxis(TL)

TLHP = Map(None, 0)
TLHP.copyAxis(TL)

SmokeMargin = Map(None, 0)
SmokeMargin.copyAxis(TL)

SMOKEMAPIATREF = 30
IAT = 30

for x in range(TL.xAxisSize):
	for y in range(TL.yAxisSize):
		rpm = TL.xAxis[x]
		pressure = TL.yAxis[y] # Ambient pressure found in TL y axis
		IQ = TL.getValue(rpm, pressure) # Requested IQ
		SOIv = SOI.interpolate(IQ, rpm) # Associated SOI with rpm/IQ combo
		InjectionQuantity = findInjectionfromSOI(durations, selector, IQ, rpm, SOIv) # Associated Injection ° from SOI, rpm, IQ
		Boostv = Boost.interpolate(IQ, rpm) # Requested boost
		TLINJ.setV(x, y, InjectionQuantity)
		TLEOI.setV(x, y, SOIv - InjectionQuantity)
		TLSOI.setV(x, y, SOIv)
		IAT_density = dry_air_density(IAT, Boostv * 100) # Density according to abs boost pressure and IAT
		VE = Map.rawInterpolate(95.94, 82.5, 900, 5000, rpm) / 100 # VE efficiency based on linear formula
		AFR = 0 if IQ == 0 else (474 * IAT_density * VE) / IQ # Calculated AFR
		
		AdjustedBoostByIAT = ((273.15 + SMOKEMAPIATREF) / (273.15 + IAT)) * Boostv # Calculate how much boost that would be at IAT (used to match smoke map that is calibrated for X°C)
		AllowedSmoke = IQByMap.interpolate(AdjustedBoostByIAT, rpm)
		# You choose PD 150 stock map related
		# Found on 038906019HK
		# For a 150, 57.5mg@1900 is 320Nm
		# 54mg @ 4000 is 150hp (263Nm)
		#TQ_coeff = 0 if rpm < 1900 or rpm > 4000 else Map.rawInterpolate(1, 0.876, 1900, 4000, rpm)
		#TQ = (320 / 57.5) * IQ * TQ_coeff
		###
		# You choose PD 130 stock map related
		# Found on 038906019NJ
		# For a 130, 59.5mg@1900 is 310Nm
		# 47,9mg @ 4000 is 130hp (228Nm)
		TQ_coeff = 0 if rpm < 1900 or rpm > 4000 else Map.rawInterpolate(1, 0.914, 1900, 4000, rpm)
		TQ = (310 / 59.5) * IQ * TQ_coeff
		###
		HP = (TQ * rpm) / 7022
		TLTQ.setV(x, y, TQ)
		TLHP.setV(x, y, HP)
		TLAFR.setV(x, y, AFR)

		# This map is showing the margin of smoke before reaching smoke map limits
		# a negative delta IQ is IQ that will never be reached due to insufficient boost for requested afr at IAT
		# a positive is IQ that can be gained just by raising TL.
		SmokeMargin.setV(x, y, AllowedSmoke - IQ)

		TLBoost.setV(x, y, Boostv - pressure)


### Experimental
SOIv = SOI.interpolate(76, 3100)
print(findInjectionfromSOI(durations, selector, 76, 3100, SOIv))
setInjectionfromSOI(durations, selector, 76, 3100, SOIv, 30)



print("TL Map")
print(TL)
print("TL Inj ° Map")
print(TLINJ)
print("TL SOI ° Map")
print(TLSOI)
print("TL EOI ° BTDC Map")
print(TLEOI)
print("TL AFR IAT=" + str(IAT) + "°C" )
print(TLAFR)
print("TL Torque")
print(TLTQ)
print("TL HP")
print(TLHP)
print("Removed TQ due to IAT (negative is pulled mg)")
print(SmokeMargin)
print("TL Boost relative Map")
print(TLBoost)

Injection = Map(None, 0)
Injection.copyAxis(SOI)

EOI = Map(None, 0)
EOI.copyAxis(SOI)


for x in range(SOI.xAxisSize):
	for y in range(SOI.yAxisSize):
		realSOIX = SOI.xAxis[x] # mg value
		realSOIY = SOI.yAxis[y] # rpm
		value = SOI.get(x, y) # value

		
		'''if realSOIX >= 45 and realSOIY >= 2000:
			bestSOI = findBestSOIToMatchEOI(realSOIY, realSOIX, value, 15, 27, 0.005, -8)
			value = bestSOI
			SOI.setV(x, y, bestSOI)'''

		realInjectionVal = findInjectionfromSOI(durations, selector, realSOIX, realSOIY, value)
		EOI.setV(x, y, value - realInjectionVal)
		Injection.setV(x, y, realInjectionVal)

print("SOI TABLE")
print(SOI)

print("INJECTION TABLE")
print(Injection)

print("EOI TABLE")
print(EOI)

#SOI.writeToFile(edc15File)