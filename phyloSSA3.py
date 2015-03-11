from random import sample, randint, uniform, random, seed, getstate
from itertools import count, permutations, combinations, izip, chain
from math import sqrt, log, exp, pi
from collections import namedtuple, defaultdict
from sys import exit

Point = namedtuple('Point', ['niche','age'])

fixed_seed = True
pretty = True

if fixed_seed: seed(0.652134937)


global nicheWidth
nicheWidth = 11

class BirthProcess(object):
	''' The Birth process is a 1st order process sitting at a Point.\n
	It create 1 new particle at Point'''
	def __init__(self, point, lattice):
		self.order = 1
		self.rate = 0.0
		self.point = point
		self.lattice = lattice

	def __del__(self):
		self.point = None
		self.lattice = None

	# def __repr__(self):
	# 	return '(birth ' + repr(self.point) + ' rate=' + repr(self.rate) + ')'

	def refresh_rate(self):
		self.rate = 1.0 * self.lattice.n(self.point)

	def do(self):
		self.lattice.create(self.point)


class DeathProcess(object):
	def __init__(self, point, lattice):
		self.order = 1
		self.rate = 0.0
		self.point = point
		self.lattice = lattice

	def __del__(self):
		self.point = None
		self.lattice = None

	# def __repr__(self):
	# 	return '(death ' + repr(self.point) + ' rate=' + repr(self.rate) + ')'
	
	def refresh_rate(self):
		self.rate = 0.5 * self.lattice.n(self.point)

	def do(self):
		self.lattice.annihilate(self.point)


class MutationProcess(object):

	global nicheWidth

	def __init__(self, point, lattice):
		self.order = 1
		self.rate = 0.0
		self.point = point
		self.lattice = lattice

	def __del__(self):
		self.point = None
		self.lattice = None
	
	def refresh_rate(self):
		self.rate = 0.03 * self.lattice.n(self.point)

	def do(self):
		new_mutant = Point(niche=(self.point.niche + sample([-1,1],1)[0] + (random() - 0.5)) % 10, age=self.lattice.t())
		#new_mutant = Point(niche=(self.point.niche + sample([-1,1],1)[0]) % nicheWidth, age=self.lattice.t())
		self.lattice.create(new_mutant)
		return [self.point, new_mutant]


class IntraCompetitionProcess(object):
	def __init__(self, point, lattice):
		self.order = 1 # Doesn't interact across niches or ages, so order 1. This is slightly inconsistent, should find a way to fix.
		self.rate = 0.0
		self.point = point
		self.lattice = lattice

	def __del__(self):
		self.point = None
		self.lattice = None

	# def __repr__(self):
	# 	return '(intracompatition in ' + repr(self.point) + ' rate=' + repr(self.rate) + ')'

	def refresh_rate(self):
		n = self.lattice.n(self.point)
		self.rate = 0.03 * n * (n-1)

	def do(self):
		self.lattice.annihilate(self.point)

class NNInterCompetitionProcess(object):
	def __init__(self, point1, point2, lattice):
		self.order = 2
		self.rate = 0.0
		self.point1 = point1
		self.point2 = point2
		self.lattice = lattice

	def __del__(self):
		self.point1 = None
		self.point2 = None
		self.lattice = None

	# @staticmethod
	# def stencile(point, lattice):
	# 	left_nn = Point(niche=point.niche-1,)

	def refresh_rate(self):
		self.rate = 0.03 * self.lattice.n(self.point1) * self.lattice.n(self.point2)

	def do(self):
		self.lattice.annihilate(self.point1)

class Bisector(object):
	def __init__(self,rate):
		self.rate = rate

class Lattice(object):
	'''The Lattice is the interface and the glue between the Metaprocess containing\n
	all active processes and the internal degrees of freedom. All creation-annihilation must go
	through the Lattice so that it controls the affect of processes that touch
	at a given site.'''
	def __init__(self, meta, zerothorder, firstorder, secondorder):
		self.sites = dict()
		self.attached_processes = dict()
		self.meta = meta
		self.classes0th = zerothorder
		self.classes1st = firstorder
		self.classes2nd = secondorder
		self.total_n = 0

	def create(self, point):
		self.total_n += 1
		if point in self.sites.keys():
			affected = self.attached_processes[point]
			self.sites[point] += 1
			for pr in affected:
				pr.refresh_rate()
		else:
			self.sites[point] = 1
			self.attached_processes[point] = []
			for ptemplate in self.classes1st:
				process = ptemplate(point, self)
				self.attached_processes[point].append(process)

			for ptemplate in self.classes2nd:
				for point2 in self.sites:
					if point2 != point:
						process1 = ptemplate(point, point2, self)
						process2 = ptemplate(point2, point, self)
						self.attached_processes[point].append(process1)
						self.attached_processes[point].append(process2)
						self.attached_processes[point2].append(process2)
						self.attached_processes[point2].append(process1)

			for pr in self.attached_processes[point]:
				pr.refresh_rate()

			self.meta._process_set |= set(self.attached_processes[point])

	def annihilate(self, point):
		self.total_n -= 1
		if point in self.sites.keys():
			affected = self.attached_processes[point]
			self.sites[point] -= 1
			for pr in affected:
				pr.refresh_rate()
			if self.sites[point] <= 0:
				del self.sites[point]
				for pr in self.attached_processes[point]:
					if pr.order == 2:
						if pr.point1 != point:
							self.attached_processes[pr.point1].remove(pr)
						if pr.point2 != point:
							self.attached_processes[pr.point2].remove(pr)
				del self.attached_processes[point]
				self.meta._process_set -= set(affected)
		else:
			raise ValueError('Orphan process, annihilation at empty site ' + repr(point))

	def n(self, point):
		if point in self.sites.keys():
			return self.sites[point]
		else:
			return 0

	def t(self):
		return meta._t




class Metaprocess(object):
	def __init__(self):
		self._process_set = set()
		self._total_rate = 0.0
		self._num_processes = 0
		self._t = 0.0
		self.tries = 0
		self.processed = 0

	def pull_processes(self, processes):
		nonnull_processes = [a for a in processes if a.rate > 0.0]
		for pr in nonnull_processes:
			self._sorted_processes.remove(pr)
		self._total_rate -= sum((pr.rate for pr in nonnull_processes))
		self._num_processes -= len(nonnull_processes)

	def push_processes(self, processes):
		nonnull_processes = [a for a in processes if a.rate > 0.0]
		self._total_rate += sum(a.rate for a in nonnull_processes)
		self._num_processes += len(nonnull_processes)
		self._sorted_processes.update(nonnull_processes)

	# This might not work because the list of process is sorted
	def step_gillespie(self):
		total_rate = sum(pr.rate for pr in self._process_set)
		if total_rate > 0.0:
			r1 = random()
			self._t += -log(r1)/total_rate
			r2 = random()
			acc = 0.0
			self.processed += 1
			for pr in self._process_set:
				acc += pr.rate
				if r2*total_rate < acc:
					self.tries += 1
					residue = pr.do()
					break
		else:
			print('Party\'s over, everybody\'s dead')
			exit()

		return residue



dark_palette = [[180,220,212],[219,65,42],[203,82,223],[104,220,73],[55,32,58],[221,173,53],[83,132,209],[222,132,172],[53,88,51],[98,225,160],[222,182,145],[73,126,144],[144,49,73],[205,220,64],[137,44,115],[146,63,32],[123,113,230],[223,51,97],[108,165,56],[146,103,82],[210,139,223],[53,43,29],[206,202,117],[217,124,57],[159,132,58],[87,160,105],[47,70,114],[212,72,177],[122,181,217],[144,155,128],[170,227,128],[195,223,173],[220,200,212],[86,28,28],[42,71,71],[94,79,27],[217,131,120],[175,162,223],[151,99,154],[78,116,40],[83,162,151],[100,223,217],[214,64,131],[114,82,172],[118,68,91],[78,49,107],[109,106,100],[181,145,165],[206,79,85],[110,113,142]]
light_palette = [[228,222,154],[224,177,237],[79,225,230],[248,150,144],[99,245,176],[248,167,82],[242,230,81],[149,205,107],[191,195,213],[173,238,199],[234,174,134],[227,205,106],[154,183,232],[215,247,102],[239,159,200],[50,243,221],[234,171,173],[100,211,241],[103,209,189],[246,194,79],[88,203,150],[163,234,233],[123,230,139],[178,195,116],[230,190,213],[151,196,156],[226,177,98],[239,156,102],[191,239,151],[142,203,133],[187,208,99],[188,228,169],[166,240,127],[188,174,211],[234,249,163],[172,215,194],[121,197,198],[229,242,122],[138,236,180],[203,181,106],[216,208,242],[171,208,222],[119,242,210],[240,221,98],[68,217,177],[218,190,75],[201,188,128],[230,213,131],[195,210,81],[212,227,145]]

iwanthue = dict()


meta = Metaprocess()
lattice = Lattice(meta, [], [BirthProcess, DeathProcess, IntraCompetitionProcess, MutationProcess], [NNInterCompetitionProcess])
#lattice = Lattice(meta, [], [BirthProcess, DeathProcess, IntraCompetitionProcess], [])

# for i in xrange(1):
# 	#initial_niche = uniform(0.0,10.0)
# 	#initial_niche = randint(0,nicheWidth)
# 	initial_niche = 5.0
# 	lattice.create(Point(niche=initial_niche, age=0.0))
for i in xrange(10):
	lattice.create(Point(niche=i, age=0.0))

print(lattice.attached_processes)

#exit()
if pretty:
	from pyprocessing import *

	def setup():
		frameRate(200)
		size(800,800)
		ellipseMode(CENTER)
		noStroke()
		fill(0)
		rect(0,0,800,800)

	def draw():
		residue = meta.step_gillespie()
		#meta.step_gillespie()
		rectMode(CORNER)
		noStroke()
		fill(0)
		rect(0,0,800,50)
		fill(200)
		textSize(15)
		text(' t = %.5f' % meta._t 
			+ ' n = ' + repr(lattice.total_n) 
			+ ' av. tries = ' + '%.1f'%(float(meta.tries)/float(meta.processed))
			+ ' processed = ' + repr(meta.processed)
			,25,25)
		fill(0,5)
		rect(0,0,800,800)
		rectMode(RADIUS)
		y = (meta._t/0.01) % 700
		for i,(point, n) in enumerate(lattice.sites.iteritems()):
			if point.niche in iwanthue.keys():
				fill(*iwanthue[point.niche])
			else:
				new_hue = sample(light_palette,1)[0]
				iwanthue[point.niche] = new_hue
				fill(*iwanthue[point.niche])

			#r = 5*sqrt(n)
			r = n
			#ellipse(50 + point.niche/10.0*700.0, 50 + y ,r,r)
			rect(50 + float(point.niche)/nicheWidth*700.0, 50 + y ,r,1)

		if residue != None:
			stroke(102,210,136)
			fill(102,210,136)

			niche_from = 50 + float(residue[0].niche)/nicheWidth*700.0
			niche_mutant = 50 + float(residue[1].niche)/nicheWidth*700.0
			line(niche_from, 50 + y, niche_mutant, 50 + y)
			ellipse(niche_mutant, 50 + y, 7, 7)
		
		# Draw a red dot in the middle of the population
		# for i,(point, n) in enumerate(lattice.sites.iteritems()):
		# 	y = (meta._t/0.01) % 700
		# 	fill(255,0,0)
		# 	rect(50 + point.niche/10.0*700.0, 50 + y, 1, 1)

	run()
else:
	for i in count():
		meta.step_gillespie()
		if meta.processed % 1000 == 0:
			print('processed = ' + repr(meta.processed) + ' | total pop = ' + repr(lattice.total_n) + ' | num. sps. = ' + repr(len(lattice.sites)) + ' | active proc. = ' + str(len(meta._process_set)) + ' | tries in 1K = ' + repr(meta.tries) + ' | avg. tries = ' + '%.1f'%(float(meta.tries)/1000) + ' | min rate = ' + '%.2f'%min(pr.rate for pr in meta._process_set) + ' | max rate = ' + '%.2f'%max(pr.rate for pr in meta._process_set))
			meta.tries = 0
