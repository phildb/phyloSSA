from sortedcontainers import SortedListWithKey, SortedList, SortedSet
from random import sample, randint, uniform, random, seed, getstate
from itertools import count, permutations, combinations, izip, chain
from math import sqrt, log, exp, pi
from collections import namedtuple, defaultdict
from sys import exit
from time import time

Point = namedtuple('Point', ['niche','age', 'vulnerability'])

fixed_seed = False
pretty = False

global nicheWidth
nicheWidth = 21

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

	def __repr__(self):
		return '\nBirth (n=' + repr(self.point.niche) + ', a=' + repr(self.point.age) + ', v=' + repr(self.point.vulnerability) + ') | rate=' + repr(self.rate)

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

	def __repr__(self):
		return '\nDeath (n=' + repr(self.point.niche) + ', a=' + repr(self.point.age) + ', v=' + repr(self.point.vulnerability) + ') | rate=' + repr(self.rate)

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
		self.rate = 0.04 * self.lattice.n(self.point)

	def do(self):
		new_mutant = Point(niche=(self.point.niche + sample([-1,1],1)[0] + (random() - 0.5)) % nicheWidth, age=self.lattice.t())
		#new_mutant = Point(niche=(self.point.niche + sample([-1,1],1)[0]) % nicheWidth, age=self.lattice.t(), vulnerability=0)
		self.lattice.create(new_mutant)
		return [self.point, new_mutant]

class StepMutationProcess(object):

	global nicheWidth

	def __init__(self, point, lattice):
		self.order = 1
		self.rate = 0.0
		self.point = point
		self.lattice = lattice

	def __del__(self):
		self.point = None
		self.lattice = None
	
	def __repr__(self):
		return '\nStepMut (n=' + repr(self.point.niche) + ', a=' + repr(self.point.age) + ', v=' + repr(self.point.vulnerability) + ') | rate=' + repr(self.rate)

	def refresh_rate(self):
		self.rate = 0.04 * self.lattice.n(self.point)

	def do(self):
		niche_step = sample([-2,-1,-1,0,1,1,2],1)[0]
		vulnerability_step = sample([-1,0,1],1)[0]
		new_mutant = Point(niche=(self.point.niche + niche_step) % nicheWidth, age=self.lattice.t(), vulnerability=0)
		#new_mutant = Point(niche=(self.point.niche + sample([-1,1],1)[0]) % nicheWidth, age=self.lattice.t())
		self.lattice.create(new_mutant)
		return ('mutation', self.point, new_mutant)


class IntraCompetitionProcess(object):
	def __init__(self, point, lattice):
		self.order = 1 # Doesn't interact across niches or ages, so order 1. This is slightly inconsistent, should find a way to fix.
		self.rate = 0.0
		self.point = point
		self.lattice = lattice

	def __del__(self):
		self.point = None
		self.lattice = None

	def __repr__(self):
		return '\nIC (n=' + repr(self.point.niche) + ', a=' + repr(self.point.age) + ', v=' + repr(self.point.vulnerability) + ') | rate=' + repr(self.rate)

	def refresh_rate(self):
		n = self.lattice.n(self.point)
		self.rate = 0.03 * n * (n-1)

	def do(self):
		self.lattice.annihilate(self.point)
		return ('intracomp', self.point)

class NNInterCompetitionProcess(object):

	global nicheWidth

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

	def __repr__(self):
		return '\nNNC (n=' + repr(self.point1.niche) + ', a=' + repr(self.point1.age) + ', v=' + repr(self.point1.vulnerability) +  ') <- (n=' + repr(self.point2.niche) + ', a=' + repr(self.point2.age) + ', v=' + repr(self.point2.vulnerability) +  ') | rate=' + repr(self.rate)

	# Starting from point, generate all NN interactions
	# consistent with creating a new individual at 'point'
	@staticmethod
	def stencile(point, lattice):
		# You MUST walk on everybody because all ages and other dimensions
		# aside from niche might be impicated in a NN competitive interaction
		for p2 in lattice.sites:
			if ((point.niche + 1) % nicheWidth == p2.niche or (point.niche - 1) % nicheWidth == p2.niche) or (point.niche == p2.niche and point != p2):
				yield (point, p2)
				yield (p2, point)


	def refresh_rate(self):
		self.rate = 0.03 * self.lattice.n(self.point1) * self.lattice.n(self.point2)

	def do(self):
		self.lattice.annihilate(self.point1)
		return ('intercomp', self.point2, self.point1)

# class RedQueenProcess(object):
# 	def __init__(self, point1, point2, lattice):
# 		self.order = 2
# 		self.rate = 0.0
# 		self.point1 = point1
# 		self.point2 = point2
# 		self.lattice = lattice

# 	def __del__(self):
# 		self.point1 = None
# 		self.point2 = None
# 		self.lattice = None

# 	# @staticmethod
# 	# def stencile(point, lattice):
# 	# 	left_nn = Point(niche=point.niche-1,)

# 	def refresh_rate(self):
# 		#self.rate = 0.01 * self.lattice.n(self.point1) * self.lattice.n(self.point2)
# 		self.rate = (0.01 + cmp()

# 	def do(self):
# 		self.lattice.annihilate(self.point1)


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
			self.meta.pull_processes(affected)
			self.sites[point] += 1
			for pr in affected:
				pr.refresh_rate()
			self.meta.push_processes([a for a in affected if a.rate > 0.0])
		else:
			self.sites[point] = 1
			self.attached_processes[point] = []
			for ptemplate in self.classes1st:
				process = ptemplate(point, self)
				self.attached_processes[point].append(process)

			for ptemplate in self.classes2nd:
				for p1, p2 in ptemplate.stencile(point, self):
					process = ptemplate(p1, p2, self)
					self.attached_processes[p1].append(process)
					self.attached_processes[p2].append(process)

			for pr in self.attached_processes[point]:
				pr.refresh_rate()

			self.meta.push_processes([a for a in self.attached_processes[point] if a.rate > 0.0])

	def annihilate(self, point):
		self.total_n -= 1
		if point in self.sites.keys():
			affected = self.attached_processes[point]
			self.meta.pull_processes(affected)
			self.sites[point] -= 1
			for pr in affected:
				pr.refresh_rate()
			if self.sites[point] > 0:
				self.meta.push_processes([a for a in affected if a.rate > 0.0])
			else:
				del self.sites[point]
				for pr in self.attached_processes[point]:
					if pr.order == 2:
						if pr.point1 != point:
							self.attached_processes[pr.point1].remove(pr)
						if pr.point2 != point:
							self.attached_processes[pr.point2].remove(pr)
				del self.attached_processes[point]
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
		self._sorted_processes = SortedSet(key=lambda pr: pr.rate)
		self._process_set = set()
		self._total_rate = 0.0
		self._num_processes = 0
		self._t = 0.0
		self.tries = 0
		self.processed = 0

	def pull_processes(self, processes):
		nonnull_processes = [a for a in processes if a.rate > 0.0]
		self._sorted_processes -= nonnull_processes
		#for pr in nonnull_processes:
		#	self._sorted_processes.remove(pr)
		self._total_rate -= sum((pr.rate for pr in nonnull_processes))
		self._num_processes -= len(nonnull_processes)

	def push_processes(self, processes):
		nonnull_processes = [a for a in processes if a.rate > 0.0]
		self._total_rate += sum(a.rate for a in nonnull_processes)
		self._num_processes += len(nonnull_processes)
		self._sorted_processes |= nonnull_processes
		#self._sorted_processes.update(nonnull_processes)

	# This might not work because the list of process is sorted
	def step_gillespie(self):
		total_rate = sum(pr.rate for pr in self._sorted_processes)
		if total_rate > 0.0:
			r1 = random()
			self._t += -log(r1)/total_rate
			r2 = random()
			acc = 0.0
			self.processed += 1
			for pr in self._sorted_processes:
				acc += pr.rate
				if r2*total_rate < acc:
					self.tries += 1
					pr.do()
					break
		else:
			print('Party\'s over, everybody\'s dead')
			exit()


	def step_purerejection(self):
		residue = None
		#if self._total_rate > 0.0:
		#if True:
		# This one's the right one, the true absorbing state; no relying on floats:
		if len(meta._sorted_processes) > 0:
			r1 = random()
			delta_t = -log(r1)/self._total_rate
			self._t += delta_t

			# For now processes of rate 0.0 might creep in and
			# Slepoy's et al algorithm cannot work without picking
			# the first nonzero rate.
			# min_rate = self._sorted_processes[0].rate
			max_rate = self._sorted_processes[-1].rate

			for i in count():
				#number_of_processes = len(self._sorted_processes)
				#ith_pr = randint(0,number_of_processes-1)
				#process = self._sorted_processes[ith_pr]
				# Test with fixed seed indicates that sampling directly _processes
				# gives the exact same process as sampling an integer 0 <= i < len(_processes)
				# and picking _processes[i]
				
				
				#print('wtf ' + repr(self._num_processes) + ' | should be ' + repr(len(self._sorted_processes)))
				#if self._num_processes != len(self._sorted_processes):
				#	print 'wtf'
				#process = sample(self._sorted_processes,1)[0]
				
				ri = randint(0, len(self._sorted_processes)-1)
				process = self._sorted_processes[ri]

				if process.rate >= uniform(0, max_rate):
					residue = process.do()
					self.processed += 1
					self.tries += i+1
					break
		else:
			print('Party\'s over, everybody\'s dead')
			exit()

		return residue



dark_palette = [[180,220,212],[219,65,42],[203,82,223],[104,220,73],[55,32,58],[221,173,53],[83,132,209],[222,132,172],[53,88,51],[98,225,160],[222,182,145],[73,126,144],[144,49,73],[205,220,64],[137,44,115],[146,63,32],[123,113,230],[223,51,97],[108,165,56],[146,103,82],[210,139,223],[53,43,29],[206,202,117],[217,124,57],[159,132,58],[87,160,105],[47,70,114],[212,72,177],[122,181,217],[144,155,128],[170,227,128],[195,223,173],[220,200,212],[86,28,28],[42,71,71],[94,79,27],[217,131,120],[175,162,223],[151,99,154],[78,116,40],[83,162,151],[100,223,217],[214,64,131],[114,82,172],[118,68,91],[78,49,107],[109,106,100],[181,145,165],[206,79,85],[110,113,142]]
light_palette = [[228,222,154],[224,177,237],[79,225,230],[248,150,144],[99,245,176],[248,167,82],[242,230,81],[149,205,107],[191,195,213],[173,238,199],[234,174,134],[227,205,106],[154,183,232],[215,247,102],[239,159,200],[50,243,221],[234,171,173],[100,211,241],[103,209,189],[246,194,79],[88,203,150],[163,234,233],[123,230,139],[178,195,116],[230,190,213],[151,196,156],[226,177,98],[239,156,102],[191,239,151],[142,203,133],[187,208,99],[188,228,169],[166,240,127],[188,174,211],[234,249,163],[172,215,194],[121,197,198],[229,242,122],[138,236,180],[203,181,106],[216,208,242],[171,208,222],[119,242,210],[240,221,98],[68,217,177],[218,190,75],[201,188,128],[230,213,131],[195,210,81],[212,227,145]]

iwanthue = dict()

if fixed_seed: seed(0.652134937)

meta = Metaprocess()
lattice = Lattice(meta, [], [BirthProcess, DeathProcess, IntraCompetitionProcess, StepMutationProcess], [NNInterCompetitionProcess])
#lattice = Lattice(meta, [], [BirthProcess, DeathProcess, IntraCompetitionProcess], [])


lattice.create(Point(niche=5,age=0.0,vulnerability=0))
lattice.create(Point(niche=5,age=0.0,vulnerability=0))
lattice.create(Point(niche=5,age=0.0,vulnerability=0))
lattice.create(Point(niche=5,age=0.0,vulnerability=0))
lattice.create(Point(niche=5,age=0.0,vulnerability=0))
lattice.create(Point(niche=5,age=0.0,vulnerability=4))
lattice.create(Point(niche=5,age=0.0,vulnerability=4))
lattice.create(Point(niche=5,age=0.0,vulnerability=4))
lattice.create(Point(niche=5,age=0.0,vulnerability=4))
lattice.create(Point(niche=5,age=0.0,vulnerability=4))
#lattice.annihilate(Point(niche=5,age=0.0,vulnerability=0))
print(meta._sorted_processes)
print(lattice.attached_processes)
#exit()

# print('---')
# for i in xrange(0,21):
# 	lattice.create(Point(niche=i%nicheWidth, age=0.0, vulnerability=0))




# for i in xrange(1):
# 	#initial_niche = uniform(0.0,10.0)
# 	#initial_niche = randint(0,nicheWidth)
# 	initial_niche = 5.0
# 	lattice.create(Point(niche=initial_niche, age=0.0))
# for i in xrange(20):
# 	lattice.create(Point(niche=i, age=0.0))

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

		for i in xrange(10):
			residue = meta.step_purerejection()
			if meta.processed % 1000 == 0:
				print('processed = ' + repr(meta.processed) + ' | total pop = ' + repr(lattice.total_n) + ' | num. sps. = ' + repr(len(lattice.sites)) + ' | active proc. = ' + str(len(meta._sorted_processes))) + ' | tries in 1K = ' + repr(meta.tries) + ' | avg. tries = ' + '%.1f'%(float(meta.tries)/1000) + ' | min rate = ' + '%.2f'%(meta._sorted_processes[0].rate) + ' | max rate = ' + '%.2f'%(meta._sorted_processes[-1].rate)
				meta.tries = 0

		#residue = meta.step_purerejection()
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
			if point in iwanthue.keys():
				fill(*iwanthue[point])
			else:
				new_hue = sample(light_palette,1)[0]
				iwanthue[point] = new_hue
				fill(*iwanthue[point])

			#r = 5*sqrt(n)
			r = n
			#ellipse(50 + point.niche/10.0*700.0, 50 + y ,r,r)
			rect(50 + float(point.niche)/nicheWidth*700.0 + 0*point.age, 50 + y + point.age, r, 1)

		if residue != None:
			if residue[0] == 'mutation':
				stroke(102,210,136)
				fill(102,210,136)
				niche_from = 50 + float(residue[1].niche)/nicheWidth*700.0 + 0*residue[1].age
				niche_mutant = 50 + float(residue[2].niche)/nicheWidth*700.0 + 0*residue[2].age
				line(niche_from, 50 + y + residue[1].age, niche_mutant, 50 + y + residue[2].age)
				ellipse(niche_mutant, 50 + y + residue[2].age, 7, 7)
			elif residue[0] == 'intracomp':
				stroke(193,79,77)
				fill(193,79,77)
				ellipse(50 + float(residue[1].niche)/nicheWidth*700.0 + 0*residue[1].age, 50 + y + residue[1].age, 7, 7)
			elif residue[0] == 'intercomp':
				stroke(193,79,77)
				fill(193,79,77)
				inter_from = 50 + float(residue[1].niche)/nicheWidth*700.0 + 0*residue[1].age
				inter_to = 50 + float(residue[2].niche)/nicheWidth*700.0 + 0*residue[2].age
				line(inter_from, 50 + y + residue[1].age, inter_to, 50 + y + residue[2].age)
				ellipse(inter_to, 50 + y + residue[2].age, 7, 7)
		# Draw a red dot in the middle of the population
		# for i,(point, n) in enumerate(lattice.sites.iteritems()):
		# 	y = (meta._t/0.01) % 700
		# 	fill(255,0,0)
		# 	rect(50 + point.niche/10.0*700.0, 50 + y, 1, 1)

	run()
else:
	old_time = time()
	for i in count():
		meta.step_purerejection()
		#meta.step_gillespie()
		if meta.processed % 1000 == 0:
			new_time = time()
			print('processed = ' + repr(meta.processed) + ' | total pop = ' + repr(lattice.total_n) + ' | num. sps. = ' + repr(len(lattice.sites)) + ' | active proc. = ' + str(len(meta._sorted_processes)) + ' | tries in 1K = ' + repr(meta.tries) + ' | avg. tries = ' + '%.1f'%(float(meta.tries)/1000) + ' | min rate = ' + '%.2f'%(meta._sorted_processes[0].rate) + ' | max rate = ' + '%.2f'%meta._sorted_processes[-1].rate + ' | delta total = ' + '%.2e'%(meta._total_rate - sum(pr.rate for pr in meta._sorted_processes)) + ' | in %.2f'%(new_time-old_time) + 'secs')
			old_time = new_time
			meta.tries = 0
		if meta.processed % 100000 == 0:
			meta._total_rate = sum(pr.rate for pr in meta._sorted_processes)
