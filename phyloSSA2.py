from sortedcontainers import SortedListWithKey, SortedList, SortedSet
from random import sample, randint, uniform, random, seed, getstate
from itertools import count, permutations, combinations, izip, chain
from math import sqrt, log, exp, pi, floor, tanh, sin, cos
from collections import namedtuple, defaultdict
from sys import exit
#import ete2
import time


# You have to be careful here, object can proliferate when calling them
# inline as parameters, creating new dictionary keys. Tuples both avoid
# the problem but create new ones by being immutable.

Point = namedtuple('Point', ['niche','age', 'vulnerability'])

# class Point(object):
#  	def __init__(self, niche, age, vulnerability):
#  		self.niche = niche
#  		self.age = age
#  		self.vulnerability = vulnerability



# class ProcessBin(object):
# 	def __init__(self, binrate=0.0, processes_in_bin=SortedSet(key=lambda pr: pr.rate)):
# 		self.binrate = binrate
# 		self.processes_in_bin = processes_in_bin
# 	def __del__(self):
# 		self.processes_in_bin = None

class ProcessBin(object):
	def __init__(self, binrate=0.0, processes_in_bin={}):
		self.binrate = binrate
		self.processes_in_bin = processes_in_bin
	def __del__(self):
		self.processes_in_bin = None


class Parameters(object):
	''' Contains program parameters'''
	def __init__(self):
		# General parameters
		self.fixed_seed = False
		
		# visualization = 'console', 'pretty', or 'pretty3d'
		#self.visualization = 'console'
		self.visualization = 'pretty3d'
		
		self.pptrait = 10.0

		self.drawinterval = 100
		self.summaryinterval = 1000
		self.recalibrateinterval = 200000

		# Step algorithm is either 'gillespie', 'rejection', or 'composition-rejection'
		#self.step_algorithm = 'gillespie'
		#self.step_algorithm = 'rejection'
		self.step_algorithm = 'composition-rejection'

		# Stochastic rates and lattice parameters
		self.birthrate = 1.0
		self.deathrate = 0.5
		self.intracomprate = 0.04
		self.nichestep = 1
		self.intercomprate = 0.005

		self.nichemutationrate = 0.03
		self.nicheKernel = [-1, 1]
		self.nicheWidth = 21

		self.vulnerabilitymutationrate = 0.07
		#self.vulnstep = 1
		self.vulnKernel = [-1, 1]
		self.vulnerabilityWidth = 21

		self.catastropherate = 0.005

		# Probabily about smallest rate above, unless some
		# processes have short-tailed kernels, in which case you
		# might want to make this dynamic
		self.pmin = 0.01

		# Pretty colors for plotting
		# Thanks IWantHue @ http://tools.medialab.sciences-po.fr/iwanthue/
		self.dark_palette = [[180,220,212],[219,65,42],[203,82,223],[104,220,73],[55,32,58],[221,173,53],[83,132,209],[222,132,172],[53,88,51],[98,225,160],[222,182,145],[73,126,144],[144,49,73],[205,220,64],[137,44,115],[146,63,32],[123,113,230],[223,51,97],[108,165,56],[146,103,82],[210,139,223],[53,43,29],[206,202,117],[217,124,57],[159,132,58],[87,160,105],[47,70,114],[212,72,177],[122,181,217],[144,155,128],[170,227,128],[195,223,173],[220,200,212],[86,28,28],[42,71,71],[94,79,27],[217,131,120],[175,162,223],[151,99,154],[78,116,40],[83,162,151],[100,223,217],[214,64,131],[114,82,172],[118,68,91],[78,49,107],[109,106,100],[181,145,165],[206,79,85],[110,113,142]]
		self.light_palette = [[228,222,154],[224,177,237],[79,225,230],[248,150,144],[99,245,176],[248,167,82],[242,230,81],[149,205,107],[191,195,213],[173,238,199],[234,174,134],[227,205,106],[154,183,232],[215,247,102],[239,159,200],[50,243,221],[234,171,173],[100,211,241],[103,209,189],[246,194,79],[88,203,150],[163,234,233],[123,230,139],[178,195,116],[230,190,213],[151,196,156],[226,177,98],[239,156,102],[191,239,151],[142,203,133],[187,208,99],[188,228,169],[166,240,127],[188,174,211],[234,249,163],[172,215,194],[121,197,198],[229,242,122],[138,236,180],[203,181,106],[216,208,242],[171,208,222],[119,242,210],[240,221,98],[68,217,177],[218,190,75],[201,188,128],[230,213,131],[195,210,81],[212,227,145]]
		self.indigo_palette = [[68,111,212],[213,85,245],[164,97,166],[231,135,239],[139,102,238],[111,109,177],[168,133,223],[162,81,201],[119,94,205],[212,115,245],[98,117,239],[134,144,233],[161,90,185],[134,100,183],[79,141,241],[153,124,235],[178,73,218],[192,133,238],[203,95,216],[207,125,215],[97,114,195]]
		self.purple_wine_palette = [[110,20,58],[177,36,151],[89,81,70],[192,35,110],[151,71,132],[65,24,32],[98,19,80],[147,74,96],[129,93,93],[146,25,90],[133,32,116],[92,37,59],[122,51,94],[76,21,62],[147,50,80],[112,63,73],[61,14,40],[169,49,132],[150,75,117],[123,47,109],[172,57,110]]
		self.blue_ocean_palette = [[85,126,191],[31,38,50],[115,115,121],[36,51,100],[86,121,230],[78,86,170],[50,100,125],[113,120,157],[47,97,146],[84,78,114],[21,50,80],[63,130,170],[65,68,77],[72,114,196],[79,91,147],[39,39,75],[71,85,110],[34,87,164],[56,127,225],[100,111,203],[32,66,118]]
		
		# Credits to Miaka @ COLOURlovers.com
		self.influenza_palette = [[48,0,48],[72,0,72],[96,24,72],[192,72,72],[240,114,65]]
		
		# Credits to mandarina81 @ COLOURlovers.com
		self.scienceNOVA_palette = [[255,246,154],[245,194,117],[207,125,127],[110,51,88],[39,5,60]]
		
		# Credits to am-y @ COLOURlovers.com
		self.muted_metals_palette = [[17,28,34],[36,61,66],[205,197,127],[237,220,204],[184,179,173]]

		# Set the palette of your choosing
		self.palette = self.influenza_palette

class BirthProcess(object):
	''' The Birth process is a 1st order process sitting at a Point.\n
	It create 1 new particle at Point'''
	def __init__(self, point, lattice, params):
		self.order = 1
		self.rate = 0.0
		self.point = point
		self.lattice = lattice
		self.params = params

	def __del__(self):
		self.point = None
		self.lattice = None
		self.params = None

	def __repr__(self):
		return('\nBirth (n=' + repr(self.point.niche) 
					+ ', a=' + repr(self.point.age) 
					+ ', v=' + repr(self.point.vulnerability) 
			   + ') | rate=' + repr(self.rate))

	def refresh_rate(self):
		self.rate = self.params.birthrate * self.lattice.n(self.point)

	def do(self):
		self.lattice.create(self.point, parent_age=self.point.age)


class DeathProcess(object):
	def __init__(self, point, lattice, params):
		self.order = 1
		self.rate = 0.0
		self.point = point
		self.lattice = lattice
		self.params = params

	def __del__(self):
		self.point = None
		self.lattice = None
		self.params = None

	def __repr__(self):
		return('\nDeath (n=' + repr(self.point.niche) 
					+ ', a=' + repr(self.point.age) 
					+ ', v=' + repr(self.point.vulnerability) 
			   + ') | rate=' + repr(self.rate))

	def refresh_rate(self):
		self.rate = self.params.deathrate * self.lattice.n(self.point)

	def do(self):
		self.lattice.annihilate(self.point)

class CatastropheProcess(object):
	def __init__(self, vuln, lattice, params):
		self.order = 0
		self.rate = params.catastropherate
		self.vulnerability = vuln
		self.lattice = lattice
		self.params = params

	def __del__(self):
		self.lattice = None
		self.params = None

	def __repr__(self):
		return('\nVuln Catastrophe | rate=' + repr(self.rate))

	@staticmethod
	def stencile(lattice, params):
		for vuln in xrange(params.vulnerabilityWidth):
			yield vuln

	def refresh_rate(self):
		self.rate = self.params.catastropherate

	def do(self):
		vulnerable_points = [p for p in self.lattice.sites if p.vulnerability == self.vulnerability]
		for p in vulnerable_points:
			self.lattice.extinction(p)



# class MutationProcess(object):

# 	def __init__(self, point, lattice, params):
# 		self.order = 1
# 		self.rate = 0.0
# 		self.point = point
# 		self.lattice = lattice
# 		self.params = params

# 	def __del__(self):
# 		self.point = None
# 		self.lattice = None
# 		self.params = None
	
# 	def refresh_rate(self):
# 		self.rate = self.params.mutationrate * self.lattice.n(self.point)

# 	def do(self):
# 		#new_mutant = Point(niche=(self.point.niche + sample([-self.params.nichestep, self.params.nichestep],1)[0] + (random() - 0.5*self.params.nichestep)) % self.params.nicheWidth, age=self.lattice.t())
# 		#new_mutant = Point(niche=(self.point.niche + sample([-1,1],1)[0]) % nicheWidth, age=self.lattice.t(), vulnerability=0)

# 		self.lattice.create(new_mutant)
# 		return [self.point, new_mutant]

class NicheStepMutationProcess(object):

	def __init__(self, point, lattice, params):
		self.order = 1
		self.rate = 0.0
		self.point = point
		self.lattice = lattice
		self.params = params

	def __del__(self):
		self.point = None
		self.lattice = None
		self.params = None
	
	def __repr__(self):
		return('\nNiStepMut (n=' + repr(self.point.niche) 
					  + ', a=' + repr(self.point.age) 
					  + ', v=' + repr(self.point.vulnerability) 
				+  ') | rate=' + repr(self.rate))

	def refresh_rate(self):
		self.rate = self.params.nichemutationrate * self.lattice.n(self.point)

	def do(self):
		#step_dbn = [-2*self.params.nichestep, -self.params.nichestep, -self.params.nichestep, 0, self.params.nichestep, self.params.nichestep, 2*self.params.nichestep]
		#step_dbn = [-self.params.nichestep, -self.params.nichestep, 0, self.params.nichestep, self.params.nichestep]
		#niche_step = sample(step_dbn, 1)[0]
		niche_step = sample(self.params.nicheKernel, 1)[0]
		new_niche = (self.point.niche + niche_step) % self.params.nicheWidth

		new_mutant = Point(niche=new_niche, age=self.lattice.t(), vulnerability=self.point.vulnerability)
		self.lattice.create(new_mutant, parent_age=self.point.age)
		return ('mutation', self.point, new_mutant)


class VulnerabilityStepMutationProcess(object):

	def __init__(self, point, lattice, params):
		self.order = 1
		self.rate = 0.0
		self.point = point
		self.lattice = lattice
		self.params = params

	def __del__(self):
		self.point = None
		self.lattice = None
		self.params = None
	
	def __repr__(self):
		return('\nVuStepMut (n=' + repr(self.point.niche) 
					  + ', a=' + repr(self.point.age) 
					  + ', v=' + repr(self.point.vulnerability) 
				+  ') | rate=' + repr(self.rate))

	def refresh_rate(self):
		self.rate = self.params.vulnerabilitymutationrate * self.lattice.n(self.point)

	def do(self):
		#vuln_dbn = [-self.params.vulnstep, 0, self.params.vulnstep]
		#vulnerability_step = sample(vuln_dbn, 1)[0]
		vulnerability_step = sample(self.params.vulnKernel, 1)[0]
		new_vuln = (self.point.vulnerability + vulnerability_step) % self.params.vulnerabilityWidth

		new_mutant = Point(niche=self.point.niche, age=self.lattice.t(), vulnerability=new_vuln)
		
		self.lattice.create(new_mutant, parent_age=self.point.age)
		return ('mutation', self.point, new_mutant)

class IntraCompetitionProcess(object):
	def __init__(self, point, lattice, params):
		self.order = 1 # Doesn't interact across niches or ages, so order 1. This is slightly inconsistent, should find a way to fix.
		self.rate = 0.0
		self.point = point
		self.lattice = lattice
		self.params = params

	def __del__(self):
		self.point = None
		self.lattice = None
		self.params = None

	def __repr__(self):
		return '\nIC (n=' + repr(self.point.niche) + ', a=' + repr(self.point.age) + ', v=' + repr(self.point.vulnerability) + ') | rate=' + repr(self.rate)

	def refresh_rate(self):
		n = self.lattice.n(self.point)
		self.rate = self.params.intracomprate * n * (n-1)

	def do(self):
		self.lattice.annihilate(self.point)
		return ('intracomp', self.point)

class NicheNNInterCompetitionProcess(object):
	def __init__(self, point1, point2, lattice, params):
		self.order = 2
		self.rate = 0.0
		self.point1 = point1
		self.point2 = point2
		self.lattice = lattice
		self.params = params

	def __del__(self):
		self.point1 = None
		self.point2 = None
		self.lattice = None
		self.params = None

	def __repr__(self):
		return '\nNNC (n=' + repr(self.point1.niche) + ', a=' + repr(self.point1.age) + ', v=' + repr(self.point1.vulnerability) +  ') <- (n=' + repr(self.point2.niche) + ', a=' + repr(self.point2.age) + ', v=' + repr(self.point2.vulnerability) +  ') | rate=' + repr(self.rate)

	# Starting from point, generate all NN interactions
	# consistent with creating a new individual at 'point'
	@staticmethod
	def stencile(point, lattice, params):
		# You MUST walk on everybody because all ages and other dimensions
		# aside from niche might be impicated in a NN competitive interaction
		# ... this is so ugly, there has to be a better way

		for p2 in lattice.sites:
			if ((point.niche + params.nichestep) % params.nicheWidth == p2.niche or (point.niche - params.nichestep) % params.nicheWidth == p2.niche) or (point.niche == p2.niche and point != p2):
				yield (point, p2)
				yield (p2, point)


	def refresh_rate(self):
		self.rate = self.params.intercomprate * self.lattice.n(self.point1) * self.lattice.n(self.point2)

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


class Lattice(object):
	'''The Lattice is the interface and the glue between the Metaprocess containing\n
	all active processes and the internal degrees of freedom. All creation-annihilation must go
	through the Lattice so that it controls the affect of processes that touch
	at a given site.'''
	def __init__(self, meta, zerothorder, firstorder, secondorder, params):
		self.sites = dict()
		self.attached_processes = dict()
		self.phylogeny = {'leaves':dict(), 'nodes':dict()}
		self.meta = meta
		self.classes0th = zerothorder
		self.classes1st = firstorder
		self.classes2nd = secondorder
		self.params = params
		self.total_n = 0

		# Now attach 0-order processes to the lattice itself
		self.attached_processes[0] = []
		for ZPrtemplate in self.classes0th:
			for vuln in ZPrtemplate.stencile(self, self.params):
				process = ZPrtemplate(vuln, self, self.params)
				self.attached_processes[0].append(process)
		self.meta.push(self.attached_processes[0])

	def append_phylogeny(self, point, parent_age):
		if self.phylogeny['leaves'].has_key(parent_age):
			if point.age != parent_age:
				self.phylogeny[point.age] = ('leaf', parent_age)
			else:
				# parent_age and point.age are the same,
				# this is just a birth event in the same species
				pass
		else: # parent_age not found in phylogeny, this is 
			pass 

	def create(self, point, parent_age):
		self.total_n += 1
		self.append_phylogeny(point, parent_age=parent_age)
		#if point in self.sites.keys():
		if self.sites.has_key(point):
			# Pull affected processes
			affected = self.attached_processes[point]
			self.meta.pull(affected)
			# Refresh affected processes
			self.sites[point] += 1
			for pr in affected:
				pr.refresh_rate()
			# Push back affected processes
			#self.meta.push_processes([a for a in affected if a.rate > 0.0])
			self.meta.push(affected)
		else:
			# Create new site
			self.sites[point] = 1
			self.attached_processes[point] = []
			# Generate and attach its 1st order processes...
			for Prtemplate in self.classes1st:
				process = Prtemplate(point, self, self.params)
				self.attached_processes[point].append(process)
			# ...also its 2nd order processes according to the stencile
			for Prtemplate in self.classes2nd:
				for p1, p2 in Prtemplate.stencile(point, self, self.params):
					process = Prtemplate(p1, p2, self, self.params)
					self.attached_processes[p1].append(process)
					self.attached_processes[p2].append(process)

			# Refresh the rate of these newly minted processes
			for pr in self.attached_processes[point]:
				pr.refresh_rate()

			# And push them into the metaprocess
			#self.meta.push_processes([a for a in self.attached_processes[point] if a.rate > 0.0])
			self.meta.push(self.attached_processes[point])

	def annihilate(self, point):
		self.total_n -= 1
		if self.sites.has_key(point):
			# # Annihilate and refresh their rates
			self.sites[point] -= 1
			self.refresh_and_delete_attached(point)
		else:
			# So you're trying to annihilate at some emtpy site? Something went horribly wrong
			raise ValueError('Orphan process, annihilation at empty site ' + repr(point))

	def extinction(self, point):
		self.total_n -= self.sites[point]
		self.sites[point] = 0
		self.refresh_and_delete_attached(point)

	def refresh_and_delete_attached(self, point):
		affected = self.attached_processes[point]
		self.meta.pull(affected)
		for pr in affected:
			pr.refresh_rate()
		if self.sites[point] > 0:
			# Push them back into the meta process if the site is still active
			#self.meta.push_processes([a for a in affected if a.rate > 0.0])
			self.meta.push(affected)
		else:
			# Otherwise delete the point from the dictionary...
			del self.sites[point]
			# ...along with its attached processes
			for pr in self.attached_processes[point]:
				# but just before that, make sure that these processes'
				# references lingering at point2 get cleaned first
				if pr.order == 2:
					if pr.point1 != point:
						self.attached_processes[pr.point1].remove(pr)
					if pr.point2 != point:
						self.attached_processes[pr.point2].remove(pr)
			# ...and delete these attached processes for real now
			del self.attached_processes[point]
			# They were pull, so nothing's left to push.

	# This is used by Process.refresh_rate methods
	def n(self, point):
		#if point in self.sites.keys():
		if self.sites.has_key(point):
			return self.sites[point]
		else:
			return 0

	# This is used for now by mutation/speciation processes to give an age
	# to new types. Processes are unaware of their metaprocess and should
	# always transact through the Lattice.
	def t(self):
		return meta._t




class Metaprocess(object):
	def __init__(self, params):
		
		# The gillespie algorithm just uses a bag of processes
		self._process_bag = set()

		# Actually the only reason we use a SortedSet instead of a simple set
		# is to easily obtain the maximum and minimum process rate > 0.0.
		# This is slightly overkill. The rejection step doesn't really need this
		# overhead if we could find a better way to get the max rate.
		self._sorted_processes = SortedSet(key=lambda pr: pr.rate)

		# For the composition-rejection step, we need to log2-bin all processes (composition)
		# and then execute a rejection within the bin. A bin is selected
		# in proportion to the total rate it contains. Each process contain
		# the bin it is currently in so that the pull-push protocol knows
		# where to pull all processes attached to a site together with the selected
		# process without having to search for them.
		# We use a dictionary of sets. The dictionary key will be the
		self._process_bins = dict()

		self._total_rate = 0.0
		self._num_processes = 0
		self._t = 0.0
		self.intervaltries = 0
		self.processed = 0
		self.p0 = params.pmin

		if params.step_algorithm == 'gillespie':
			self.pull = self.pull_processes
			self.push = self.push_processes
			self.step = self.step_gillespie
		elif params.step_algorithm == 'rejection':
			self.pull = self.pull_processes
			self.push = self.push_processes
			self.step = self.step_rejection
		elif params.step_algorithm == 'composition-rejection':
			self.pull = self.pull_binned_processes
			self.push = self.push_binned_processes
			self.step = self.step_compositionrejection

	@property
	def number_of_active_processes(self):
	    return max(len(self._process_bag),
	    	len(self._sorted_processes),
	    	sum([len(b.processes_in_bin) for b in self._process_bins.itervalues()]))
	
	@property
	def sum_of_process_rates(self):
		return max(sum(pr.rate for pr in self._sorted_processes), 
		sum(pr.rate for pr in self._process_bag), 
		sum(pr.rate for pr in chain.from_iterable([b.processes_in_bin for b in self._process_bins.itervalues()])))


	# Used by direct Gillespie algorithm (simple accumulator)
	def pull_processes(self, processes):
		active_processes = {pr for pr in processes if pr.rate > 0.0}
		self._process_bag -= active_processes
		self._total_rate -= sum((pr.rate for pr in active_processes))
		self._num_processes -= len(active_processes)		

	def push_processes(self, processes):
		active_processes = {pr for pr in processes if pr.rate > 0.0}
		self._total_rate += sum(pr.rate for pr in active_processes)
		self._num_processes += len(active_processes)
		self._process_bag |= active_processes

	# Used by rejection algorithm (might want to do away with the sortedcontainer)
	def pull_sorted_processes(self, processes):
		active_processes = [pr for pr in processes if pr.rate > 0.0]
		self._sorted_processes -= active_processes
		#for pr in active_processes:
		#	self._sorted_processes.remove(pr)
		self._total_rate -= sum((pr.rate for pr in active_processes))
		self._num_processes -= len(active_processes)		

	def push_sorted_processes(self, processes):
		active_processes = [pr for pr in processes if pr.rate > 0.0]
		self._total_rate += sum(pr.rate for pr in active_processes)
		self._num_processes += len(active_processes)
		self._sorted_processes |= active_processes
		#self._sorted_processes.update(active_processes)

	# Used by the composition-rejection algorithm
	def pull_binned_processes(self, processes):
		active_processes = [(floor(log(pr.rate/self.p0, 2)), pr) for pr in processes if pr.rate > 0.0]
		for li, pr in active_processes:
			# Silently do nothing if it's not in process set. This occurs if the last time
			# this process refreshed its rate went to 0.0
			self._process_bins[li].processes_in_bin -= {pr}
			self._process_bins[li].binrate -= pr.rate
			if len(self._process_bins[li].processes_in_bin) == 0:
				del self._process_bins[li].processes_in_bin
				del self._process_bins[li]
		self._total_rate -= sum(pr.rate for li, pr in active_processes)
		self._num_processes -= len(active_processes)

	def push_binned_processes(self, processes):
		active_processes = [(floor(log(pr.rate/self.p0, 2)), pr) for pr in processes if pr.rate > 0.0]
	
		self._total_rate += sum(pr.rate for li, pr in active_processes)
		self._num_processes += len(active_processes)
		for li, pr in active_processes:
			#if li in self._process_bins.keys():
			if self._process_bins.has_key(li):
				self._process_bins[li].binrate += pr.rate
				self._process_bins[li].processes_in_bin |= {pr}
			else:
				#self._process_bins[li] = ProcessBin(binrate=pr.rate, processes_in_bin = SortedSet({pr},key=lambda pr: pr.rate))
				self._process_bins[li] = ProcessBin(binrate=pr.rate, processes_in_bin = {pr})

	# Random bag of processes used by the gillespie algorithm


	def step_gillespie(self):
		if len(self._process_bag) > 0:
			r1 = random()
			self._t += -log(r1)/self._total_rate
			r2 = random()
			acc = 0.0
			for pr in self._process_bag:
				acc += pr.rate
				self.intervaltries += 1

				if r2*self._total_rate < acc:
#					self.intervaltries += 1
					residue = pr.do()
					break
			self.processed += 1

		else:
			print('Party\'s over, everybody\'s dead')
			exit()

		return residue


	def step_rejection(self):
		#residue = None
		#if self._total_rate > 0.0:
		#if True:
		# This one's the right one, the true absorbing state; no relying on floats:
		if len(meta._process_bag) > 0:
			r1 = random()
			delta_t = -log(r1)/self._total_rate
			self._t += delta_t

			# For now processes of rate 0.0 might creep in and
			# Slepoy's et al algorithm cannot work without picking
			# the first nonzero rate.
			# min_rate = self._sorted_processes[0].rate
			#max_rate = self._sorted_processes[-1].rate
			#max_rate = max(pr.rate for pr in self._sorted_processes)
			max_rate = max((pr.rate for pr in self._process_bag))

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
				
				#ri = randint(0, len(self._sorted_processes)-1)
				#process = self._sorted_processes[ri]
				process = sample(self._process_bag, 1)[0]

				if process.rate >= uniform(0, max_rate):
					residue = process.do()
					self.processed += 1
					self.intervaltries += i+1
					break

		else:
			print('Party\'s over, everybody\'s dead')
			exit()

		return residue


	def step_compositionrejection(self):
		residue = None
		if len(self._process_bins) > 0:
			r1 = random()
			delta_t = -log(r1)/self._total_rate
			self._t += delta_t

			rG = random()
			accG = 0.0
			# Get at the logarithmic bin in a direct Gillespie manner (composition)
			for logbin in self._process_bins.itervalues(): # Don't iterate otherwise the dict change
				accG += logbin.binrate
			#	self.intervaltries += 1
				if rG*self._total_rate < accG:
					selected_bin = logbin
					break
			# Enter the bin...
			max_bin_rate = max(pr.rate for pr in selected_bin.processes_in_bin)
			for i in count():
				process = sample(selected_bin.processes_in_bin,1)[0]
				if process.rate >= uniform(0, max_bin_rate):
					residue = process.do()
					self.processed += 1
					self.intervaltries += i+1
					break

		else:
			print('Party\'s over, everybody\'s dead')
			exit()

		return residue


def summary(old_time, new_time):
	return ('processed = ' + repr(meta.processed) 
		+ ' | time = ' + '%.1f'%meta._t
		+ ' | total pop = ' + repr(lattice.total_n) 
		+ ' | num. sps. = ' + repr(len(lattice.sites)) 
		+ ' | active proc. = ' + repr(meta.number_of_active_processes) 
		+ ' | tries in 1K = ' + repr(meta.intervaltries) 
		+ ' | avg. tries = ' + '%.1f'%(float(meta.intervaltries)/params.summaryinterval) 
		#+ ' | min rate = ' + '%.2f'%(meta._sorted_processes[0].rate) 
		#+ ' | max rate = ' + '%.2f'%meta._sorted_processes[-1].rate 
		+ ' | delta total = ' + '%.2e'%(meta._total_rate - meta.sum_of_process_rates) 
		+ ' | in %.2f'%(new_time-old_time) + 'secs')

params = Parameters()

if params.fixed_seed: seed(params.fixed_seed)

meta = Metaprocess(params)

lattice = Lattice(meta, 
	[CatastropheProcess], 
	[BirthProcess, DeathProcess, IntraCompetitionProcess, NicheStepMutationProcess, VulnerabilityStepMutationProcess], 
	[NicheNNInterCompetitionProcess], 
	params)
#lattice = Lattice(meta, [], [BirthProcess, DeathProcess, IntraCompetitionProcess, VulnerabilityStepMutationProcess], [NicheNNInterCompetitionProcess], params)
#lattice = Lattice(meta, [], [BirthProcess, DeathProcess, IntraCompetitionProcess], [], params)

# p1 = Point(niche=15,age=0.0,vulnerability=0)
# p4 = Point(niche=5,age=0.0,vulnerability=4)

# for i in xrange(5):
# 	lattice.create(p1)
# 	lattice.create(p4)

# for i in xrange(0,21,3):
# 		lattice.create(Point(niche=i, age=0.0, vulnerability=0))
# 		lattice.create(Point(niche=i, age=0.0, vulnerability=6))
# 		lattice.create(Point(niche=i, age=0.0, vulnerability=12))
# 		lattice.create(Point(niche=i, age=0.0, vulnerability=18))

for i in xrange(20):
	lattice.create(Point(niche=10.0, age=0.0, vulnerability=10.0), parent_age=0.0)

# print(meta.__dict__)

# exit()

if params.visualization == 'pretty':
	from pyprocessing import *
	iwanthue = dict()

	def setup():
		global old_time, new_time
		old_time = time.time()
		#frameRate(100)
		size(800,800)
		smooth()
		ellipseMode(CENTER)
		noStroke()
		fill(0)
		rect(0,0,800,800)

	def draw():
		# Fugly
		global old_time, new_time

		for i in xrange(params.drawinterval):
			residue = meta.step()
			if meta.processed % params.summaryinterval == 0:
				new_time = time.time()
				print(summary(old_time, new_time))
				old_time = new_time
				meta.intervaltries = 0


		rectMode(CORNER)
		noStroke()
		fill(0)
		rect(0,0,800,50)
		fill(200)
		textSize(15)
		text(' t = %.5f' % meta._t 
			+ ' n = ' + repr(lattice.total_n) 
			+ ' av. tries = ' + '%.1f'%(float(meta.intervaltries)/float(params.summaryinterval))
			+ ' processed = ' + repr(meta.processed)
			,25,25)

		fill(0,5)
		rect(0,0,800,800)
		rectMode(RADIUS)
		y = (meta._t/0.01) % 700
		for i,(point, n) in enumerate(lattice.sites.iteritems()):
			#if point in iwanthue.keys():
			if iwanthue.has_key(point):
				fill(*iwanthue[point])
			else:
				new_hue = sample(params.palette,1)[0]
				iwanthue[point] = new_hue
				fill(*iwanthue[point])

			#r = 5*sqrt(n)
			r = n
			#ellipse(50 + point.niche/10.0*700.0, 50 + y ,r,r)
			rect(50 + float(point.niche)/params.nicheWidth*700.0 + 0*point.age, 50 + y + point.age, r, 1)

		if residue != None:
			if residue[0] == 'mutation':
				stroke(102,210,136)
				fill(102,210,136)
				niche_from = 50 + float(residue[1].niche)/params.nicheWidth*700.0 + 0*residue[1].age
				niche_mutant = 50 + float(residue[2].niche)/params.nicheWidth*700.0 + 0*residue[2].age
				line(niche_from, 50 + y + residue[1].age, niche_mutant, 50 + y + residue[2].age)
				ellipse(niche_mutant, 50 + y + residue[2].age, 7, 7)
			elif residue[0] == 'intracomp':
				stroke(193,79,77)
				fill(193,79,77)
				ellipse(50 + float(residue[1].niche)/params.nicheWidth*700.0 + 0*residue[1].age, 50 + y + residue[1].age, 7, 7)
			elif residue[0] == 'intercomp':
				stroke(193,79,77)
				fill(193,79,77)
				inter_from = 50 + float(residue[1].niche)/params.nicheWidth*700.0 + 0*residue[1].age
				inter_to = 50 + float(residue[2].niche)/params.nicheWidth*700.0 + 0*residue[2].age
				line(inter_from, 50 + y + residue[1].age, inter_to, 50 + y + residue[2].age)
				ellipse(inter_to, 50 + y + residue[2].age, 7, 7)
		# Draw a red dot in the middle of the population
		# for i,(point, n) in enumerate(lattice.sites.iteritems()):
		# 	y = (meta._t/0.01) % 700
		# 	fill(255,0,0)
		# 	rect(50 + point.niche/10.0*700.0, 50 + y, 1, 1)

	run()
elif params.visualization == 'pretty3d':
	from pyprocessing import *
	iwanthue = dict()
	h = 0.0

	def setup():
		size(1200, 840)
		#size(1900,1000)
		#size(displayWidth, displayHeight)
		#frameRate(100)
		background(20)
		ellipseMode(CENTER)

	def draw():
		global h
		for i in xrange(params.drawinterval):
			residue = meta.step()

		background(20)

		fill(150)
		textSize(15)
		text(' t = %.1f' % meta._t 
			+ ' n = ' + repr(lattice.total_n) 
			+ ' processed = ' + repr(meta.processed)
			,25,25)


		h = max(max(p.age for p in lattice.sites.iterkeys()), h)

		#camera(105, -130, -80+h, 105, 20, -40+h, 0, 0, 1);
		camera(105 + 200*sin(frame.count*0.01), 105 - 200*cos(frame.count*0.01), -80+h, 105, 105, -40+h, 0, 0, 1);
	
		#translate(-100, -100, 0)

		stroke(50)
		strokeWeight(1)
		line(0.0, 0.0, h, params.nicheWidth*params.pptrait, 0.0, h)
		line(0.0, 0.0, h, 0.0, params.vulnerabilityWidth*params.pptrait, h)
		line(0.0, params.vulnerabilityWidth*params.pptrait, h, params.nicheWidth*params.pptrait, params.vulnerabilityWidth*params.pptrait, h)
		line(params.nicheWidth*params.pptrait, 0, h, params.nicheWidth*params.pptrait, params.vulnerabilityWidth*params.pptrait, h)


		oldest = min(p.age for p in lattice.sites.iterkeys())
		youngest = max(p.age for p in lattice.sites.iterkeys())

		line(params.nicheWidth*params.pptrait, 0, -1, params.nicheWidth*params.pptrait, 0, oldest)
		stroke(150,50,50)
		strokeWeight(2)
		line(params.nicheWidth*params.pptrait, 0, oldest, params.nicheWidth*params.pptrait, 0, youngest)
		stroke(50)
		strokeWeight(1)
		line(params.nicheWidth*params.pptrait, 0, youngest, params.nicheWidth*params.pptrait, 0, youngest+100)

		for i, (point, n) in enumerate(lattice.sites.iteritems()):

			if not iwanthue.has_key(point):
				new_hue = sample(params.palette,1)[0]
				iwanthue[point] = new_hue

			r = 1.0*n
			y = 1.0*point.niche*params.pptrait
			x = 1.0*point.vulnerability*params.pptrait
			z = 1.0*point.age

			noFill()
			stroke(*iwanthue[point])
			line(x,y,h,x,y,z)

			pushMatrix()
			translate(x,y,h)
			stroke(50)
			ellipse(0.0, 0.0, r, r)
			noStroke()
			translate(0.0, 0.0, z-h)
			fill(*iwanthue[point])
			ellipse(0.0, 0.0, r, r)
			popMatrix()
			
			# dh = 20
			# for z in [-2*dh, -dh, 0, dh, 2*dh]:
			# 	gridz = floor(h/dh)*dh + z
			# 	alphaz = 255*(tanh(2*dh + gridz - h) + tanh(2*dh - gridz + h))/2.0
			# 	stroke(200, 200, 200, alphaz)
			# 	line(0, 0, gridz, 200, 0, gridz)
			# 	line(200, 0, gridz, 200, 200, gridz)
			# 	line(200, 200, gridz, 0, 200, gridz)
			# 	line(0, 200, gridz,0, 0, gridz)
	
	run()

elif params.visualization == 'console':
	old_time = time.time()
	for i in count():

		meta.step()

		# Return a summary 
		if meta.processed % params.summaryinterval == 0:
			new_time = time.time()
			print(summary(old_time, new_time))
			old_time = new_time
			meta.intervaltries = 0
		# We've substracted and added from and to the total_rate variable and group rates of log-bined process sets a lot.
		# So to make sure that floating point underflow don't bit us back, we recompute all of these every few hundred-thousand steps.
		if meta.processed % params.recalibrateinterval == 0:
			if params.step_algorithm == 'gillespie' or params.step_algorithm == 'rejection':
				meta._total_rate = sum(pr.rate for pr in meta._sorted_processes)
			elif params.step_algorithm == 'composition-rejection':
				meta._total_rate = sum(pr.rate for pr in chain.from_iterable(logbin.processes_in_bin for logbin in meta._process_bins.itervalues()))
			
