from random import random, seed, sample, uniform
from itertools import count, izip, permutations, combinations
from collections import namedtuple
from math import sqrt, exp, pi, log
from sortedcontainers import SortedList

def gaussian(x, mu, sigma):
	return exp(-0.5*((x-mu)/sigma)**2) / (sqrt(2*pi)*sigma)

	
class BirthProcess(object):
	@staticmethod
	def per_capita_birth_rate(incipience, niche):
		return 1.0

	def __init__(self, population1, metaprocess1):
		self._population = population1
		self._population.attach_process(self)
		self._metaprocess = metaprocess1
		self._rate = 0.0
		self._oldrate = 0.0
		self.refresh_rate()
		#self.rate = self._population.n * BirthProcess.per_capita_birth_rate(self._population.incipience, self._population.niche)

	def do(self):
		print('Birth!')
		self._population.n += 1

	def rate_and_self(self):
		return [self.rate, self]

	def refresh_rate(self):
		self._oldrate = self.rate
		self.rate = self._population.n * BirthProcess.per_capita_birth_rate(self._population.incipience, self._population.niche)
		if self._population.n <= 0:
			self._metaprocess.remove_process(self)
			self._population = None
			self._metaprocess = None

	@property
	def rate(self):
		return self._rate
	@rate.setter
	def rate(self, value):
		self._rate = value
		print('..birth rate refreshed!')
		self._metaprocess.update_meta_rates(self._oldrate, self._rate)
	
class DeathProcess(object):
	@staticmethod
	def per_capita_death_rate(incipience, niche):
		return 0.9

	def __init__(self, population1, metaprocess1):
		self._population = population1
		self._population.attach_process(self)
		self._metaprocess = metaprocess1
		self._rate = 0.0
		self._oldrate = 0.0
		self.refresh_rate()
		#self.rate = self._population.n * DeathProcess.per_capita_death_rate(self._population.incipience, self._population.niche)

	def do(self):
		print('Death!')
		self._population.n -= 1

	def rate_and_self(self):
		return [self.rate, self]

	def refresh_rate(self):
		self._oldrate = self.rate
		self.rate = self._population.n * DeathProcess.per_capita_death_rate(self._population.incipience, self._population.niche)
		# Rates are refreshed before check if the process should be removed
		# so that the meta community total rate stay consitent.
		# Oh boy, this is getting creepy if it were not already
		# This looks like SQL triggers! NOOooooo...
		if self._population.n <= 0:
			self._metaprocess.remove_process(self)
			self._population = None
			self._metaprocess = None

	@property
	def rate(self):
		return self._rate
	@rate.setter
	def rate(self, value):
		self._rate = value
		print('..death rate refreshed!')
		self._metaprocess.update_meta_rates(self._oldrate, self._rate)

#class DiversificationProcess(object):


# class CompetitionProcess(object):
# 	@staticmethod
# 	def competition_kernel():

# 	def __init__(self, population1, population2):



class Population(object):
	def __init__(self, n0=1, incipience0=0.0, niche0=0.0):
		self._n = n0
		self._incipience = incipience0
		self._niche = niche0
		self._attached_process = []

	def attach_process(self, process):
		self._attached_process.append(process)

	@property
	def n(self):
	    return self._n
	@n.setter
	def n(self, value):
		self._n = value
		for proc in self._attached_process:
			proc.refresh_rate()
		if self._n <= 0:
			self._attached_process = None


	@property
	def incipience(self):
	    return self._incipience
	@incipience.setter
	def incipience(self, value):
	    self._incipience = value
	
	@property
	def niche(self):
	    return self._niche
	@niche.setter
	def niche(self, value):
	    self._niche = value
	
class Metaprocess(object):
	def __init__(self):
		self._total_rate = 0.0
		self._max_rate = 0.0
		self._min_rate = 0.0
		self._t = 0.0
 		#self._populations = []
		self._processes = []

	def add_population(self, n0, niche0):
		new_population = Population(n0, self._t, niche0)
		#self._populations.append(new_population)
		self._processes.append(BirthProcess(new_population,self))
		self._processes.append(DeathProcess(new_population,self))

	def update_meta_rates(self, old_rate, new_rate):
		#self._total_rate += new_rate - old_rate
		#self._max_rate = max(new_rate, self._max_rate)

		# Ok, the whole garbage collection screwed my day, we'll do it
		# old-school slow to begin with, back to O(N)
		# a sortedcontainer might be reasonable, it would
		# give us the max and min for free, 
		self._total_rate = sum([p.rate for p in self._processes])
		self._max_rate = max([p.rate for p in self._processes])
		self._min_rate = min([p.rate for p in self._processes])
		print('...Meta rates updated')

	def remove_process(self, process):
		self._processes.remove(process)

	def pick_time_step(self):
		return -log(random())/self._total_rate

	def pick_process(self):
		for i in count():
			p = sample(self._processes,1)[0]
			if p.rate >= uniform(0,self._max_rate):
				print('\t' + repr(i) + ' rejections')
				break
		return p

	def do_next_process(self):
		self._t += self.pick_time_step()
		self.pick_process().do()

		

seed(652134927)

metacomm = Metaprocess()
metacomm.add_population(1, 0.0)

print(metacomm.__dict__)
for i in xrange(6):
	metacomm.do_next_process()
	print(metacomm.__dict__)


# pop1 = Population(45,3.0,0.1)
# process1 = BirthProcess(pop1)
# print('pop1 looks like\t\t' + repr(pop1) + ' = ' +repr(pop1.__dict__))
# print('process1 looks like\t\t' + repr(process1) + ' = ' + repr(process1.__dict__))
# process1.process()
# process1.process()
# print(pop1.__dict__)
# print(process1.__dict__)
# print(process1.rate_and_self())