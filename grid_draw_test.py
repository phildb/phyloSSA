from pyprocessing import *
from math import cos, sin
from random import sample, seed, randint

global nicheWidth, vulnerabilityWidth
nicheWidth = 21
vulnerabilityWidth = 21

def setup():
	size(1100, 700);
	noFill();

c2 = [[0,0,255],[181,30,203]]


def draw():
	background(0);
	#perspective()
	#camera(150.0*cos(0.01*frame.count), -150, 150.0*sin(0.01*frame.count), 0.0, 0.0, 0.0, 
	#       0.0, 1.0, 0.0);
	camera(50, -50+frame.count*0.2, 50, 0.0, 0.0, 0.0, 
	       0.0, -1.0, 0.0);
	seed(1)
	for yi in xrange(10):
		strokeWeight(1)
		#translate(0,-10,10)
		translate(0,-10)
		pushMatrix()
		for xi in xrange(10):
			stroke(*sample(c2,1)[0])
			strokeWeight(randint(1,4))
			translate(10,0,0)
			line(-55,55,-55,-55,55,55)
		popMatrix()
		pushMatrix()
		strokeWeight(1)
		for zi in xrange(10):
			stroke(*sample(c2,1)[0])
			strokeWeight(randint(1,4))
			translate(0,0,10)
			line(-55,55,-55,55,55,-55)
		popMatrix()


run()