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
	camera(10.0, 150, -150.0, 0.0, 0.0, 0.0, 
	       0.0, -1.0, 0.0);
	translate(frame.count,0*frame.count,0*frame.count)
	pushMatrix()
	for xi in xrange(10):
		stroke(*c2[0])
		strokeWeight(xi)
		translate(10,0,0)
		line(-55,55,-55,-55,55,55)
	popMatrix()
	pushMatrix()
	for yi in xrange(10):
		stroke(*c2[1])
		strokeWeight(yi)
		translate(0,0,10)
		line(-55,55,-55,55,55,-55)
	popMatrix()

run()