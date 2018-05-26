#!/usr/bin/env python2

from __future__ import print_function

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

import os, subprocess, re, sys
import argparse
import numpy as np

from Bio import SeqIO

def hmmtop(sequence):
	p = subprocess.Popen(['hmmtop', '-if=--', '-pi=spred', '-sf=FAS', '-is=pseudo'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	fasta = '>{0.name}\n{0.seq}'.format(sequence)
	out, err = p.communicate(input=str(fasta))

	indices = re.findall('(?: IN| OUT)((?:\s*(?:[0-9]+))+)', out.strip())[0].strip().split()
	indices = [int(i) for i in indices[1:]]
	if not indices: return []

	seq = ''
	for c in sequence: seq += c
	seq = re.sub('[A-Z\-]', '', seq)
	for i, c in enumerate(seq):
		if c in 'ACDEFGHIKLMNPQRSTVWY': pass
		else:
			for j, n in enumerate(indices):
				if (n - 1) >= i: indices[j] += 1
	spans = [[indices[i], indices[i+1]] for i in range(0, len(indices), 2)]
	return spans


class Plot(object):
	def __init__(self):
		self.fig = Figure()

		self.canvas = FigureCanvas(self.fig)

		self.ax = self.fig.add_subplot(111)
		self.width, self.height = None, None

		self.xlim = [None, None]
		self.ylim = [None, None]

		self.xticks, self.yticks = None, None
		self.axeslabels = ['seq1', 'seq2']

		self.title = ''

	def render(self):
		#title = 'THREAT plot for seq1 and seq2'

		if self.title: self.ax.set_title(self.title)
		else: self.ax.set_title('unnamed plot')

		self.fig.gca().invert_yaxis()
		

	def save(self, filename, dpi=80, format='png'): self.fig.savefig(filename, dpi=dpi, format=format)

class Entity(object):
	def __init__(self): pass

	def draw(self, plot): raise NotImplementedError

def average_hydropathies(sequence, window=19):
	kyledolittle = {
		'A': 1.8,
		'C': 2.5,
		'D':-3.5,
		'E':-3.5,
		'F': 2.8,
		'G':-0.4,
		'H':-3.2,
		'I': 4.5,
		'K':-3.9,
		'L': 3.8,
		'M': 1.9,
		'N':-3.5,
		'P':-1.6,
		'Q':-3.5,
		'R':-4.5,
		'S':-0.8,
		'T':-0.7,
		'V': 4.2,
		'W':-0.9,
		'Y':-1.3,
	}
	rawhydro = []
	for i, aa in enumerate(sequence):
		try: rawhydro.append(kyledolittle[aa])
		except KeyError: rawhydro.append(np.nan)

	hydro = []
	for i, aa in enumerate(sequence):
		if not (window//2 < i < (len(sequence) - window//2)): continue
		hydro.append(np.nanmean(rawhydro[i-window//2:i-window//2+window]))
	return np.array(hydro)

class CorrelGraph(object):
	def __init__(self, seq1, seq2, style='inferno', window=19):
		self.seq1 = seq1
		self.seq2 = seq2
		self.style = style
		self.window = window

		self.power = 4

		dx, dy = 1, 1
		self.y, self.x = np.mgrid[slice(self.window//2, len(self.seq2)+self.window//2-window, dy), \
			slice(self.window//2, len(self.seq1)+self.window//2-window, dx)]
		self.z = None
		self.compute_hydropathies()

	def compute_hydropathies(self):
		self.hydro1 = average_hydropathies(self.seq1, window=self.window)
		self.hydro2 = average_hydropathies(self.seq2, window=self.window)
		self.correls = np.zeros([len(self.hydro1), len(self.hydro2)])

		def truncate(x):
			if x < 0: return 0.0
			elif x > 1: return 1.0
			else: return x**self.power

		for i1, h1 in enumerate(self.hydro1):
			if not (self.window//2 < i1 < (len(self.hydro1)-self.window//2)): continue
			v1 = self.hydro1[i1-self.window//2:i1-self.window//2+self.window]

			for i2, h2 in enumerate(self.hydro2):
				if not (self.window//2 < i2 < (len(self.hydro2)-self.window//2)): continue
				v2 = self.hydro2[i2-self.window//2:i2-self.window//2+self.window]
				#print(
				#	len(self.hydro1[i1-window//2:i1-window//2+window]), 
				#	len(self.hydro2[i2-window//2:i2-window//2+window])
				#)
				self.correls[i1, i2] = truncate(np.corrcoef(v1, v2)[0,1])
				#print('1:', v1)
				#print('2:', v2)

	def multiply_in_tmss(self, k=0.0):
		for span in hmmtop(self.seq1): self.correls[span[0]:span[1],:] *= k
		for span in hmmtop(self.seq2): self.correls[:,span[0]:span[1]] *= k

	def multiply_not_tmss(self, k=0.0):
		def invert(spans, length):
			newspans = []
			for i, span in enumerate(spans):

				if i == 0: newspans.append([1, span[0]-1])
				elif i != (len(spans)-1): newspans.append([spans[i][1]+1, spans[i+1][0]-1])
				else: 
					newspans.append([span[1]+1, length-1])
					return newspans

			return newspans
		for span in invert(hmmtop(self.seq1), len(self.seq1)): self.correls[span[0]:span[1],:] *= k
		for span in invert(hmmtop(self.seq2), len(self.seq2)): self.correls[:,span[0]:span[1]] *= k

	def draw(self, plot): 
		#print(self.x.shape)
		#print(self.y.shape)
		#print(self.correls.shape)
		plot.ax.pcolormesh(self.x, self.y, self.correls.T, cmap=self.style)

def main(fn1, fn2, **kwargs):
	plt = Plot()

	if 'title'in kwargs: plt.title = kwargs['title']

	if 'outfile'in kwargs and kwargs['outfile'] is not None: outfile = kwargs['outfile']
	else: outfile = 'test.png'

	if 'dpi' in kwargs and kwargs['dpi']: dpi = kwargs['dpi']
	else: dpi = 80

	seq1 = SeqIO.read(fn1, 'fasta')
	seq2 = SeqIO.read(fn2, 'fasta')
	x = CorrelGraph(seq1, seq2)
	if 'in_tms'in kwargs and kwargs['in_tms'] != 1: x.multiply_in_tmss(k=kwargs['in_tms'])
	if 'not_tms'in kwargs and kwargs['not_tms'] != 1: x.multiply_not_tmss(k=kwargs['not_tms'])
	x.draw(plt)

	plt.render()
	plt.save(outfile, dpi=dpi)

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('seq1')
	parser.add_argument('seq2')
	parser.add_argument('-o', '--outfile')
	parser.add_argument('--title')
	parser.add_argument('-r', '--dpi', type=int, default=80)
	parser.add_argument('--in-tms', type=float, default=1.0)
	parser.add_argument('--not-tms', type=float, default=1.0)

	args = parser.parse_args()

	main(args.seq1, args.seq2, title=args.title, outfile=args.outfile, dpi=args.dpi, in_tms=args.in_tms, not_tms=args.not_tms)
