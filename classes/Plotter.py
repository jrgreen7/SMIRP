import random
import numpy
from collections import namedtuple
import matplotlib.pyplot as plt
import scipy
import scipy.stats

class Plotter:

	def plot_roc(self, data, title, outName):
		plt.scatter([d[1] for d in data], [d[0] for d in data])
		plt.xlabel('1-Specificity')
		plt.ylabel('Sensitivity')
		plt.axis([0, 1, 0, 1])
		plt.title(title)
		# plt.show()
		plt.savefig(outName)
		plt.clf()

	def plot_pr(self, data, title, imbalance, outName):
		plt.scatter([d[1] for d in data], [d[0] for d in data])
		plt.xlabel('Recall')
		plt.ylabel('Precision (1:'+str(imbalance)+')')
		plt.axis([0, 1, 0, 1])
		plt.title(title)
		# plt.show()
		plt.savefig(outName)

	def plot_scatter(self, data, outName):
		plt.scatter([d[0] for d in data], [d[1] for d in data])
		plt.savefig(outName)

	def plot_heatmap(self, data, outName):
		return