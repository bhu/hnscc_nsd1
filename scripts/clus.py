#!/usr/bin/env python

import sys
import hdbscan
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

samp, minpts, minsamps = sys.argv[1:] 
d = np.genfromtxt('misc/%s.csv' % samp, delimiter=',')

clusterer = hdbscan.HDBSCAN(min_cluster_size = int(minpts), min_samples = int(minsamps), core_dist_n_jobs = 1).fit(d)
color_palette = sns.color_palette('deep', 8)
cluster_colors = [color_palette[x] if x >= 0
                  else (0.5, 0.5, 0.5)
                  for x in clusterer.labels_]
cluster_member_colors = [sns.desaturate(x, p) for x, p in
                         zip(cluster_colors, clusterer.probabilities_)]

np.savetxt('misc/clus.%s.%s.%s.txt' % (samp, minpts, minsamps), clusterer.labels_.astype(int), fmt = '%i')
plt.scatter(*d.T, s=50, linewidth=0, c=cluster_member_colors, alpha=0.25)
plt.xlabel('KO')
plt.ylabel('PA')
plt.savefig('misc/clus.%s.%s.%s.png' % (samp, minpts, minsamps), transparent=True)

