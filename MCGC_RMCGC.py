#! /usr/bin/env python

import sys 
import os
import math
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import random
import seaborn as sns
import time


def build_Poisson_graph(N, z):
    
    adj = np.zeros((N, N))
    for i in range(N):
        for j in range(i):
            if random.uniform(0, 1) < 1.0*z/N:
                adj[i, j] = 1
                adj[j, i] = 1                
    return nx.from_numpy_matrix(adj)  

def build_SF_graph(N, k_min, gamma):

    norm = 0
    theta = np.zeros(N)
    adj = np.zeros((N, N))
    for i in range(N):
        theta[i] = k_min * random.uniform(0, 1)**(1/(1 - gamma))
        while theta[i] > math.sqrt(N):
            theta[i] = k_min * random.uniform(0, 1)**(1/(1 - gamma))
        norm += theta[i]

    for i in range(N):
        for j in range(i):
            if random.uniform(0, 1) < theta[i]*theta[j]/norm:
                adj[i, j] = 1
                adj[j, i] = 1

    return nx.from_numpy_matrix(adj)    

def ER_SF_RMCGC(Nsim, N, kmin, gamma, prob, z, graphtype):

	S_av_RMCGC = []
	for n in range(Nsim):
		print n
		if graphtype == 'SF':
			G_a = build_SF_graph(N, kmin, gamma)
			G_d = build_SF_graph(N, kmin, gamma)
			G_u = build_SF_graph(N, kmin, gamma)
		elif graphtype == 'ER':
			G_a = nx.gnp_random_graph(N, prob)
			G_d = nx.gnp_random_graph(N, prob)
			G_u = nx.gnp_random_graph(N, prob)
		elif graphtype == 'PO':
			G_a = build_Poisson_graph(N, z)
			G_d = build_Poisson_graph(N, z)
			G_u = build_Poisson_graph(N, z)			
			
		
		init = np.zeros((N, 2))
		for i in range(1, N+1):
			init[i-1, 0] = i
			init[i-1, 1] = random.uniform(0, 1)

		S = np.zeros(len(np.arange(0, 1.0, 0.01)))
		for ix, p in enumerate(np.arange(0, 1.0, 0.01)):

			#make copies of the original network
			G_a_temp = G_a.copy()
			G_d_temp = G_d.copy()
			G_u_temp = G_u.copy()

			#get the nodes to be removed
			to_be_removed = init[init[:, 1]<p, 0].astype(int)
			
			if len(to_be_removed)==0:
				G_a_LCC = max(nx.connected_component_subgraphs(G_a_temp), key=len)
				G_d_LCC = max(nx.connected_component_subgraphs(G_d_temp), key=len)
				G_u_LCC = max(nx.connected_component_subgraphs(G_u_temp), key=len)
				S[ix] = len(set(G_a_LCC.nodes()) & set(G_d_LCC.nodes()) & set(G_u_LCC.nodes()))/float(N)
			else:
				while len(to_be_removed)>0:
					#initial damage on all networks 
					G_a_temp.remove_nodes_from(to_be_removed)
					G_d_temp.remove_nodes_from(to_be_removed)
					G_u_temp.remove_nodes_from(to_be_removed)

					to_be_removed = []
					if (len(G_a_temp)==0) | (len(G_d_temp)==0) | (len(G_u_temp)==0):
						S[ix] = 0.0
						break
					else:
						#calculate LCC of each network
						G_a_LCC = max(nx.connected_component_subgraphs(G_a_temp), key=len)
						G_d_LCC = max(nx.connected_component_subgraphs(G_d_temp), key=len)
						G_u_LCC = max(nx.connected_component_subgraphs(G_u_temp), key=len)

						#compare LCCs and get the extra nodes to be deleted
						to_be_removed.extend(set(G_a_LCC) - (set(G_d_LCC) | set(G_u_LCC)))
						to_be_removed.extend(set(G_d_LCC) - (set(G_a_LCC) | set(G_u_LCC)))
						to_be_removed.extend(set(G_u_LCC) - (set(G_d_LCC) | set(G_a_LCC)))

				S[ix] = len(set(G_a_LCC.nodes()) & set(G_d_LCC.nodes()) & set(G_u_LCC.nodes()))/float(N)

		S_av_RMCGC.append(S)    

	return S_av_RMCGC

def ER_SF_MCGC(Nsim, N, kmin, gamma, prob, z, graphtype):

	S_av_MCGC = []
	for n in range(Nsim):
		print n
		if graphtype == 'SF':
			G_a = build_SF_graph(N, kmin, gamma)
			G_d = build_SF_graph(N, kmin, gamma)
			G_u = build_SF_graph(N, kmin, gamma)
		elif graphtype == 'ER':
			G_a = nx.gnp_random_graph(N, prob)
			G_d = nx.gnp_random_graph(N, prob)
			G_u = nx.gnp_random_graph(N, prob)
		elif graphtype == 'PO':
			G_a = build_Poisson_graph(N, z)
			G_d = build_Poisson_graph(N, z)
			G_u = build_Poisson_graph(N, z)	
			
		init = np.zeros((N, 2))
		for i in range(1, N+1):
			init[i-1, 0] = i
			init[i-1, 1] = random.uniform(0, 1)

		S = np.zeros(len(np.arange(0, 1.0, 0.01)))
		for ix, p in enumerate(np.arange(0, 1.0, 0.01)):

			#make copies of the original network
			G_a_temp = G_a.copy()
			G_d_temp = G_d.copy()
			G_u_temp = G_u.copy()

			#get the nodes to be removed
			to_be_removed = init[init[:, 1]<p, 0].astype(int)
			
			if len(to_be_removed)==0:
				G_a_LCC = max(nx.connected_component_subgraphs(G_a_temp), key=len)
				G_d_LCC = max(nx.connected_component_subgraphs(G_d_temp), key=len)
				G_u_LCC = max(nx.connected_component_subgraphs(G_u_temp), key=len)
				S[ix] = len(set(G_a_LCC.nodes()) & set(G_d_LCC.nodes()) & set(G_u_LCC.nodes()))/float(N)
			else:
				while len(to_be_removed)>0:
					#initial damage on all networks 
					G_a_temp.remove_nodes_from(to_be_removed)
					G_d_temp.remove_nodes_from(to_be_removed)
					G_u_temp.remove_nodes_from(to_be_removed)

					to_be_removed = []
					if (len(G_a_temp)==0) | (len(G_d_temp)==0) | (len(G_u_temp)==0):
						S[ix] = 0.0
						break
					else:
						#calculate LCC of each network
						G_a_LCC = max(nx.connected_component_subgraphs(G_a_temp), key=len)
						G_d_LCC = max(nx.connected_component_subgraphs(G_d_temp), key=len)
						G_u_LCC = max(nx.connected_component_subgraphs(G_u_temp), key=len)

						#compare LCCs and get the extra nodes to be deleted
						to_be_removed.extend(set(G_a_LCC) - (set(G_a_LCC) & set(G_d_LCC) & set(G_u_LCC)))
						to_be_removed.extend(set(G_d_LCC) - (set(G_a_LCC) & set(G_d_LCC) & set(G_u_LCC)))
						to_be_removed.extend(set(G_u_LCC) - (set(G_a_LCC) & set(G_d_LCC) & set(G_u_LCC)))

				S[ix] = len(set(G_a_LCC.nodes()) & set(G_d_LCC.nodes()) & set(G_u_LCC.nodes()))/float(N)

		S_av_MCGC.append(S)    

	return S_av_MCGC

if __name__ == "__main__":

	start = time.time()
	Nsim = int(sys.argv[1])
	N = int(sys.argv[2])
	kmin = int(sys.argv[3])
	gamma = float(sys.argv[4])
	prob = float(sys.argv[5])
	z = float(sys.argv[6])
	graphtype = sys.argv[7]
	S_av_RMCGC = ER_SF_RMCGC(Nsim, N, kmin, gamma, prob, z, graphtype)
	print 'RMCGC done'
	S_av_MCGC = ER_SF_MCGC(Nsim, N, kmin, gamma, prob, z, graphtype)
	print 'MCGC done'
	print 'Time elapsed: ', time.time() - start  
	
	sns.set(font_scale=2) 
	sns.set_style("whitegrid")
	plt.figure(figsize=(8, 6))
	plt.plot((1.0 - np.arange(0, 1.0, 0.01)), np.mean(np.array(S_av_RMCGC), 0), 'v', color='dodgerblue', label='RMCGC')	
	plt.plot((1.0 - np.arange(0, 1.0, 0.01)), np.mean(np.array(S_av_MCGC), 0), '<', color='tomato', label='MCGC')
	plt.legend(loc=0)
	plt.ylim(0, 1.0)
	plt.xlabel('p')
	plt.ylabel('S')
	plt.savefig('/Users/arda/Desktop/%s_gamma%s_%sav.pdf' % (graphtype, gamma, Nsim), format='pdf')
	plt.show()
