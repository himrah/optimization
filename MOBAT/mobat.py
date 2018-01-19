import networkx as nx
import numpy as np
from collections import Counter, defaultdict
import time
#from sklearn import metrics
#from sklearn.metrics import normalized_mutual_info_score as NMI
import operator
from operator import itemgetter
import matplotlib
#matplotlib.use("Qt4Agg")
import matplotlib.pyplot as plt
import random

class Bats:
	def __init__(self):
		self.bats_fitness_value=[]
		self.bats_fitness_values={}
		self.zReference_point = 0
		self.NS=[]
		#self.best_positions = []
		#self.best_p = []
		self.ns_size=0
		self.n_distance=[]
		self.zstr1=0
		self.zstr2=0
		self.gbest=[]
		self.gbest_mod=0
		#self.pos_modularity = 0
		self.pbest=[]
		self.fitness_value=-1
		self.velocity = []
		self.G = nx.Graph()
		self.bats = []
		self.file_name = 'dolphin.txt'
		self.number_of_bats = 50
		self.modularity = 0
		self.iteration = 10
		self.R=0.5
		self.Z1d = 0
		self.Z2d = 0
		self.A=0.5

	def Input_Graph(self):
		temp=open(self.file_name,'r').read().split('\n')
		graph=[]
		for i in temp:
		 	#pass i in temp:
			t=[]
			for j in i.split():
				if(j):
					t.append(int(j))
			if(t):
				graph.append(tuple(t))
		self.G.add_edges_from(graph)
		j=1
		for i in self.G:
			self.G.node[i]={'pos':j}
			j+=1
		#nx.set_edge_attributes(self.G,'n_dis',0)	


	def weight_init(self):
		for i in range(self.number_of_bats):
			w1=round(float(i)/float((len(self.bats)-1)),3)
			w2=round(1.0-w1,3);
			nx.set_node_attributes(self.bats[i],'weight',[w1,w2])

	def best_kkm(self):
		rs=[]
		for i in self.bats:
			rs.append(self.KKM(i))
		return max(rs)

	def best_rc(self):
		rs=[]
		for i in self.bats:
			rs.append(self.KKM(i))
		return max(rs)

	def z_reference(self):
		self.zstr1=self.best_kkm()
		self.zstr2=self.best_rc()



	def update_reference_point(self):
		if self.Z1d<self.zstr1:
			self.zstr1=self.Z1d
		if self.Z2d<self.zstr2:
			self.zstr2=self.Z2d		


	def update_pbest(self,graph):
		nob=graph.number_of_nodes()
		sum=0
		counter=0
		better=1
		pbest_sol=graph.copy()
		j=0
		for i in pbest_sol.node:
			pbest_sol.node[i]['pos']=self.pbest[j]
			j+=1

		for i in range(nob):
			if(self.KKM(graph)<=self.KKM(pbest_sol)):
				sum+=1
			elif(self.KKM(graph)<self.KKM(pbest_sol)):
				counter+=1

		if sum==nob:
			better = 0
		else:
			if sum == 0:
				better = 1
			elif counter == 1:
				t1 = ((self.KKM(graph)*(graph.node[graph.nodes()[0]]['weight'][0]))+(self.RC(graph)*(graph.node[graph.nodes()[0]]['weight'][1])))
				t2 = ((self.KKM(pbest_sol)*(pbest_sol.node[pbest_sol.nodes()[0]]['weight'][0]))+(self.RC(pbest_sol)*(pbest_sol.node[pbest_sol.nodes()[0]]['weight'][1])))
				if t1<t2:
					better = 0

		j=0
		if better == 0:
			for i in graph:
				self.pbest[j]=graph.node[i]['pos']
				j+=1
		#self.G = graph.copy()		


	def updatepos(self,graph):    #simple function to update position with base on thier neighbors frequency
		j=0
		#self.new=graph.copy()
		for i in graph:
			n=[]
			if(self.velocity[j]):
				temp=graph.neighbors(i)
				for k in temp:
					n.append(graph.node[k]['pos'])
					#print(n)
				if(len(n)!=len(set(n))):
					p=Counter(n).most_common(1)[0][0]
					graph.node[i]['pos']=p
				else:
					if(graph.node[i]['pos'] in n):
						pass
					else:
						p=np.random.choice(n)
						graph.node[i]['pos']=p
			j+=1
		self.G = graph.copy()
		return graph.copy()

	def updatevelocity(self,graph):
		fmin=0
		fmax=1
		v1=[]
		v2=[]
		v3=[]
		j=0
		p=self.pbest
		for i in graph:
			beta=round(np.random.uniform(0,1),2)
			f=fmin+(fmax-fmin)*beta
			v1.append(int((p[j]==graph.node[i]['pos']) and '0' or '1'))
			v2.append(round((self.velocity[j]+v1[j]*f),2))			
			v3.append(round(float(1/round(1+round(float(np.exp(-v2[j])),2),2)),2))
			self.velocity[j]=int((round(np.random.uniform(0,1),2)<v3[j]) and '1' or '0')
			j+=1
			

	def bats_init(self): #based on random neighborhood
		copy=self.G.copy()
		for j in range(self.number_of_bats):
			for i in copy:
				n=[]
				temp=copy.neighbors(i)
				for k in temp:
					n.append(copy.node[k]['pos'])
				if(len(n)!=len(set(n))):
					p=Counter(n).most_common(1)[0][0]
					self.G.node[i]['pos']=p
				else:
					p=np.random.choice(n)
					self.G.node[i]['pos']=p
			self.bats.append(self.G.copy())
		
	def fitness(self,graph):
		m=graph.number_of_edges()
		l=1.0/(2.0*m)
		temp=0
		for j in graph:	
			for i in graph:
				A=int(i in graph.neighbors(j))
				k1=len(list(graph.neighbors(j)))
				k2=len(list(graph.neighbors(i)))
				gama=int(graph.node[j]['pos']==graph.node[i]['pos'])
				temp+=((A-(k1*k2)/(2*m))*gama)
			
		mod=temp*l
		return np.round(mod,4)


	def KKM(self,graph):
		community=defaultdict(list)
		for i in graph:
			community[graph.node[i]['pos']].append(i)	# community{label[node,node,node...],label[node,node,node...]}
		md=0
		for com in community:
			nd=[]
			for i in community[com]:
				temp=0
				nd.append(i)
				n=graph.neighbors(i)      
				for j in community[com]:  
					if(j in n):
						pass
					else:
						temp+=1       #////////// this is L(v1-v1')
						#print(temp)
			v1=(graph.subgraph(nd).number_of_edges())*2 #  this is just like L(v1,v1) function
			v2=(graph.subgraph(nd).number_of_nodes())
			md+=float(v1)/float(v2)
		
		
		n=graph.number_of_nodes()
		k=len(community)

		res=2*(n-k)-md
		return round(res,3)    	



	def RC(self,graph):
		community=defaultdict(list)
		for i in graph:
			community[graph.node[i]['pos']].append(i)	# community{label[node,node,node...],label[node,node,node...]}
		md=0
		for com in community:
			nd=[]
			for i in community[com]:
				temp=0
				nd.append(i)
				n=graph.neighbors(i)      
				for j in community[com]:  
					if(j in n):
						pass
					else:
						temp+=1       #////////// this is L(v1-v1')
						#print(temp)
			#v1=(graph.subgraph(nd).number_of_edges())*2 #  this is just like L(v1,v1) function
			v2=(graph.subgraph(nd).number_of_nodes())
			md+=float((temp))/float(v2)
			#print(md)
		#print(round(md,2))
		return round(md,3)



	def NRA(self,graph):
    		community=defaultdict(list)
		for i in graph:
			community[graph.node[i]['pos']].append(i)	# community{label[node,node,node...],label[node,node,node...]}
		md=0
		for com in community:
			nd=[]
			for i in community[com]:
				temp=0
				nd.append(i)
				n=graph.neighbors(i)      
				for j in community[com]:  
					if(j in n):
						pass
					else:
						temp+=1       #////////// this is L(v1-v1')
						#print(temp)
			v1=(graph.subgraph(nd).number_of_edges())*2 #  this is just like L(v1,v1) function
			#v2=(graph.subgraph(nd).number_of_nodes())
			md+=float(v1)/float(v2)
			#print(md)
		#print(round(md,2))
		return round(md,3)




	def Euclidean_distance(self):
		#n_distance=[]
		i=0
		for graph in self.bats:
			temp=[]
			for g in self.bats:
				temp.append(g.copy())
			dic={}	
			del(temp[i])	
			#remain=temp.remove(graph)
			for g in temp:
				A = graph.node[graph.nodes()[0]]['weight'][0]
				B = graph.node[graph.nodes()[0]]['weight'][1]

				a = g.node[g.nodes()[0]]['weight'][0]
				b = g.node[g.nodes()[0]]['weight'][1]			

				ed = round(np.sqrt(np.power((B-A),2) + np.power((b-a),2)),1)
				dic.update({g:ed})
			sorted_dic = sorted(dic.items(), key=itemgetter(1))  #   sort dictionary
			#nx.set_node_attributes(self.bats[i],'n_distance',sorted_dic[:self.ns_size])


			nx.set_node_attributes(self.bats[i],'n_distance',[s[0] for s in sorted_dic[:self.ns_size]])
			#self.n_distance.append(sorted_dic[:self.ns_size])
			#n_distance.append(ed)
			i+=1	


	def scalar_func(fun,point,weight):
		max_fun=-1
		for i in range(1,3):
			diff = abs(fun[i]-point[i])
			if weight[i] is 0:
				feval = 0.00001 * diff
			else:
				feval = diff * weight[i]
			if feval > max_fun:
				max_fun = feval
		feval = max_fun 
		return feval			


	def update_neighborhood_solution(self,graph,child):

		#child=graph()

		ns=nx.get_node_attributes(graph,'n_distance')[1]

		for i in ns:

			point = [self.Z1d,self.Z2d]
			fun = [self.KKM(i),self.RC(i)]
			weight = nx.get_node_attributes(i,'weight')[1]
			#weight = []
			f1 = self.scalar_func(point,fun,weight)

			fun = [self.KKM(child),self.RC(RC)]
			weight = child.node[1]['weight']

			f2 = self.scalar_func(point,fun,weight)

			if f1>f2:
				graph
				nx.set_node_attributes(graph,'n_distance',child)

		self.G = graph.copy()
		#if max()

		#if max(self.KKM(new_sol)-self.zstr1,self.RC(new_sol)-self.zstr2)>max(self.KKM(old_sol)-self.zstr1,self.RC(old_sol)-self.zstr2):
			#old_sol=new_sol
			#self.modularity=self.modular_density(new_sol)


	def modular_density(self,graph):
		community=defaultdict(list)
		for i in graph:
			community[graph.node[i]['pos']].append(i)	# community{label[node,node,node...],label[node,node,node...]}
		md=0
		for com in community:
			nd=[]
			for i in community[com]:
				temp=0
				nd.append(i)
				n=graph.neighbors(i)      
				for j in community[com]:  
					if(j in n):
						pass
					else:
						temp+=1       #////////// this is L(v1-v1')
						#print(temp)
			v1=(graph.subgraph(nd).number_of_edges())*2 #  this is just like L(v1,v1) function
			v2=(graph.subgraph(nd).number_of_nodes())
			md+=(v1-temp)/v2
			#print(md)
		#print(round(md,2))
		return round(md,2)	


	def pbest_init(self):
		for i in self.bats:
			t=self.fitness(i)
			if(t>self.fitness_value):
				self.fitness_value=t
				self.pbest = [i.node[j]['pos'] for j in i]



	def gbest_init(self,bats):
		best=[]
		for i in range(bats[0].number_of_nodes()):
			best.append(0)
		f = -1
		for p in bats:
			t=self.fitness(p)
			if(t>f):
				j=0
				f=t
				for k in p:
					best[j]=p.node[k]['pos']
					j+=1
		self.gbest_mod=t
		self.gbest=best


	
	def turbulance_operation(self,graph):
		pm=.7
		temp_graph=graph.copy()
		r=np.random.choice([0,1])
		for i in graph:
			if(r<pm):
				for j in graph.neighbors(i):
					temp_graph.node[j]['pos']=graph.node[i]['pos']
		self.G=temp_graph.copy()


	def Z1dZ2d(self):
		self.Z1d=self.KKM(self.G)
		self.Z2d=self.RC(self.G)


	def draw(self):	

		pf=list((set(self.bats)))
		M = [self.KKM(i) for i in pf]
		N = [self.RC(i) for i in pf]

		F = [self.fitness(i) for i in pf]

		p1 = sorted(M)
		p1_min = min(M)
		p1_max = max(M)

		p2 = sorted(N)
		p2_min = min(N)
		p2_max = max(N)

		#print(p1)
		#print(p2)
		h=[]
		l=[]

		

		for i in range(len(M)):
			h.append(round(float((p1[i]-p1_min))/float((p1_max - p1_min)),5))

			l.append(round(float((p2[i]-p2_min))/float((p2_max - p2_min)),5))

		fig, ax = plt.subplots()

		plt.figure(1)

		ax.scatter(l, h)
		for i, txt in enumerate(F):
		    ax.annotate(txt, (l[i],h[i]))
		#plt.show()    

		#plt.plot(h, l,'ro')
		#print(h)
		#print(l)
		plt.figure(2)
		nx.draw(self.G,node_color=[self.G.node[i]['pos'] for i in self.G])

		plt.show()


	def optimize(self):

		self.__init__()
		startTime = time.time()
		self.Input_Graph()
		self.bats_init()
		t=self.G.number_of_nodes()
		vel=[]
		for ll in range(t):
			vel.append(0)

		self.velocity=vel
		self.weight_init()
		self.pbest_init()
		self.z_reference()
		self.Euclidean_distance()
		for i in range(self.iteration):
			#print('af',self.velocity)
			#print("Iteration : %d"%(i+1),end='\r')
			print("iteration : ",(i+1))
			inc=0
			for p in self.bats:
				#self.velocity=vel
				#print(p.node)
				#print('bf',self.best_positions)
				

				self.updatevelocity(p)
				#print('af',self.best_positions)
				
				t1=self.updatepos(p)
				
				#y=self.rearrange(t1)
				self.turbulance_operation(t1)
				self.Z1dZ2d()

				self.update_neighborhood_solution(t1,p)
				self.update_reference_point()
				self.update_pbest(t1)

				#print(p.node)
			print(len(set([self.G.node[i]['pos'] for i in self.G])))
			#print(len(set(self.pbest)))
		
		print("\n\n**********************************************************")	
		#print()
		
		print('\nThe script take {0} second '.format(np.round((time.time() - startTime),2)))
		#print("\nModularity is : ",fmax)
		print("\nFitness : ",self.fitness(self.G))
		print("\nNumber of Communites : ",len(set([self.G.node[i]['pos'] for i in self.G])))
		print("\nGlobal Best Position : ",[self.G.node[i]['pos'] for i in self.G])
		#nx.draw(self.G,node_color=[self.G.node[i]['pos'] for i in self.G])
		#nx.draw(self.G, node_color=pos)
		#plt.show()

		self.draw()
if __name__=='__main__':
    
	f=Bats()
	f.optimize()
