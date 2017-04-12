import networkx as nx
import numpy as np
from collections import Counter, defaultdict
import time
import operator
import matplotlib.pyplot as plt

class Cat:
	def __init__(self):
		self.pbest = []
		self.SMP = 5
		self.SRD = .2
		#self.gbest = []
		#self.gbest_mod = 0
		self.velocity = []
		self.G = nx.Graph()
		self.cats = []
		self.file_name = 'l.txt'
		self.number_of_cats = 10
		self.modularity = 0
		self.iteration = 20

	def Input_Graph(self):
		temp=open(self.file_name,'r').read().split('\n')
		graph=[]
		for i in temp:
			t=[]
			for j in i.split():
				if(j):
					t.append(int(j))
			if(t):
				graph.append(tuple(t))
		self.G.add_edges_from(graph)
		j=0
		for i in self.G:
			self.G.node[i]={'pos':j}
			j+=1


	def seek(self,graph):
			


	def trace(self,graph):	






	def updatepos_simple(self,graph): #update position based on most common neighbors

		for i in graph:
			graph.node[i]['pos']+=self.velocity[i]
		return graph.copy()	


	def iweight(self,k):
		wmax=0.9
		wmin=0.4
		kmax=self.iteration
		w=((wmax-wmin)*((kmax-k)/kmax))+wmin
		return w


	def updatevelocity(self,graph):


		e=2.7182

		c1=np.random.choice([1,2])
		r=np.random.uniform(0,1)
		#for i in graph:
		#	self.velocity[i]+=(r1)*(c1)*(self.gbest[i]-graph.node[i]['pos'])

		d1=[]
		d0=[]
		#j=0
		for i in self.gbest:
			if(i==1):
				#d1[j]=r1*c1
				d1.append(r1*c1)
				d0.append(-r1*c1)
				j+=1
			else
				d1.append(-r1*c1)
				d0.append(r1*c1)
				j+=1
		v1=[0,0,0,0,0]
		v0=[0,0,0,0,0]
		j=0
		for i in range(len(v1)):
			v1[j]=w*v1[j]+d1[j]
			v0[j]=w*v0[j]+d0[j]
			j+=1

		j=0	
		for i in graph:
			if(graph.node[i][pos]==1):
				self.velocity[j]=v0[j]
				j+=1
			else:
				self.velocity[j]=v1[j]
				j+=1	

		t=[0,0,0,0]
		j=0
		for i in self.velocity:
			t[j]=1/(1+np.power(e,-(self.velocity[i])))
			j+=1
		j=0	
		for i in graph:
			if(np.random.choice([0,1])>t[j])
				#j+=0
				graph[i]['pos']=self.best[j]
				j+=1


			

		"""c1=c2=1.494
		w=self.iweight(k)
		v1=[]
		v2=[]
		v3=[]
		v4=[]
		v5=[]
		j=0
		for i in graph:
			v1.append(int((self.pbest[j]==graph.node[i]['pos']) and '0' or '1'))
			v2.append(int((self.gbest[j]==graph.node[i]['pos']) and '0' or '1'))
			r1=float(np.round(np.random.uniform(0.1,0.9),3))
			r2=float(np.round(np.random.uniform(0.1,0.9),3))
			R1=c1*r1
			R2=c2*r2
			v3.append(v1[j]*R1)
			v4.append(v2[j]*R2)
			v5.append(v3[j]+v4[j]+(self.velocity[j]*w))
			self.velocity[j]=(int((v5[j]>=1) and '1' or '0'))
			j+=1"""


	# initialization base on random neighbors		
	def cats_init(self):
		a=self.G.nodes()
		l=np.random.randint(0,2,len(a)).tolist()
		self.pbest=l

		copy=self.G.copy()
		#print(copy)
		for j in range(self.number_of_particles):
			k=0
			temp=np.random.randint(0,2,len(a)).tolist()
			for i in copy:
				self.G.node[i]['pos']=temp[k]
				k+=1
			self.particle.append(self.G.copy())
"""			for i in copy:
				n=[]
				temp=copy.neighbors(i)
				for k in temp:
					n.append(copy.node[k]['pos'])
				if(len(n)!=len(set(n))):
						p=Counter(n).most_common(1)[0][0]
						self.G.node[i]['pos']=p
				else:
						p=np.random.choice(n)
						self.G.node[i]['pos']=p"""



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



	def rearrange(self,graph):
		pos=[]
		node=graph.nodes()
		for i in graph:
			pos.append(graph.node[i]['pos'])
		new_pos=[]
		single=list(set(pos))
		for i in single:
			if(i in node):
				node.remove(i)
				#print(node)
				f=True
			else:
				f=False	
			num=np.random.choice(node)
			new_pos.append(num)
			node.remove(num)
			if(f is True):
				node.append(i)
		for i in graph:
			t=graph.node[i]['pos']
			d=single.index(t)
			graph.node[i]['pos']=new_pos[d]
			#print(graph.node[i]['pos'])
		return graph.copy()			



	def gbest_init(self,particle):
		best=[]
		for i in range(particle[0].number_of_nodes()):
			best.append(0)
		f=-1
		for p in particle:
			t=self.modular_density(p)
			if(t>f):
				j=0
				f=t
				for k in p:
					best[j]=p.node[k]['pos']
					j+=1
		self.gbest_mod=t
		self.gbest=best


	def optimize(self):
		startTime = time.time()
		self.Input_Graph()
		self.cats_init()
		#print(self.particle)
		#self.gbest_init(self.particle)
		t=self.G.number_of_nodes()
		vel=[]
		for ll in range(t):
			vel.append(0)
		self.velocity=vel
		
		for i in range(self.iteration):
			print("Iteration : %d"%(i+1),end='\r')
			#print(self.particle)
			for p in self.particle:
				#print(p.node)
				self.updatevelocity(p,i+1)
				#t1=self.updatepos(p)
				print(p.node)
				if(np.mod(i,2)):
					t1=self.updatepos_ion(p)
				else:
					t1=self.updatepos_simple(p)	
				#print(t1.node)
				t2=self.rearrange(t1)
				
				for r in p:
					p.node[r]['pos']=t2.node[r]['pos']
				m=self.modular_density(t2)
				#print("best modularity : %f"%m,end='\r')
				if(m>self.modularity):
					self.modularity=m
					k=0
					for d in t2:
						self.pbest[k]=t2.node[d]['pos']
						k+=1

			if(self.modularity>self.gbest_mod):
				self.gbest_mod=self.modularity
				
				self.gbest=self.pbest
		
		j=0		
		for i in self.G:
			self.G.node[i]['pos']=self.gbest[j]
			j+=1
		print("\n\n**********************************************************")	
		print('\nThe script take {0} second '.format(np.round((time.time() - startTime),2)))
		print("\nModularity is : ",self.gbest_mod)
		print("\nNumber of Communites : ",len(set(self.gbest)))
		print("\nGlobal Best Position : ",self.gbest)
		print("\nGraph : ",self.G.node)
		nx.draw(self.G,node_color=[self.G.node[i]['pos'] for i in self.G])
		plt.show()

if __name__=='__main__':
	f=Particle()
	f.optimize()

