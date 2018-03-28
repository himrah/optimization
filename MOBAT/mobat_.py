import networkx as nx
import numpy as np
from collections import Counter, defaultdict
import time
from operator import itemgetter
import random
import matplotlib.pyplot as plt

class Bats:	
	def __init__(self):
		self.ns_size=10
		self.zstr1=0
		self.zstr2=0
		self.gbest_mod=0
		self.pm = 0.15  #probability mutation
		self.fitness_value=-1
		self.velocity = []
		self.G = nx.Graph()
		self.bats = []
		self.file_name = 'kara.txt'
		self.number_of_bats = 20
		self.modularity = 0
		self.iteration = 10
		self.Z1d = 0
		self.Z2d = 0
		self.R=0.5
		self.A=0.5


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
		j=1
		for i in self.G:
			self.G.node[i]={'pos':j}
			j+=1


	def weight_init(self):
		for i in range(self.number_of_bats):
			w1=round(float(i)/float((len(self.bats)-1)),3)
			w2=round(1.0-w1,3);
			nx.set_node_attributes(self.bats[i],'weight',[w1,w2])

	def best_kkm(self):
		rs=[]
		for i in self.bats:
			rs.append(self.KKM(i))
		return min(rs)

	def best_rc(self):
		rs=[]
		for i in self.bats:
			rs.append(self.RC(i))
		return min(rs)

	def z_reference(self):# ideal point 
		self.zstr1=self.best_kkm()
		self.zstr2=self.best_rc()


	def update_reference_point(self):
		if self.Z1d<self.zstr1:
			self.zstr1=self.Z1d
		if self.Z2d<self.zstr2:
			self.zstr2=self.Z2d		


	def update_pbest(self,graph,child):
		nob=2
		sum=0
		counter=0
		better=1
		pbest_sol=graph.copy()
		for i in range(nob):
			if(self.KKM(child)<=self.KKM(graph)):
				sum+=1
			elif(self.KKM(child)<self.KKM(graph)):
				counter+=1
		if sum==nob:
			better = 0
		else:
			if sum == 0:
				better = 1
			elif counter == 1:
				t1 = ((self.KKM(child)*(child.node[child.nodes()[0]]['weight'][0]))+(self.RC(child)*(child.node[child.nodes()[0]]['weight'][1])))
				t2 = ((self.KKM(graph)*(graph.node[graph.nodes()[0]]['weight'][0]))+(self.RC(graph)*(graph.node[graph.nodes()[0]]['weight'][1])))
				if t1<t2:
					better = 0
		j=0
		if better == 0:
			nx.set_node_attributes(graph,'pbest',nx.get_node_attributes(child,'pbest').values()[1])


	def updatepos(self,g):    #simple function to update position with base on thier neighbors frequency
		j=0
		graph=g.copy()
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
		#self.G = graph
		return graph


	def updatevelocity(self,graph,gbest):
		fmin=0
		fmax=1
		v1=[]
		v2=[]
		v3=[]
		j=0
		for i in graph:
			beta=round(np.random.uniform(0,1),2)
			f=fmin+(fmax-fmin)*beta
			v1.append(int((gbest.node[i]['pos']==graph.node[i]['pos']) and '0' or '1'))
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
			v2=(graph.subgraph(nd).number_of_nodes())
			md+=float((temp))/float(v2)
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
			v1=(graph.subgraph(nd).number_of_edges())*2 #  this is just like L(v1,v1) function
			v2=(graph.subgraph(nd).number_of_nodes())
			md+=float(v1)/float(v2)
		return round(md,3)




	def Euclidean_distance(self):
		i=0
		t=[]
		for graph in self.bats:
			dic={}	
			f=self.bats[i]
			for g in self.bats:
				if g != f:
					A = graph.node[graph.nodes()[0]]['weight'][0]
					B = graph.node[graph.nodes()[0]]['weight'][1]
					a = g.node[g.nodes()[0]]['weight'][0]
					b = g.node[g.nodes()[0]]['weight'][1]			
					ed = round(np.sqrt(np.power((a-A),2) + np.power((b-B),2)),1)
					dic.update({g.copy():ed})
			sorted_dic = sorted(dic.items(), key=itemgetter(1))
			neigbhorhood_set=[s[0] for s in sorted_dic[:self.ns_size]]
			nx.set_node_attributes(self.bats[i],'n_distance',neigbhorhood_set)
			i+=1


	def scalar_func(self,point,fun,weight):
		max_fun=-1
		for i in range(2):
			diff = abs(fun[i]-point[i])
			if weight[i] is 0:
				feval = 0.00001 * diff
			else:
				feval = diff * weight[i]
			if feval > max_fun:
				max_fun = feval
		feval = max_fun 
		return feval			


	def UpdateInPop(self,neighbor,child):
		for i in range(len(self.bats)):
			if self.bats[i].node==neighbor.node:
				child_pos=nx.get_node_attributes(child,'pos')
				nx.set_node_attributes(self.bats[i],'pos',child_pos) 	
			i+=1


	def update_neighborhood_solution(self,ns,child):
		for i in ns:
			point = [self.zstr1,self.zstr2]
			fun = [self.KKM(i),self.RC(i)]
			weight = nx.get_node_attributes(i,'weight')[1]
			f1 = self.scalar_func(point,fun,weight)
			fun = [self.KKM(child),self.RC(child)]
			f2 = self.scalar_func(point,fun,weight)
			if f1>f2:
				self.UpdateInPop(i,child)


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
		return round(md,2)	



	def pbest_init(self):
		for i in self.bats:
			nx.set_node_attributes(i,'pbest',nx.get_node_attributes(i,'pos').values())
	
	def turbulance_operation(self,graph):
		temp_graph=graph.copy()
		#r=round(np.random.random(),2)
		r=round(np.random.uniform(0,1),2)
		for i in graph:
			if(r<self.pm):
				for j in graph.neighbors(i):
					temp_graph.node[j]['pos']=graph.node[i]['pos']
		self.G=temp_graph.copy()


	def Z1dZ2d(self,child):
		self.Z1d=self.KKM(child)
		self.Z2d=self.RC(child)


	def draw(self):	
		pf=list((set(self.bats)))		
		M = [self.KKM(i) for i in pf]
		N = [self.RC(i) for i in pf]
		F = [self.fitness(i) for i in pf]
		p1=M
		p1_min = min(M)
		p1_max = max(M)
		p2=N
		p2_min = min(N)
		p2_max = max(N)
		h=[]
		l=[]
		for i in range(len(M)):
			h.append(round(float((p1[i]-p1_min))/float((p1_max - p1_min)),5))
			l.append(round(float((p2[i]-p2_min))/float((p2_max - p2_min)),5))
		fig, ax = plt.subplots()
		plt.figure(1)
		ax.scatter(h, l)
		for i, txt in enumerate(F):
		    ax.annotate(txt, (h[i],l[i]))
		plt.figure(2)
		nx.draw(self.G,node_color=[self.G.node[i]['pos'] for i in self.G])
		plt.show()

	def eq_4(self,g,gbest):
		#new_graph=graph
		av=np.mean(self.A)
		graph=g.copy()
		av1=round(av,2)
		for i in graph:
			e=round(np.random.uniform(0,1),2)
			if(e>av1):
                
				#new_graph.node[i]['pos']=self.best_positions[k]+np.mean(self.A)*e
				
				neighbor=graph.neighbors(i)
				temp=np.array(list(neighbor))
				#new_graph[i]['pos']=self.best_positions[k]+np.random.choice(neighbor)
				graph.node[i]['pos']=np.random.choice(temp)
			else:
				graph.node[i]['pos']=gbest.node[i]['pos']
#			if(c>np.mean(self.A))
		return graph

	def increase_r(self,itr):
		#if(itr==0):
		r0=0.5
		#else:
		#1	r0=self.R
		gama=0.03
		x=round(float(1-round(np.exp(-gama*itr),2)),2)
		x1=r0*x
		x2=round(float(x1),2)
		self.R=x2
		

	def decrease_a(self):
		alpha=0.98
		a=self.A * alpha
		a1=round(float(a),2)
		self.A=a1

	def optimize(self):
		print("sdfd")
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
			print("iteration : ",(i+1))
			inc=0
			for p in self.bats:
				nbrs=nx.get_node_attributes(p,'n_distance')[1]
				gbst=nbrs[np.random.randint(len(nbrs))]

				#print("gbest ",gbst)
				#print("p",p.node[1]['pos'])
				self.updatevelocity(p,gbst)
				child=self.updatepos(p)
				#print("c",child.node[1]['pos'])
				#if(child==p):
				#	print True
				if(round(np.random.uniform(0,1),2)>self.R):
					child=self.eq_4(p,gbst) #equation 9 in your docx paper

				if(i<self.iteration*self.pm):
					self.turbulance_operation(child)
				
				self.Z1dZ2d(child)
				#self.update_neighborhood_solution(p,child)


				rd=round(np.random.uniform(0,1),2)

				fi = round(self.fitness(p),3)
				fnew = round(self.fitness(child),3)

				#print(rd,fi,fnew)

				if(rd<self.A and fi<fnew):
					p=child	# accept new solution in full population			
					self.increase_r(i) # increase pulse emission rate
					self.decrease_a() # decrease loudness
					
				#ns=nx.get_node_attributes(p,'n_distance')[1]
				self.update_neighborhood_solution(nbrs,child)
				self.update_reference_point()
				self.update_pbest(p,child)

			print(self.fitness(self.G))
			print(len(set([self.G.node[i]['pos'] for i in self.G])))
		print("\n\n**********************************************************")	
		print('\nThe script take {0} second '.format(np.round((time.time() - startTime),2)))
		print("\nFitness : ",self.fitness(self.G))
		print("\nNumber of Communites : ",len(set([self.G.node[i]['pos'] for i in self.G])))
		print("\nGlobal Best Position : ",[self.G.node[i]['pos'] for i in self.G])

		self.draw()
if __name__=='__main__':
    
	f=Bats()
	f.optimize()