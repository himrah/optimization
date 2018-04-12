import networkx as nx
import numpy as np
from collections import Counter, defaultdict
import time
from operator import itemgetter
import random
import matplotlib.pyplot as plt

class Bats:
	def __init__(self):
		self.ns_size=40
		self.pm = 0.15  #probability mutation
		self.velocity = []
		self.G = nx.Graph()
		self.bats = []
		self.file_name = 'dolphin.txt'
		self.number_of_bats = 100
		self.nbrs=[]
		self.modularity = 0
		self.iteration = 10
		self.Z1d = 0
		self.Z2d = 0
		self.R = self.A = [0.5]*self.number_of_bats
		#self.R=0.5
		#self.A=0.5


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
		temp={}
		dic={}
		index=0
		for i in self.bats:
			temp.update({index:self.KKM(i)})
			index+=1
			#rs.append(self.KKM(i))

		rs = min(temp.items(),key=itemgetter(1))

		dic.update({"index":rs[0],"KKM":rs[1],"RC":self.RC(self.bats[rs[0]])})
		
		return dic
		#return min(rs)

	def best_rc(self):
		rs=[]
		temp={}
		dic={}
		index=0
		for i in self.bats:
			temp.update({index:self.RC(i)})
			index+=1
		rs = min(temp.items(),key=itemgetter(1))
		dic.update({"index":rs[0],"RC":rs[1],"KKM":self.KKM(self.bats[rs[0]])})		
		return dic

	def z_reference(self):# ideal point 
		self.zstr1=self.best_kkm()
		self.zstr2=self.best_rc()


	def update_reference_point(self):

		if self.Z1d < self.zstr1['KKM']:
			self.zstr1['KKM'] = self.Z1d
			self.zstr1['RC'] = self.Z2d

		if self.Z2d < self.zstr2['RC']:
			self.zstr2['RC'] = self.Z2d
			self.zstr2['KKM'] = self.Z1d


	def Z1dZ2d(self,child):
		self.Z1d=(self.KKM(child))
		self.Z2d=(self.RC(child))


	def update_pbest(self,graph,child):
		nob=2
		sum=0
		counter=0
		better=1
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
			v1.append(int((self.bats[gbest].node[i]['pos']==graph.node[i]['pos']) and '0' or '1'))
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
		for graph in self.bats:
			dic={}
			f=self.bats[i]
			index=0
			for g in self.bats:
				if g != f:
					A = graph.node[graph.nodes()[0]]['weight'][0]
					B = graph.node[graph.nodes()[0]]['weight'][1]
					a = g.node[g.nodes()[0]]['weight'][0]
					b = g.node[g.nodes()[0]]['weight'][1]			
					ed = round(np.sqrt(np.power((a-A),2) + np.power((b-B),2)),3)
					dic.update({index:ed})
				index+=1		
			sorted_dic = sorted(dic.items(), key=itemgetter(1))
			neigbhorhood_set=[s[0] for s in sorted_dic[:self.ns_size]]
			self.nbrs.append(tuple(neigbhorhood_set))
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
		child_pos=nx.get_node_attributes(child,'pos')
		nx.set_node_attributes(self.bats[neighbor],'pos',child_pos)
		#self.G=self.bats[neighbor].copy()


	def update_neighborhood_solution(self,ns,child):
		for i in ns:
			point = [self.zstr1['KKM'],self.zstr2['RC']]
			fun = [self.KKM(self.bats[i]),self.RC(self.bats[i])]
			weight = nx.get_node_attributes(self.bats[i],'weight')[1]
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


	def draw(self):
		bats=list((set(self.bats)))		
		M = [self.KKM(i) for i in bats]
		N = [self.RC(i) for i in bats]
		F = [self.fitness(i) for i in bats]
		print(F)
		ind = F.index(max(F))

		print("\nFitness : ",self.fitness(bats[ind]))
		print("\nNumber of Communites : ",len(set([bats[ind].node[i]['pos'] for i in bats[ind]])))
		print("\nGlobal Best Position : ",[bats[ind].node[i]['pos'] for i in bats[ind]])

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
		ax.scatter(M, N)
		for i, txt in enumerate(F):
		    ax.annotate(txt, (M[i],N[i]))
		plt.figure(2)
		nx.draw(bats[ind],node_color=[bats[ind].node[i]['pos'] for i in bats[ind]])
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
				graph.node[i]['pos']=self.bats[gbest].node[i]['pos']
#			if(c>np.mean(self.A))
		return graph

	def increase_r(self,index):
		#if(itr==0):
		r0=0.5
		#else:
		#1	r0=self.R
		gama=0.03
		x=round(float(1-round(np.exp(-gama*index),2)),2)
		x1=r0*x
		x2=round(float(x1),2)
		self.R[index]=x2
		

	def decrease_a(self,index):
		alpha=0.98
		a=self.A[index] * alpha
		a1=round(float(a),2)
		self.A[index]=a1

	def optimize(self):
		self.__init__()
		startTime = time.time()
		self.Input_Graph()
		self.bats_init()
		t=self.G.number_of_nodes()
		self.velocity = [0]*self.G.number_of_nodes() # put 0 in velocity t(number of nodes) times
		self.weight_init()
		self.pbest_init()
		self.z_reference()
		self.Euclidean_distance()
		
		for i in range(self.iteration):
			print("iteration : ",(i+1))
			index=0
			for p in self.bats:
				nbrs = self.nbrs[index]
				gbst = np.random.choice(nbrs)
				self.updatevelocity(p,gbst)
				child=self.updatepos(p)
				if(round(np.random.uniform(0,1),2)>self.R[index]):
					child=self.eq_4(p,gbst) #equation 9 in your docx paper

				if(i<self.iteration*self.pm):
					self.turbulance_operation(child)
				
				self.Z1dZ2d(child)

				rd=round(np.random.uniform(0,1),2)

				fi = round(self.fitness(p),3)
				fnew = round(self.fitness(child),3)

				if(rd<self.A[index] and fi<fnew):
					p=child	# accept new solution in full population
					self.G=child.copy()
					self.increase_r(i) # increase pulse emission rate
					self.decrease_a(i) # decrease loudness
					
				self.update_neighborhood_solution(nbrs,child)
				self.update_reference_point()
				self.update_pbest(p,child)
				index += 1

		self.draw()
if __name__=='__main__':
    
	f=Bats()
	f.optimize()