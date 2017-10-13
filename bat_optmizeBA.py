import networkx as nx
import numpy as np
from collections import Counter, defaultdict
import time
#from sklearn import metrics
#from sklearn.metrics import normalized_mutual_info_score as NMI
import operator
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
		self.best_p = []
		self.global_best=[]
		self.pos_modularity = 0
		self.velocity = []
		self.G = nx.Graph()
		self.bats = []
		self.file_name = 'kara.txt'
		self.number_of_bats = 100
		self.modularity = 0
		self.iteration = 100
		self.R=0.5
		self.A=0.5

	def Input_Graph(self):
		temp=open(self.file_name,'r').read().split('\n')
		graph=[]
		for x in xrange(1,10):
		 	pass i in temp:
			t=[]
			for j in i.split():
				if(j):
					t.append(int(j))
			if(t):
				graph.append(tuple(t))
		self.G.add_edges_from(graph)
		j=1
		for i in self.G:
			self.G.node[i]={'pos':j,'w':[0,0],'n_dis':[]}
			j+=1


	def weight_init(self):
    	for i in self.G:
    		t=float(i)/(G.number_of_nodes()-1)
			self.G.node[i]['w'][0]=round(t,1)
			self.G.node[i]['w'][1]=1-round(t,1)




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
		return graph.copy()

	def updatevelocity(self,graph):
		fmin=0
		fmax=1
		v1=[]
		v2=[]
		v3=[]
		j=0
		p=self.best_positions
		#print('bf',self.best_positions)
		for i in graph:
			beta=round(np.random.uniform(0,1),2)
			f=fmin+(fmax-fmin)*beta
			v1.append(int((p[j]==graph.node[i]['pos']) and '0' or '1'))
			#bf=self.best_positions[j]
			#print('bf',self.best_positions[j])
			#print(self.best_positions[j])
			v2.append(round((self.velocity[j]+v1[j]*f),2))
			
			v3.append(round(float(1/round(1+round(float(np.exp(-v2[j])),2),2)),2))
			#print('af',self.best_positions[j])
			#print('bf',self.velocity)
			self.velocity[j]=int((round(np.random.uniform(0,1),2)<v3[j]) and '1' or '0')
			#print('Af',self.velocity)
			#af=self.best_positions[j]
			#print(bf==af)
			#if(round(np.random.uniform(0,1),2)< v3[j]):
			#	self.velocity[j]=1
			#else:
			#	self.velocity[j]=0
			j+=1
			
			#print('l',self.best_positions)
		#print('af',self.best_positions)	


	def bats_init(self): #based on random neighborhood
		#a=self.G.nodes()
		#l=np.random.randint(1,self.G.number_of_nodes(),self.G.number_of_nodes()).tolist()
		#self.newposition=l
		#self.best_positions=l
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
				#print(self.G.node[i]['pos'])
                                #t1=[]
				#t1.append(self.G.node[i]['pos'])
				#self.best_positions=t1
			#print('\n')
			self.bats.append(self.G.copy())
			#print(self.best_positions)
		
	def fitness(self,graph):
		m=graph.number_of_edges()
		l=1/(2*m)
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
			md+=(v1)/v2
			#print(md)
		
		
		n=G.number_of_nodes()
		k=community
		res=2(n-k)-md
		#print(round(md,2))
		return round(res,2)    	



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
			md+=(temp)/v2
			#print(md)
		#print(round(md,2))
		return round(md,2)



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
			md+=(v1)/v2
			#print(md)
		#print(round(md,2))
		return round(md,2)



	def zReference_init(self,graph):
		zr=[]
    	self.zReference_point=min(self.KKM(graph),self.RC(graph))
    	for i in graph:
    		



	def Euclidean_distance(self,graph):
		#round(np.sqrt(np.power((w2[0]-w1[0]),2) + np.power((w2[1]-w1[1]),2)),1)
		for i in graph:
			#graph.node[i]['w']
			nodes=graph.nodes()
			remain_node=nodes.remove(i)
			n_distance=[]
			for r in remain_node:
				w1=[]
				w2=[]
				w1.append(graph.nodes[i]['w'][0])
				w1.append(graph.nodes[i]['w'][1])

				w2.append(graph.nodes[i]['w'][0])
				w2.append(graph.nodes[i]['w'][1])

				#w1=[0,0]
				#w2=[0,0]
				#w1[0]=graph.nodes[i]['w'][0]
				#w1[1]=graph.nodes[i]['w'][1]

				#w2[0]=graph.nodes[i]['w'][0]
				#w2[0]=graph.nodes[i]['w'][1]
				n_distance.append(round(np.sqrt(np.power((w2[0]-w1[0]),2) + np.power((w2[1]-w1[1]),2)),1))
			graph.node[i].update({'n_dis':n_distance})




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
		node=list(graph.nodes())
		for i in graph:
			pos.append(graph.node[i]['pos'])
		new_pos=[]
		single=list(set(pos))
		for i in single:
			if(i in node):
				list(node).remove(i)
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


	def eq_4(self):
		new_graph=self.G.copy()           
		k=0
		av=np.mean(self.A)
		av1=round(av,2)
		for i in new_graph:
			e=round(np.random.uniform(0,1),2)
			if(e>av1):
                
				#new_graph.node[i]['pos']=self.best_positions[k]+np.mean(self.A)*e
				
				neighbor=new_graph.neighbors(i)
				temp=np.array(list(neighbor))
				#new_graph[i]['pos']=self.best_positions[k]+np.random.choice(neighbor)
				new_graph[i]['pos']=np.random.choice(temp)
			else:
				new_graph[i]['pos']=self.best_positions[k]
#			if(c>np.mean(self.A))
			k+=1
		return new_graph.copy()
		



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

	def best_bats(self):
		fit=-1
		for i in self.bats:
			t=self.fitness(i)
			self.bats_fitness_value.append(t)
			self.best_p.append([i.node[j]['pos'] for j in i])
			#self.bats_fitness_values.update{t:[i.node[j]['pos'] for j in i]}
			if(t>fit):
				fit=t
				#print(i.node)
				self.best_positions = [i.node[j]['pos'] for j in i]

	
	def turbulance_operation(self,graph):
		pm=1
		temp_graph=graph.copy()
		r=np.random.choice([0,1])
		for i in graph:
			if(r<pm):
				for j in graph.neighbors(i):
					temp_graph.node[j]['pos']=graph.node[i]['pos']
		self.G=temp_graph.copy()

	
	def optimize(self):
		self.__init__()
		startTime = time.time()
		self.Input_Graph()
		self.bats_init()
		#for i in self.bats:
		#	print(i.node)
		t=self.G.number_of_nodes()
		vel=[]
		for ll in range(t):
			vel.append(0)

		self.velocity=vel
		#print(self.velocity)
		#self.best_positions=vel;
		#self.best_p=vel
		#print('bf',self.velocity)
		self.best_bats()
		
		#print(self.velocity)
		#print(self.best_positions)
		#print('af',self.velocity)
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
				self.turbulance_operation(p)

				if(round(np.random.uniform(0,1),2)>self.R):
					#t2=self.eq_4()
					#y=self.rearrange(t2)
					y=self.eq_4()
				
				rd=round(np.random.uniform(0,1),2)
				fnew=self.fitness(y)
				#print(fnew)
				fi=self.bats_fitness_value[inc]
				bp=[]
				if(rd<self.A and fi<fnew):
					#print("p ",p.node)
					#print("y ",y.node)
					p=y.copy()					
					self.bats_fitness_value[inc] = fnew
					#self.best_p = 
					#self.bats.bats_fitness_values[(self.bats_fitness_values.keys())[inc]]=fnew

					self.best_positions = [y.node[nd]['pos'] for nd in y]
					self.best_p[inc] = [y.node[nd]['pos'] for nd in y] # this is an array of positions
					#self.best_p[inc] = []

					#if (self.pos_modularity<fnew):   # here comparison with old fitness if good then position vector and 
					#	self.pos_modularity=fnew     # modularity is update
						#fit = self.fitness(y)
						#print()
					#	self.best_p = [y.node[nd]['pos'] for nd in y]
						#print(self.best_p)

					#print(self.best_positions)
					#print("communites : ",len(set(self.best_positions)))
					#print(self.best_positions)
					self.increase_r(i)
					self.decrease_a()
				inc+=1
                        
			fmax=max(self.bats_fitness_value)
			#print(' max ',self.pos_modularity)
			#print(self.best_p)
			#fmax=max(self.bats_fitness_values.keys())
			#print(self.bats_fitness_value)
			#print(self.best_positions)
			pos=self.best_p[self.bats_fitness_value.index(fmax)] 
			self.positions=pos
			#print("communites : ",len(set(pos)))
			#print("Position : ",self.best_positions)
			print(fmax)
		
		print("\n\n**********************************************************")	
		print('\nThe script take {0} second '.format(np.round((time.time() - startTime),2)))
		print("\nModularity is : ",fmax)
		print("\nNumber of Communites : ",len(set(pos)))
		print("\nGlobal Best Position : ",pos)
		#nx.draw(self.G,node_color=[self.G.node[i]['pos'] for i in self.G])
		nx.draw(self.G, node_color=pos)
		plt.show()
if __name__=='__main__':
    
	f=Bats()
	f.optimize()
