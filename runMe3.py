from Controller import *
import numpy
import numpy as np
import matplotlib.pyplot as plt
import pylab
import mpl_toolkits.mplot3d.axes3d as p3
from time import time
###########
##@summary: This is to trace back the connectivity in the network saved in **.obj, which is generated during training by rumMe.py
##@author: Akihiro Eguchi
##Aug 13, 2013

nTrainings = 2001
ndim = 30
dim_stim = 10

controller = Controller(dim_stim,ndim)
resultFolderName = "results";
controller.resultFolderName = resultFolderName;
transOn = 1
singleColor = 1
alteringInput = 0#not fully implemented yet
itr = 2000

class InputStim:
    r=0
    g=0
    b=0
    def setRGB(self,r,g,b):
        self.r= r
        self.g= g
        self.b= b

inputStim = [[InputStim() for x in xrange(dim_stim)] for x in xrange(dim_stim)]

controller.loadWeightsAndDelays("Network_"+str(itr)+".obj",resultFolderName)
# controller.loadWeightsAndDelays("Network_"+str(itr)+".obj","results_2000_3rd_colorInputSaved")
controller.setLearningStates(0)#stop synaptic modifications
controller.variables.tstop = 300;
nThreads = 7;
controller.pc = h.ParallelContext()
controller.pc.nthread(nThreads)

if(alteringInput):
    controller.AlteringInputInit();


weightTemp_LtoL4 = []
weightTemp_LtoL4_s = []
weightTemp_LtoL4_e = []
weightTemp_LtoL4_se = []
weightTemp_C1toL4 = []
weightTemp_C1toL4_s = []
weightTemp_C1toL4_e = []
weightTemp_C1toL4_se = []
weightTemp_C2toL23 = []
weightTemp_C2toL23_e = []
weightTemp_C2toL23_se = []
weightTemp_C2toL23_s = []
weightTemp_L4toL23 = []
weightTemp_L4toL23_n = []
weightTemp_L4toL23_e = []
weightTemp_L4toL23_s = []
weightTemp_L4toL23_w = []
weightTemp_L23toL5 = []
weightTemp_L23toL5_n = []
weightTemp_L23toL5_e = []
weightTemp_L23toL5_s = []
weightTemp_L23toL5_w = []
weightTemp_L_r = []
weightTemp_L_g = []
weightTemp_C1_r = []
weightTemp_C1_g = []
weightTemp_C2_r = []
weightTemp_C2_g = []
weightTemp_C2_b = []

maxrgbWeights = 0.01
minrgbWeights = 0.0025

colPow=2;
weightPow1 = 2;
weightPow2 = 2;
weightPow3 = 1;

print len(controller.NetCons)
#for index in range(len(controller.NetCons_STDP_LtoL4)):
for index in range(ndim*ndim):
    weightTemp_LtoL4.append(pow(controller.NetCons_STDP_LtoL4[index].weight[0],weightPow1))
    weightTemp_LtoL4_s.append(pow(controller.NetCons_STDP_LtoL4[(index)*3+(ndim*ndim)].weight[0],weightPow1))
    weightTemp_LtoL4_e.append(pow(controller.NetCons_STDP_LtoL4[(index)*3+(ndim*ndim)+1].weight[0],weightPow1))
    weightTemp_LtoL4_se.append(pow(controller.NetCons_STDP_LtoL4[(index)*3+(ndim*ndim)+2].weight[0],weightPow1))
        
    weightTemp_C1toL4.append(pow(controller.NetCons_STDP_C1toL4[index].weight[0],weightPow1))
    weightTemp_C1toL4_s.append(pow(controller.NetCons_STDP_C1toL4[(index)*3+(ndim*ndim)].weight[0],weightPow1))
    weightTemp_C1toL4_e.append(pow(controller.NetCons_STDP_C1toL4[(index)*3+(ndim*ndim)+1].weight[0],weightPow1))
    weightTemp_C1toL4_se.append(pow(controller.NetCons_STDP_C1toL4[(index)*3+(ndim*ndim)+2].weight[0],weightPow1))
    
    weightTemp_C2toL23.append(pow(controller.NetCons_STDP_C2toL23[index].weight[0],weightPow1))
    weightTemp_C2toL23_s.append(pow(controller.NetCons_STDP_C2toL23[(index)*3+(ndim*ndim)].weight[0],weightPow1))
    weightTemp_C2toL23_e.append(pow(controller.NetCons_STDP_C2toL23[(index)*3+(ndim*ndim)+1].weight[0],weightPow1))
    weightTemp_C2toL23_se.append(pow(controller.NetCons_STDP_C2toL23[(index)*3+(ndim*ndim)+2].weight[0],weightPow1))
    
    weightTemp_L4toL23.append(pow(controller.NetCons_STDP_L4toL23[index].weight[0],weightPow2))
    weightTemp_L4toL23_s.append(pow(controller.NetCons_STDP_L4toL23[(index)*4+(ndim*ndim)].weight[0],weightPow2))
    weightTemp_L4toL23_e.append(pow(controller.NetCons_STDP_L4toL23[(index)*4+(ndim*ndim)+1].weight[0],weightPow2))
    weightTemp_L4toL23_n.append(pow(controller.NetCons_STDP_L4toL23[(index)*4+(ndim*ndim)+2].weight[0],weightPow2))
    weightTemp_L4toL23_w.append(pow(controller.NetCons_STDP_L4toL23[(index)*4+(ndim*ndim)+3].weight[0],weightPow2))

    weightTemp_L23toL5.append(pow(controller.NetCons_STDP_L23toL5[index].weight[0],weightPow2))
    weightTemp_L23toL5.append(pow(controller.NetCons_STDP_L23toL5[(index)*4+(ndim*ndim)].weight[0],weightPow2))
    weightTemp_L23toL5.append(pow(controller.NetCons_STDP_L23toL5[(index)*4+(ndim*ndim)+1].weight[0],weightPow2))
    weightTemp_L23toL5.append(pow(controller.NetCons_STDP_L23toL5[(index)*4+(ndim*ndim)+2].weight[0],weightPow2))
    weightTemp_L23toL5.append(pow(controller.NetCons_STDP_L23toL5[(index)*4+(ndim*ndim)+3].weight[0],weightPow2))

    weightTemp_L_r.append(pow(controller.NetCons[(index)*7].weight[0],weightPow3))
    weightTemp_L_g.append(pow(controller.NetCons[(index)*7+1].weight[0],weightPow3))
    weightTemp_C1_r.append(pow(controller.NetCons[(index)*7+2].weight[0],weightPow3))
    weightTemp_C1_g.append(pow(controller.NetCons[(index)*7+3].weight[0],weightPow3))
    weightTemp_C2_r.append(pow(controller.NetCons[(index)*7+4].weight[0],weightPow3))
    weightTemp_C2_g.append(pow(controller.NetCons[(index)*7+5].weight[0],weightPow3))
    weightTemp_C2_b.append(pow(controller.NetCons[(index)*7+6].weight[0],weightPow3))
    

totalWeights = [1 for x in xrange(ndim*ndim)]
decodedColors = numpy.zeros(3*ndim*ndim).reshape((ndim*ndim, 3))

#decodedColors = [[0 for x in xrange(3)] for x in xrange(ndim*ndim)] 


def calcL4toL23(L5i, L4i, preWeight):
    x = L4i%ndim
    y = L4i/ndim 
    i = L4i
    #C1
#     decodedColors[L5i,0] += preWeight*weightTemp_C1toL4[i]*pow(weightTemp_C1_r[i],colPow);
#     decodedColors[L5i,1] -= preWeight*weightTemp_C1toL4[i]*pow(weightTemp_C1_g[i],colPow);
    decodedColors[L5i,0] += preWeight * weightTemp_C1toL4[i] * max(0,(weightTemp_C1_r[i] - weightTemp_C1_g[i])) * pow(weightTemp_C1_r[i],colPow);
    decodedColors[L5i,1] += preWeight * weightTemp_C1toL4[i] * max(0,(weightTemp_C1_r[i] - weightTemp_C1_g[i])) * pow(maxrgbWeights+minrgbWeights-weightTemp_C1_g[i],colPow)/2;
    #0.0025 - 0.01

    #L
#     decodedColors[L5i,0] += preWeight*weightTemp_LtoL4[i]*pow(weightTemp_L_r[i],colPow);
#     decodedColors[L5i,1] += preWeight*weightTemp_LtoL4[i]*pow(weightTemp_L_g[i],colPow);
    decodedColors[L5i,0] += preWeight * weightTemp_LtoL4[i] * (weightTemp_L_r[i] + weightTemp_L_g[i]) * pow(weightTemp_L_r[i],colPow);
    decodedColors[L5i,1] += preWeight * weightTemp_LtoL4[i] * (weightTemp_L_r[i] + weightTemp_L_g[i]) * pow(weightTemp_L_g[i],colPow);
    
    x_e = (x+1)%ndim
    y_e = y
    i_e = y_e*ndim+x_e
#     decodedColors[L5i,0] += preWeight*weightTemp_C1toL4_e[i]*pow(weightTemp_C1_r[i_e],colPow);
#     decodedColors[L5i,1] -= preWeight*weightTemp_C1toL4_e[i]*pow(weightTemp_C1_g[i_e],colPow);
    decodedColors[L5i,0] += preWeight*weightTemp_C1toL4_e[i]*max(0,(weightTemp_C1_r[i_e] - weightTemp_C1_g[i_e]))*pow(weightTemp_C1_r[i_e],colPow);
    decodedColors[L5i,1] += preWeight*weightTemp_C1toL4_e[i]*max(0,(weightTemp_C1_r[i_e] - weightTemp_C1_g[i_e]))*pow(maxrgbWeights+minrgbWeights-weightTemp_C1_g[i_e],colPow)/2;

#     decodedColors[L5i,0] += preWeight*weightTemp_LtoL4_e[i]*pow(weightTemp_L_r[i_e],colPow);
#     decodedColors[L5i,1] += preWeight*weightTemp_LtoL4_e[i]*pow(weightTemp_L_g[i_e],colPow);
    decodedColors[L5i,0] += preWeight*weightTemp_LtoL4_e[i]*(weightTemp_L_r[i_e] + weightTemp_L_g[i_e])*pow(weightTemp_L_r[i_e],colPow);
    decodedColors[L5i,1] += preWeight*weightTemp_LtoL4_e[i]*(weightTemp_L_r[i_e] + weightTemp_L_g[i_e])*pow(weightTemp_L_g[i_e],colPow);
    
    x_s = x
    y_s = (y+1)%ndim
    i_s = y_s*ndim+x_s
#     decodedColors[L5i,0] += preWeight*weightTemp_C1toL4_s[i]*pow(weightTemp_C1_r[i_s],colPow);
#     decodedColors[L5i,1] -= preWeight*weightTemp_C1toL4_s[i]*pow(weightTemp_C1_g[i_s],colPow);
    decodedColors[L5i,0] += preWeight*weightTemp_C1toL4_s[i]*max(0,(weightTemp_C1_r[i_s]-weightTemp_C1_g[i_s]))*pow(weightTemp_C1_r[i_s],colPow);
    decodedColors[L5i,1] += preWeight*weightTemp_C1toL4_s[i]*max(0,(weightTemp_C1_r[i_s]-weightTemp_C1_g[i_s]))*pow(maxrgbWeights+minrgbWeights-weightTemp_C1_g[i_s],colPow)/2;
    
#     decodedColors[L5i,0] += preWeight*weightTemp_LtoL4_s[i]*pow(weightTemp_L_r[i_s],colPow);
#     decodedColors[L5i,1] += preWeight*weightTemp_LtoL4_s[i]*pow(weightTemp_L_g[i_s],colPow);
    decodedColors[L5i,0] += preWeight*weightTemp_LtoL4_s[i]*(weightTemp_L_r[i_s] + weightTemp_L_g[i_s])*pow(weightTemp_L_r[i_s],colPow);
    decodedColors[L5i,1] += preWeight*weightTemp_LtoL4_s[i]*(weightTemp_L_r[i_s] + weightTemp_L_g[i_s])*pow(weightTemp_L_g[i_s],colPow);

    
    x_se = (x+1)%ndim
    y_se = (y+1)%ndim
    i_se = y_se*ndim+x_se
#     decodedColors[L5i,0] += preWeight*weightTemp_C1toL4_se[i]*pow(weightTemp_C1_r[i_se],colPow);
#     decodedColors[L5i,1] -= preWeight*weightTemp_C1toL4_se[i]*pow(weightTemp_C1_g[i_se],colPow);
    decodedColors[L5i,0] += preWeight*weightTemp_C1toL4_se[i]*max(0,(weightTemp_C1_r[i_se] - weightTemp_C1_g[i_se]))*pow(weightTemp_C1_r[i_se],colPow);
    decodedColors[L5i,1] += preWeight*weightTemp_C1toL4_se[i]*max(0,(weightTemp_C1_r[i_se] - weightTemp_C1_g[i_se]))*pow(maxrgbWeights+minrgbWeights-weightTemp_C1_g[i_se],colPow)/2;
    
#     decodedColors[L5i,0] += preWeight*weightTemp_LtoL4_se[i]*pow(weightTemp_L_r[i_se],colPow);  
#     decodedColors[L5i,1] += preWeight*weightTemp_LtoL4_se[i]*pow(weightTemp_L_g[i_se],colPow);  
    decodedColors[L5i,0] += preWeight*weightTemp_LtoL4_se[i]*(weightTemp_L_r[i_se]+weightTemp_L_g[i_se])*pow(weightTemp_L_r[i_se],colPow);  
    decodedColors[L5i,1] += preWeight*weightTemp_LtoL4_se[i]*(weightTemp_L_r[i_se]+weightTemp_L_g[i_se])*pow(weightTemp_L_g[i_se],colPow);  
    
    


def calcL23toL5(L5i, L23i, preWeight):
    multi = 1.5
    x = L23i%ndim
    y = L23i/ndim
    i = L23i
    calcL4toL23(L5i, i, preWeight*weightTemp_L4toL23[i]*multi)
#     decodedColors[L5i,0] += preWeight*weightTemp_C2toL23[i]*pow(weightTemp_C2_r[i],colPow);
#     decodedColors[L5i,1] += preWeight*weightTemp_C2toL23[i]*pow(weightTemp_C2_g[i],colPow);
#     decodedColors[L5i,2] -= preWeight*weightTemp_C2toL23[i]*pow(weightTemp_C2_b[i]/2,colPow);
    decodedColors[L5i,0] += preWeight*weightTemp_C2toL23[i]*max(0,(weightTemp_C2_r[i]+weightTemp_C2_g[i]-weightTemp_C2_b[i]/2))*pow(weightTemp_C2_r[i],colPow);
    decodedColors[L5i,1] += preWeight*weightTemp_C2toL23[i]*max(0,(weightTemp_C2_r[i]+weightTemp_C2_g[i]-weightTemp_C2_b[i]/2))*pow(weightTemp_C2_g[i],colPow);
    decodedColors[L5i,2] += preWeight*weightTemp_C2toL23[i]*max(0,(weightTemp_C2_r[i]+weightTemp_C2_g[i]-weightTemp_C2_b[i]/2))*pow(maxrgbWeights+minrgbWeights-(weightTemp_C2_b[i]/2),colPow)/2;
    
    x_n = x 
    y_n = (y - 1)%ndim
    i_n = y_n*ndim+x_n
    calcL4toL23(L5i, i_n, preWeight*weightTemp_L4toL23_n[i]*multi)

    x_w = (x-1)%ndim
    y_w = y
    i_w = y_w*ndim+x_w
    calcL4toL23(L5i, i_w, preWeight*weightTemp_L4toL23_w[i]*multi)

    
    x_e = (x+1)%ndim
    y_e = y
    i_e = y_e*ndim+x_e
    calcL4toL23(L5i, i_e, preWeight*weightTemp_L4toL23_e[i]*multi)
#     decodedColors[L5i,0] += preWeight*weightTemp_C2toL23_e[i]*pow(weightTemp_C2_r[i_e],colPow);
#     decodedColors[L5i,1] += preWeight*weightTemp_C2toL23_e[i]*pow(weightTemp_C2_g[i_e],colPow);
#     decodedColors[L5i,2] -= preWeight*weightTemp_C2toL23_e[i]*pow(weightTemp_C2_b[i_e]/2,colPow);
    decodedColors[L5i,0] += preWeight*weightTemp_C2toL23_e[i]*max(0,(weightTemp_C2_r[i_e]+weightTemp_C2_g[i_e]-weightTemp_C2_b[i_e]/2))*pow(weightTemp_C2_r[i_e],colPow);
    decodedColors[L5i,1] += preWeight*weightTemp_C2toL23_e[i]*max(0,(weightTemp_C2_r[i_e]+weightTemp_C2_g[i_e]-weightTemp_C2_b[i_e]/2))*pow(weightTemp_C2_g[i_e],colPow);
    decodedColors[L5i,2] += preWeight*weightTemp_C2toL23_e[i]*max(0,(weightTemp_C2_r[i_e]+weightTemp_C2_g[i_e]-weightTemp_C2_b[i_e]/2))*pow(maxrgbWeights+minrgbWeights-(weightTemp_C2_b[i_e]/2),colPow)/2;
    
            
    x_s = x
    y_s = (y+1)%ndim
    i_s = y_s*ndim+x_s
    calcL4toL23(L5i, i_s, preWeight*weightTemp_L4toL23_s[i]*multi)
#     decodedColors[L5i,0] += preWeight*weightTemp_C2toL23_s[i]*pow(weightTemp_C2_r[i_s],colPow);
#     decodedColors[L5i,1] += preWeight*weightTemp_C2toL23_s[i]*pow(weightTemp_C2_g[i_s],colPow);
#     decodedColors[L5i,2] -= preWeight*weightTemp_C2toL23_s[i]*pow(weightTemp_C2_b[i_s]/2,colPow);
    decodedColors[L5i,0] += preWeight*weightTemp_C2toL23_s[i]*max(0,(weightTemp_C2_r[i_s]+weightTemp_C2_g[i_s]-weightTemp_C2_b[i_s]/2))*pow(weightTemp_C2_r[i_s],colPow);
    decodedColors[L5i,1] += preWeight*weightTemp_C2toL23_s[i]*max(0,(weightTemp_C2_r[i_s]+weightTemp_C2_g[i_s]-weightTemp_C2_b[i_s]/2))*pow(weightTemp_C2_g[i_s],colPow);
    decodedColors[L5i,2] += preWeight*weightTemp_C2toL23_s[i]*max(0,(weightTemp_C2_r[i_s]+weightTemp_C2_g[i_s]-weightTemp_C2_b[i_s]/2))*pow(maxrgbWeights+minrgbWeights-(weightTemp_C2_b[i_s]/2),colPow)/2;
        
    x_se = (x+1)%ndim
    y_se = (y+1)%ndim
    i_se = y_se*ndim+x_se
#     decodedColors[L5i,0] += preWeight*weightTemp_C2toL23_se[i]*pow(weightTemp_C2_r[i_se],colPow); 
#     decodedColors[L5i,1] += preWeight*weightTemp_C2toL23_se[i]*pow(weightTemp_C2_g[i_se],colPow); 
#     decodedColors[L5i,2] -= preWeight*weightTemp_C2toL23_se[i]*pow(weightTemp_C2_b[i_se]/2,colPow); 

    decodedColors[L5i,0] += preWeight*weightTemp_C2toL23_se[i]*max(0,(weightTemp_C2_r[i_se]+weightTemp_C2_g[i_se]-weightTemp_C2_b[i_se]/2))*pow(weightTemp_C2_r[i_se],colPow); 
    decodedColors[L5i,1] += preWeight*weightTemp_C2toL23_se[i]*max(0,(weightTemp_C2_r[i_se]+weightTemp_C2_g[i_se]-weightTemp_C2_b[i_se]/2))*pow(weightTemp_C2_g[i_se],colPow); 
    decodedColors[L5i,2] += preWeight*weightTemp_C2toL23_se[i]*max(0,(weightTemp_C2_r[i_se]+weightTemp_C2_g[i_se]-weightTemp_C2_b[i_se]/2))*pow(maxrgbWeights+minrgbWeights-(weightTemp_C2_b[i_se]/2),colPow)/2; 

for L5i in range(ndim*ndim):
    L5x = L5i%ndim
    L5y = L5i/ndim
    
    x = L5x
    y = L5y
    i = L5i 
    calcL23toL5(L5i, i, weightTemp_L23toL5[i])
    
    x_n = x 
    y_n = (y - 1)%ndim
    i_n = y_n*ndim+x_n
    calcL23toL5(L5i, i_n, weightTemp_L23toL5[i_n])
    
    x_w = (x-1)%ndim
    y_w = y
    i_w = y_w*ndim+x_w
    calcL23toL5(L5i, i_w, weightTemp_L23toL5[i_w])
    
    x_e = (x+1)%ndim
    y_e = y
    i_e = y_e*ndim+x_e
    calcL23toL5(L5i, i_e, weightTemp_L23toL5[i_e])
    
    x_s = x
    y_s = (y+1)%ndim
    i_s = y_s*ndim+x_s
    calcL23toL5(L5i, i_s, weightTemp_L23toL5[i_s])
    
#     decodedColors[L5i][0]/=counting[L5i][0];
#     decodedColors[L5i][1]/=counting[L5i][1];
#     decodedColors[L5i][2]/=counting[L5i][2];
#    decodedColors[L5i][1]/=100;
    
    print decodedColors[L5i]



maxr = np.max(decodedColors[:,0])
minr = np.min(decodedColors[:,0])
maxg = np.max(decodedColors[:,1])
ming = np.min(decodedColors[:,1])
maxb = np.max(decodedColors[:,2])
minb = np.min(decodedColors[:,2])

# max_ = max(maxr,max(maxg,maxb));
# maxr = max_
# maxg = max_
# maxb = max_

# print minr, ming, minb
print "maxr = " + str(maxr)
print "maxg = " + str(maxg)
print "maxb = " + str(maxb)
minr = 0
ming = 0
minb = 0


# maxr = np.max(controller.variables.trainingColors[:,0])
# minr = np.min(controller.variables.trainingColors[:,0])
# maxg = np.max(controller.variables.trainingColors[:,1])
# ming = np.min(controller.variables.trainingColors[:,1])
# maxb = np.max(controller.variables.trainingColors[:,2])
# minb = np.min(controller.variables.trainingColors[:,2])



fig=plt.figure()
ax = p3.Axes3D(fig)


for i in range(ndim*ndim):
#     ax.scatter((max(0,decodedColors[i,0])-minr)/(maxr-minr), (max(0,decodedColors[i,1])-ming)/(maxg-ming), (max(0,decodedColors[i,2])-minb)/(maxb-minb),s=80,c=((max(0,decodedColors[i,0])-minr)/(maxr-minr),(max(0,decodedColors[i,1])-ming)/(maxg-ming),(max(0,decodedColors[i,2])-minb)/(maxb-minb)))
    ax.scatter((decodedColors[i,0]-minr)/(maxr-minr), (decodedColors[i,1]-ming)/(maxg-ming), (decodedColors[i,2]-minb)/(maxb-minb),s=80,c=(math.sqrt((decodedColors[i,0]-minr)/(maxr-minr)),math.sqrt((decodedColors[i,1]-ming)/(maxg-ming)),math.sqrt((decodedColors[i,2]-minb)/(maxb-minb))))
#     r =(decodedColors[i,0]-minr)/(maxr-minr);
#     r_rev = (maxr - decodedColors[i,0])/(maxr-minr);
#     g = (decodedColors[i,1]-ming)/(maxg-ming);
#     g_rev = (maxg - decodedColors[i,1])/(maxg-ming);
#     b = (decodedColors[i,2]-minb)/(maxb-minb);
#     b_rev = (maxb-decodedColors[i,2])/(maxb-minb);
#     ax.scatter(r/(r+g+b), g/(r+g+b), b/(r+g+b),s=30*pow(r+g+b,1),c=(min(1,r/(r+g+b)*r*2+r_rev*0.5),min(1,g/(r+g+b)*g*2+g_rev*0.5),min(1,b/(r+g+b)*b*2+b_rev*0.5)))

# for i in range(2000):
# #     ax.scatter(decodedColors[i,0], decodedColors[i,1], decodedColors[i,2],c=((decodedColors[i,0]-minr)/(maxr-minr),(decodedColors[i,1]-ming)/(maxg-ming),(decodedColors[i,2]-minb)/(maxb-minb)/2))
# #    ax.scatter(decodedColors[i,0], decodedColors[i,1], decodedColors[i,2],c=((decodedColors[i,0]-minr)/(maxr-minr),(decodedColors[i,1]-ming)/(maxg-ming),(decodedColors[i,2]-minb)/(maxb-minb)))
#     ax.scatter(controller.variables.trainingColors[i,0], controller.variables.trainingColors[i,1], controller.variables.trainingColors[i,2],s=30,c=(controller.variables.trainingColors[i,0],controller.variables.trainingColors[i,1],controller.variables.trainingColors[i,2]))
#     print controller.variables.trainingColors[i,0], controller.variables.trainingColors[i,1], controller.variables.trainingColors[i,2]



ax.set_xlabel('Red')
ax.set_ylabel('Green')
ax.set_zlabel('Blue')
fig.add_axes(ax)
# ax.set_xlim([0, 5.0e-08])
# ax.set_ylim([0, 5.7e-08])
# ax.set_zlim([-5.6e-08, 0])
# ax.set_xlim([0, maxr])
# ax.set_ylim([0, maxg])
# ax.set_zlim([0, maxb])
ax.set_xlim([0, 1])
ax.set_ylim([0, 1])
ax.set_zlim([0, 1])

# ax.set_zlim([minb, 0])

plt.show()


#output
f = open("colors.txt", "a")
for i in range(ndim*ndim):
    f.write(str((max(0,decodedColors[i,0])-minr)/(maxr-minr)) + ',' + str((max(0,decodedColors[i,1])-ming)/(maxg-ming)) + ',' + str((max(0,decodedColors[i,2])-minb)/(maxb-minb)))
    if(i != ndim*ndim-1):
        f.write("\n")
f.close()

# f = open("colors_inputs.txt", "a")
# for i in range(2000):
#     f.write(str(controller.variables.trainingColors[i,0]) + ',' + str(controller.variables.trainingColors[i,1]) + ',' + str(controller.variables.trainingColors[i,2]))
#     if(i != 2000-1):
#         f.write("\n")
# f.close()

# f = open("colors_mod.txt", "a")
# for i in range(ndim*ndim):
#     r =(decodedColors[i,0]-minr)/(maxr-minr);
#     r_rev = (maxr - decodedColors[i,0])/(maxr-minr);
#     g = (decodedColors[i,1]-ming)/(maxg-ming);
#     g_rev = (maxg - decodedColors[i,1])/(maxg-ming);
#     b = (decodedColors[i,2]-minb)/(maxb-minb);
#     b_rev = (maxb-decodedColors[i,2])/(maxb-minb);
#     f.write(str(min(1,r/(r+g+b)*r*2)) + ',' + str(min(1,g/(r+g+b)*g*2)) + ',' + str(min(1,b/(r+g+b)*b*2)))
#     if(i != ndim*ndim-1):
#         f.write("\n")
# f.close()



# 
# if(singleColor==1):
#     fig1 = plt.gcf()
#     plt.clf()
#     for r in range(2):
#         for g in range(2):
#             for b in range(2):
#                 for y in range(dim_stim):
#                     for x in range(dim_stim):
#                         inputStim[y][x].setRGB(r,g,b)
#                 if(alteringInput):
#                     controller.setAlteringInput(inputStim, 0.25)
#                 else:
#                     controller.setInput(inputStim,0.8)
#                 controller.recordVols()
#                 #controller.recordChannelVols()
#                 controller.run()
#                                 
#                 controller.updateSpikeCount()
#                 #controller.outputFR(itr)
#                 #controller.saveSpikeDetails(r,g,b,itr)
#                 #controller.saveChannelSpikeDetails(r,g,b,itr)
#                 
#                 
#              
# 
#                 
#                 plt.subplot(8,4,r*4+g*2+b+1)
#                 plt.imshow(controller.spikeCount_L5,cmap=pylab.gray())
#                 plt.colorbar()
#                 
#                 plt.subplot(8,4,r*4+g*2+b+13)
#                 plt.imshow(controller.spikeCount_L23,cmap=pylab.gray())
#                 plt.colorbar()
#                 
#                 plt.subplot(8,4,r*4+g*2+b+25)
#                 plt.imshow(controller.spikeCount_L4,cmap=pylab.gray())
#                 plt.colorbar()
#                 
#                 
#                 
#                 
#                 
#                 if(transOn):
#                     controller.outputFR_trans(r,g,b,itr)
#                     
#                     #transformation: varies input with similar colours
#                     modVal = 0.05
#                     if r ==  0:
#                         rMod = r+modVal
#                     else:
#                         rMod = r-modVal
#                     if g ==  0:
#                         gMod = g+modVal
#                     else:
#                         gMod = g-modVal
#                     if b ==  0:
#                         bMod = b+modVal
#                     else:
#                         bMod = b-modVal
#                     
#                     for y in range(dim_stim):
#                         for x in range(dim_stim):
#                             inputStim[y][x].setRGB(rMod,g,b)
#                     controller.setInput(inputStim)
#                     controller.recordVols()
#                     controller.run()
#                     controller.updateSpikeCount()
#                     controller.outputFR_trans(r,g,b,itr)
#                     
#                     for y in range(dim_stim):
#                         for x in range(dim_stim):
#                             inputStim[y][x].setRGB(r,gMod,b)
#                     controller.setInput(inputStim)
#                     controller.recordVols()
#                     controller.run()
#                     controller.updateSpikeCount()
#                     controller.outputFR_trans(r,g,b,itr)
#                     
#                     for y in range(dim_stim):
#                         for x in range(dim_stim):
#                             inputStim[y][x].setRGB(r,g,bMod)                
#                     controller.setInput(inputStim)
#                     controller.recordVols()
#                     controller.run()
#                     controller.updateSpikeCount()
#                     controller.outputFR_trans(r,g,b,itr)
#                     
#                     for y in range(dim_stim):
#                         for x in range(dim_stim):
#                             inputStim[y][x].setRGB(rMod,gMod,b)
#                     controller.setInput(inputStim)
#                     controller.recordVols()
#                     controller.run()
#                     controller.updateSpikeCount()
#                     controller.outputFR_trans(r,g,b,itr)
#                     
#                     for y in range(dim_stim):
#                         for x in range(dim_stim):
#                             inputStim[y][x].setRGB(r,gMod,bMod)
#                     controller.setInput(inputStim)
#                     controller.recordVols()
#                     controller.run()
#                     controller.updateSpikeCount()
#                     controller.outputFR_trans(r,g,b,itr)
#     
#                     for y in range(dim_stim):
#                         for x in range(dim_stim):
#                             inputStim[y][x].setRGB(rMod,g,bMod)                
#                     controller.setInput(inputStim)
#                     controller.recordVols()
#                     controller.run()
#                     controller.updateSpikeCount()
#                     controller.outputFR_trans(r,g,b,itr)
#     
#                     for y in range(dim_stim):
#                         for x in range(dim_stim):
#                             inputStim[y][x].setRGB(rMod,gMod,bMod)                
#                     controller.setInput(inputStim)
#                     controller.recordVols()
#                     controller.run()
#                     controller.updateSpikeCount()
#                     controller.outputFR_trans(r,g,b,itr)
#     fig1.savefig("results/"+str(itr),dpi=100)
#     plt.show()
# else:
# #     controller.variables.tstop = 300
#     fig1 = plt.gcf()
#     plt.clf()
#     b = 0
#     for r in range(2):
#         for g in range(2):
#             if (r==g):
#                 continue
#             for y in range(dim_stim):
#                 for x in range(dim_stim):
#                     if y>dim_stim/3 and y<dim_stim*2/3 and x>dim_stim/3 and x>dim_stim*2/3:
#                         inputStim[y][x].setRGB(1-r,1-g,b)
#                     else:
#                         inputStim[y][x].setRGB(r,g,b)
#             controller.setInput(inputStim,0.3)
#             controller.recordVols()
#             controller.run()
#             controller.updateSpikeCount()
#             controller.outputFR(itr)    
#             
#             plt.subplot(2,1,r+1)
#             plt.imshow(controller.spikeCount_L4,cmap=pylab.gray())
#             plt.colorbar()
#             
#             #controller.drawGraph()
#             controller.saveSpikeDetails(r,g,b,111110);
#             
#     fig1.savefig("results/multiColTest300_normal"+str(itr),dpi=100)
            
controller.setLearningStates(1)#start synaptic modifications



    
# raw_input("Press Enter to exit...")