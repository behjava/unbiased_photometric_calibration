#This code uses the selected star sample for the iterative 3-parameter calibration: r_cal=aL+b(g-r)+c. The calibration starts with a first estimation of a,b, and c and rejects outliers from a given standard deviation (2 or 3). The fitting and the outlier rejection continues until there is no outlier. The final fitting is performed with higher precision (i.e. finer parameter grids). The final results of the calibration will be written in the calib_result_final.txt
#Behnam Javanmardi (Nov. 2016)

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import pylab


mag_limit= 19 #maxmimum desired magnitude limit of stars
sigma_limit= 2 #the level above which a star is considered as an outlier e.g. 2 means stars that deviate more than 2\sigma from the best fit of each step will be rejected and will not be used in the next step  


# importing the selected stars
num,ra,dec,l,l_err,r,g_r=np.genfromtxt("number_RA_DEC_L_L-err_r_g-r.txt",unpack=True)

def rcalc(a,b,c,flux,color):
    return a*flux+b*(color)+c

print "first estimate:"
# estimating the intercept (i.e. the parameter c)
slope_rl, intercept_rl, r_value, p_value, std_err = stats.linregress(l,r)

# estimating the color correction (i.e. the parameter b)

r_estimate=np.zeros(len(num))
r_estimate=slope_rl*l+intercept_rl

slope_gr, intercept_gr, r_value2, p_value2, std_err2 = stats.linregress((r-r_estimate),g_r)

    
#parameter_range    
a_min=slope_rl-0.6
a_max=slope_rl+0.6
b_min=slope_gr-1
b_max=slope_gr+1
c_min=(intercept_rl+intercept_gr)-5
c_max=(intercept_rl+intercept_gr)+5

#bin size
bin1=0.1


a_num=int((a_max-a_min)/bin1)
b_num=int((b_max-b_min)/bin1)+1
c_num=int((c_max-c_min)/bin1)

  
chi=np.zeros((a_num,b_num,c_num))





rc=np.zeros(len(num))
for i in range(len(num)):
    rc[i]=rcalc(slope_rl,slope_gr,intercept_rl+intercept_gr,l[i],g_r[i])

stdv1=np.std(r-rc)

print "initial standard deviation of the residuals = " +str(stdv1)
print "end of first estimate"

rloop=r
grloop=g_r
numloop=num
raloop=ra
decloop=dec
lloop=l
lerrloop=l_err

print "outlier rejection:"

stdv2=0
dim=len(num)
while stdv2<stdv1:
	stdv1=np.std(rloop-rc)
	n=0
	for i in range(dim):
	    if abs(rloop[i]-rc[i])<sigma_limit*stdv1 and rloop[i]<mag_limit:
	        n=n+1

		numnew=np.zeros(n)
	ranew=np.zeros(n)
	decnew=np.zeros(n)
	lnew=np.zeros(n)
	lerrnew=np.zeros(n)
	rnew=np.zeros(n)
	grnew=np.zeros(n)

	n=0
	for i in range(dim):
	    if abs(rloop[i]-rc[i])<sigma_limit*stdv1 and rloop[i]<mag_limit:
	        n=n+1
	        rnew[n-1]=rloop[i]
	        grnew[n-1]=grloop[i]
	        numnew[n-1]=numloop[i]
	        ranew[n-1]=raloop[i]
	        decnew[n-1]=decloop[i]
	        lnew[n-1]=lloop[i]
 	       	lerrnew[n-1]=lerrloop[i]

	for j in range(a_num):
	    a=a_min+(j*bin1)
	    for k in range(b_num):
	        b=(b_min+k*bin1)
	        for m in range(c_num):
	            c=c_min+(m*bin1)
 	            for i in range(len(numnew)):
	                chi[j,k,m]=chi[j,k,m]+(rnew[i]-rcalc(a,b,c,lnew[i],grnew[i]))**2/lerrnew[i]**2




	aa,bb,cc=np.unravel_index(chi.argmin(), chi.shape)
	aa=(aa*bin1)+a_min
	bb=((bb*bin1)+b_min)
	cc=(cc*bin1)+c_min


	rc=np.zeros(len(numnew))
	for i in range(len(numnew)):
 	   rc[i]=rcalc(aa,bb,cc,lnew[i],grnew[i])

	stdv2=np.std(rnew-rc)
	print "standard deviation of the residuals = " +str(stdv2)
	dim=len(numnew)
	rloop=rnew
	grloop=grnew
	numloop=numnew
	raloop=ranew
	decloop=decnew
	lloop=lnew
	lerrloop=lerrnew

print "end of outlier rejection"

f = open("calib_results_after_outlier_rejection_new_code_v2.txt", "w")
f.write('#a' + "\n")
f.write('#b' + "\n")
f.write('#c' + "\n")
f.write(str(aa) + "\n")
f.write(str(bb) + "\n")
f.write(str(cc) + "\n")

f.close()   

np.savetxt('after_outlier_rejection.txt',np.column_stack((numloop,raloop,decloop,lloop,lerrloop,rloop,grloop)),delimiter='    ') 


print "final calibration with higher precision:"
print "please wait ..."

#new parameter_range (smaller) 
a_min=aa-0.15
a_max=aa+0.15
b_min=bb-0.15
b_max=bb+0.15
c_min=cc-2
c_max=cc+2

#new bin size (smaller)
bin1=0.01


a_num=int((a_max-a_min)/bin1)
b_num=int((b_max-b_min)/bin1)+1
c_num=int((c_max-c_min)/bin1)

  
chi=np.zeros((a_num,b_num,c_num))

for j in range(a_num):
    a=a_min+(j*bin1)
    for k in range(b_num):
        b=(b_min+k*bin1)
        for m in range(c_num):
            c=c_min+(m*bin1)
            for i in range(len(numloop)):
                chi[j,k,m]=chi[j,k,m]+(rloop[i]-rcalc(a,b,c,lloop[i],grloop[i]))**2/lerrloop[i]**2
            
            



aa,bb,cc=np.unravel_index(chi.argmin(), chi.shape)
aa=(aa*bin1)+a_min
bb=((bb*bin1)+b_min)
cc=(cc*bin1)+c_min

rc=np.zeros(len(numloop))
for i in range(len(numloop)):
    rc[i]=rcalc(aa,bb,cc,lloop[i],grloop[i])

stdv=np.std(rloop-rc)



f = open("calib_result_final.txt", "w")
f.write('a = ' +str(aa) + "\n")
f.write('b = ' +str(bb) +"\n")
f.write('c = ' + str(cc) +"\n")
f.write('stdv = ' +str(stdv) +"\n")
f.write('final number of stars = ' +str(len(numloop)) +"\n")

f.close() 

print "end of calibration"

print "final standard deviation of the residuals = " +str(stdv)



  
