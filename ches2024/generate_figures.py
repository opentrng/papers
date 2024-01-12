# -*- coding: utf-8 -*-

import colorednoise as cn
import numpy as np
import matplotlib.pyplot as plt
import math


#Least squares normalized Error Regression algorithm from [Grantham & Bailey 2006]

def LSNE(x,y,a):
    H=np.zeros((np.size(x),np.size(a)))
    for r in range(np.size(x)):
        for c in range (np.size(a)):
            H[r,c]=x[r]**a[c]
    
    C=np.zeros((np.size(a),1))
    F=np.zeros((np.size(a),np.size(x)))
    for r in range(np.size(x)):
        for c in range (np.size(a)):
            F[c,r]=x[r]**a[c]/(y[r]**2)
    
    C=np.matmul(np.matmul(np.linalg.inv(np.matmul(F,H)),F),y)
    return C

#proportionality coefficients between colored noise and allan variance 

#input values
f = 1               # signal frequency
T = 100000000       # samples
n = int(T * f)      # total number of samples
#generate perfect series of t
t = np.linspace(0, T, n, endpoint=True)
samples = np.size(t)  # number of samples to generate (time series extension)

#generate noisy time series
di_thermal = cn.powerlaw_psd_gaussian(0, samples)
di_flicker = cn.powerlaw_psd_gaussian(1, samples)

dti_norm_th=1/f+di_thermal
dti_norm_fl=1/f+di_flicker

ti_norm_th=np.cumsum(dti_norm_th)
ti_norm_fl=np.cumsum(dti_norm_fl)

#Calculate Allan variance and proportionality factors
NoSamples=1000
Nvalues=np.array([1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,
                  200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,
                  12000,14000,16000,18000,20000,22000,24000,26000,28000,30000,32000,34000,36000,38000,40000,
                  42000,44000,46000,48000,50000,52000,54000,56000,58000,60000,62000,64000,66000,68000,70000,
                  72000,74000,76000,78000,80000,82000,84000,86000,88000,90000,92000,94000,96000,98000,100000,
                  120000,140000,160000,180000,200000,220000,240000,260000,280000,300000,320000,340000,360000,380000,400000,
                  420000,440000,460000,480000,500000,520000,540000,560000,580000,600000,620000,640000,660000,680000,700000,
                  720000,740000,760000,780000,800000,820000,840000,860000,880000,900000,920000,940000,960000,980000,1000000])

sN_norm_th=np.array([])
sN_norm_fl=np.array([])
for N in Nvalues:    
    sN_norm_th=np.append(sN_norm_th,np.var(np.diff(np.diff(ti_norm_th[0:((NoSamples+1)*N-1):N]))))
    sN_norm_fl=np.append(sN_norm_fl,np.var(np.diff(np.diff(ti_norm_fl[0:((NoSamples+1)*N-1):N]))))

poly_norm_th=np.zeros(3)
poly_norm_fl=np.zeros(3)
poly_norm_th[1]=LSNE(Nvalues,sN_norm_th,[1])
poly_norm_fl[0]=LSNE(Nvalues,sN_norm_fl,[2])

plt.figure()
p1=plt.plot(Nvalues,sN_norm_th,marker='o')
p2=plt.plot(Nvalues,np.polyval(poly_norm_th,Nvalues),linestyle='--',color='red')
p1=plt.plot(Nvalues,sN_norm_fl,marker='o')
p2=plt.plot(Nvalues,np.polyval(poly_norm_fl,Nvalues),linestyle='--',color='red')

plt.legend(['thermal noise','fit thermal noise','flicker noise','fit flicker noise'])
plt.title('Allan variance on normalized noise')
plt.yscale('log')
plt.xscale('log')
plt.ylabel('Allan variance (s²)')
plt.xlabel("accumulation time (s)")
plt.savefig("allan_variance_normalized_noise.pdf") 

factor_th=poly_norm_th[1]
factor_fl=poly_norm_fl[0]

#Allan variance from measurement results (all in seconds)

#parameters of the measured curves 
Nvalues_measurement=np.array([1,2,3,4,5,6,7,8,9,
         10,20,30,40,50,60,70,80,90,100,
         200,300,400,500,600,700,800,900,1000])

#sN_measurement_factor, where factor means sampling to signal frequency ratio
sN_measurement_10000=np.array([2.841881638875585970e-19,3.998914558142675538e-19,5.203638197797840649e-19,7.290925681926167425e-19,8.008383643442806245e-19,9.631163923195517231e-19,1.120585198572292624e-18,1.304127342183394220e-18,1.500829393475229912e-18,1.655452639647892064e-18,4.049978021883548450e-18,7.441104101810376871e-18,1.316175405657757989e-17,1.788307220858829079e-17,2.361016109317395456e-17,3.197381247561899276e-17,4.291411740177915851e-17,5.436051000894076872e-17,6.668424418418818530e-17,2.458698391379283347e-16,5.551362852071704166e-16,1.004404659621934921e-15,1.635147505752994960e-15,2.434490890555127810e-15,3.593861020561086983e-15,5.134622857670429137e-15,7.000510948829927250e-15,8.739506776149941947e-15])
sN_measurement_2000=np.array([8.437025796615816395e-19,8.678350139259610602e-19,9.385426716773394948e-19,1.010517144786125692e-18,1.105389891575551543e-18,1.118808432898057220e-18,1.287603570988638769e-18,1.473913585302100671e-18,1.581966383604705806e-18,1.555656927337188604e-18,3.719629976064230479e-18,6.064627158223283126e-18,1.031687957877347539e-17,1.504844478703817207e-17,2.025766292844892386e-17,2.658620549371996842e-17,3.479420479960331150e-17,4.355211066848127579e-17,5.224140164298456278e-17,1.927034701341616697e-16,4.059059573917002554e-16,7.524092005066064525e-16,1.094493526145383158e-15,1.579419732901781831e-15,2.141935456752761646e-15,2.753608857015976641e-15,3.520707189641217071e-15,4.509443032025994794e-15])
sN_measurement_1000=np.array([1.614567526560164307e-18,2.066082265667357982e-18,2.500758725209207432e-18,2.838185083926241530e-18,2.755477763410453748e-18,2.681696904337584853e-18,3.037168141780383653e-18,2.967622664573631627e-18,3.118397085625135319e-18,3.374421637537298157e-18,5.171483621849479331e-18,8.125305892125413974e-18,1.138795656385816726e-17,1.588482221812833765e-17,2.262286972615949923e-17,2.886935604198830628e-17,3.488111489066605599e-17,4.078875271037570251e-17,5.324648896363061424e-17,2.013695221788571103e-16,4.292710066381756614e-16,7.670710582178261159e-16,1.172811627327798097e-15,1.705150560322697365e-15,2.276542295780546796e-15,2.932770761878471074e-15,3.753874915358839711e-15,4.663421822689921046e-15])
sN_measurement_500=np.array([3.706626201330611575e-18,5.232629272493106228e-18,4.880121396187617639e-18,5.436893204122116721e-18,5.585160202593114913e-18,6.050303556291688917e-18,6.428318584206542995e-18,7.611497852281875185e-18,8.014571949226467261e-18,8.183992640421792814e-18,1.048163616937108248e-17,1.427618445604298192e-17,1.829016468101865577e-17,2.176188901480090923e-17,2.779676146726926066e-17,3.523814200601273944e-17,4.122237801879325085e-17,4.947307267137520084e-17,5.829027192989253545e-17,1.856806061617351235e-16,3.770124707258463369e-16,6.824902793890548602e-16,1.030645787148293725e-15,1.461622626488669270e-15,2.005526374678629834e-15,2.573165823851822343e-15,3.210714596542315390e-15,4.093495671857174945e-15])
sN_measurement_200=np.array([3.621648963062981880e-17,7.269049224476777796e-17,9.097111087038338063e-17,5.663430420687169913e-17,2.344013490694108542e-17,1.682567215914899865e-17,5.017691283513182533e-17,8.075539568415082003e-17,7.896166568822564952e-17,4.931002759914639435e-17,7.263969171373882451e-17,5.905964942491188892e-17,5.923000986620523504e-17,6.626984127700280151e-17,7.174119452585240670e-17,7.457627119415067114e-17,8.681308699480006398e-17,8.698688678007649339e-17,1.025049095962463135e-16,2.447093399371772445e-16,4.206652209474059361e-16,6.994950518040046192e-16,1.039696051670686735e-15,1.467878379649995627e-15,1.969089276438463748e-15,2.586355372405730066e-15,3.092217222701692496e-15,4.000864923486900880e-15])
sN_measurement_100=np.array([3.279715719142902484e-16,1.432230143464601097e-16,1.884671445432236317e-16,2.834951456298521928e-16,4.451939291741827113e-17,3.639199072680228999e-16,1.111504424780493150e-16,2.104313311952178814e-16,2.677592310557717233e-16,6.550137994622254094e-17,1.156069364189063152e-16,1.520078354562776645e-16,2.021717670264812539e-16,2.142857142805532186e-16,2.177110467657960871e-16,2.249252243101792112e-16,2.289706297639587519e-16,2.342338334366249678e-16,2.581162324561040611e-16,3.657603222384619355e-16,5.737866643728368550e-16,8.718466198100217582e-16,1.311110702880911143e-15,1.734908682679121778e-15,2.124030201191050877e-15,2.892016733202349490e-15,3.507837976056911296e-15,4.241455607375863763e-15])
avgT_RO_measurement=2e-6


poly_VA_norm_2000=LSNE(Nvalues_measurement*avgT_RO_measurement,sN_measurement_2000,[2,1,0])

avgT_RO_emulator=2e-9
dti_emulator_2000_1=avgT_RO_emulator*np.ones(samples)+di_thermal*np.sqrt(poly_VA_norm_2000[1]*avgT_RO_emulator/factor_th)+di_flicker*np.sqrt(poly_VA_norm_2000[0]*(avgT_RO_emulator**2)/factor_fl)

ti_emulator_2000_1=np.cumsum(dti_emulator_2000_1)

#apply artificial sampling
def sampling(matrix,sampling_step):
    sampled_matrix=sampling_step*np.round(matrix/sampling_step)
    return sampled_matrix

#Sampling step 1ns 
dti_sampling=1e-9

ti_emulator_2000_1=sampling(ti_emulator_2000_1,dti_sampling)

sN_emulator_2000_1=np.array([])
for N in Nvalues:    
    sN_emulator_2000_1=np.append(sN_emulator_2000_1,np.var(np.diff(np.diff(ti_emulator_2000_1[0:((NoSamples+1)*N-1):N]))))
    
poly_emulator_2000_1=LSNE(Nvalues[25:]*avgT_RO_emulator,sN_emulator_2000_1[25:],[2,1,0])

plt.figure()
p1=plt.plot(avgT_RO_measurement*Nvalues_measurement,sN_measurement_2000,marker='o')
p2=plt.plot(avgT_RO_measurement*Nvalues_measurement,np.polyval(poly_VA_norm_2000,avgT_RO_measurement*Nvalues_measurement),linestyle='--',color='red')
p1=plt.plot(Nvalues[25:]*avgT_RO_emulator,sN_emulator_2000_1[25:],marker='o')
p1=plt.plot(Nvalues[25:]*avgT_RO_emulator,np.polyval(poly_emulator_2000_1,Nvalues[25:]*avgT_RO_emulator),linestyle='--',color='red')

plt.legend(['measurement','fit measurement','emulator','fit emulator'])
plt.title('Allan variance measurement vs emulator')
plt.yscale('log')
plt.xscale('log')
plt.ylabel('Allan variance (s²)')
plt.xlabel("time (s)")
plt.savefig("allan_variance_measurements_vs_emulator.pdf") 

print('poly 1GS/s measurement     '+str(poly_VA_norm_2000))
print('poly 1GS/s emulator   '+str(poly_emulator_2000_1))
print('difference:                   '+str(round(100*abs(poly_VA_norm_2000[2]-poly_emulator_2000_1[2])/poly_VA_norm_2000[2],2))+' %        '
                                      +str(round(100*abs(poly_VA_norm_2000[1]-poly_emulator_2000_1[1])/poly_VA_norm_2000[1],2))+' %        '
                                      +str(round(100*abs(poly_VA_norm_2000[0]-poly_emulator_2000_1[0])/poly_VA_norm_2000[0],2))+' %      ')


    
################################################################################################
########################              FPGA results          ####################################
################################################################################################

counter_full=np.loadtxt('fpga_counter.txt')[:110000]

plt.figure()
p1=plt.hist(counter_full)

plt.legend(['mean value='+str(np.mean(counter_full))])
plt.title('Histogram counter')
plt.yscale('linear')
plt.xscale('linear')
plt.ylabel('samples')
plt.xlabel("counter")
plt.savefig("fpga_histogram_counter.pdf")

Nmin=100
N_bar=np.mean(counter_full)
NoSamples=1000 

Nvalues_counter_i=(np.linspace(1,np.int_(np.size(counter_full)/NoSamples))).astype(int)
Nvalues_counter=(int(N_bar)*np.linspace(1,np.int_(np.size(counter_full)/NoSamples))).astype(int)
      


sN_counter=np.array([])
for N in Nvalues_counter_i:    
    sN_counter=np.append(sN_counter,np.var(np.diff(np.convolve(counter_full,np.ones(N, dtype=int),mode='valid')[0:((NoSamples+1)*N-1):N])))

poly_counter=LSNE(Nvalues_counter,sN_counter,[2,1,0])

plt.figure()
p1=plt.plot(Nvalues_counter,sN_counter,marker='o')
p2=plt.plot(Nvalues_counter,np.polyval(poly_counter,Nvalues_counter),color='red')

plt.legend(['Allan variance counter','fit'])
plt.title('Allan variance counter')
plt.yscale('log')
plt.xscale('log')
plt.ylabel('variance counter [1]')
plt.xlabel("N*k ")
plt.savefig("allan_variance_counter.pdf")

samples = 10000000  # number of samples to generate (time series extension)

di_thermal = cn.powerlaw_psd_gaussian(0, samples)
di_flicker = cn.powerlaw_psd_gaussian(1, samples)

a2=poly_counter[0]
a1=poly_counter[1]

factor_th=2
factor_fl=0.135

avgT=N_bar

dcounter_emulator_baseline=avgT*np.ones(samples)+di_thermal*np.sqrt(a1*avgT/factor_th)+di_flicker*np.sqrt(a2*(avgT**2)/factor_fl)

counter_emulator_baseline_quantized=np.round(np.diff(np.cumsum(dcounter_emulator_baseline)))

NoSamples=1000
sN_counter_emulated=np.array([])
for N in Nvalues_counter_i:    
    sN_counter_emulated=np.append(sN_counter_emulated,np.var(np.diff(np.convolve(counter_emulator_baseline_quantized,np.ones(N, dtype=int),mode='valid')[0:((NoSamples+1)*N-1):N])))

poly_counter_emulator=LSNE(Nvalues_counter,sN_counter_emulated,[2,1,0])

plt.figure()
p1=plt.plot(Nvalues_counter,sN_counter,marker='o')
p2=plt.plot(Nvalues_counter,np.polyval(poly_counter,Nvalues_counter),color='red')
p1=plt.plot(Nvalues_counter,sN_counter_emulated,marker='o')
p2=plt.plot(Nvalues_counter,np.polyval(poly_counter_emulator,Nvalues_counter),color='red')

plt.legend(['Allan variance counter FPGA','fit FPGA','Allan variance counter emulator','fit emulator'])
plt.title('Allan variance counter')
plt.yscale('linear')
plt.xscale('linear')
plt.ylabel('variance counter')
plt.xlabel("accumulation in periods of RO0")
plt.savefig("allan_variance_counter_fpga.pdf")

plt.figure()
p1=plt.hist(counter_full,alpha=0.7)
p1=plt.hist(counter_emulator_baseline_quantized[:1*np.size(counter_full):1],alpha=0.7)

plt.legend(['counter FPGA','counter emulator'])
plt.title('Comparisson FPGA/emulator')
plt.yscale('linear')
plt.xscale('linear')
plt.ylabel('occurences')
plt.xlabel("counter value")
plt.savefig("counter_fpga_vs_emulator.pdf")

################################################################################################
##################       Entropy and autocorrelation ERO-TRNG    ###############################
################################################################################################


#bit ourput for ERO-TRNG
def ERO_bits(T1,T2,Ath,Afl,N,size):
    #N=frequency divider factor, T1/T2=period, 
    #RO2 with jitter RO1 perfect - transfer of jitter 
    #example: ERO_bits(2e-9,2e-9,Ath,Afl,1000,int(1e7))
       
    #generate noise 
    di_thermal = cn.powerlaw_psd_gaussian(0 ,size)
    di_flicker = cn.powerlaw_psd_gaussian(1,size)
    
    dti_emulator=N*T2*np.ones(size)+di_thermal*np.sqrt(Ath*T2*N/factor_th)+di_flicker*np.sqrt(Afl*((N*T2)**2)/factor_fl)

    ti_emulator=np.cumsum(dti_emulator)
    bits=np.round((ti_emulator/T1)%1)
    
    return bits

 

def H_Baudet(T1,T2,a1,N):
    #si a1 obtenu en avar t vs t alors T1=T2=1
    #if a1 obtained from Allan variance vs accumulation time (time vs time) then T1 & T2 cancel each other
    #Hmin=1-(4/(math.pi**2*math.log(2)))*math.exp((-2*math.pi**2)*a1*N/(T1**2*T2))
    Hmin=1-(4/(math.pi**2*math.log(2)))*math.exp((-2*math.pi**2)*a1*N/(T2))
    return Hmin

def entropy(bit_series, order):
    # Convert bit series to a sequence of integers (0 or 1)
    sequence = np.array([int(b) for b in bit_series])
    # Initialize count array for each possible state
    state_count = np.zeros(2**order)
    # Count occurrences of each state
    for i in range(order, len(sequence)):
        state = (np.packbits(sequence[i-order:i])/2**(8-order)).astype(int)
        state_count[state] += 1
    # Calculate probabilities and entropy
    prob = state_count / np.sum(state_count)
    prob=prob[np.nonzero(prob)]
    entropy = -np.sum(prob * np.log2(prob)/order)
    return entropy

def entropy_16bits(bit_series, order=16):
    # Convert bit series to a sequence of integers (0 or 1)
    sequence = np.array([int(b) for b in bit_series])
    # Initialize count array for each possible state
    state_count = np.zeros(2**order)
    # Count occurrences of each state
    for i in range(order, len(sequence)):
        state = (np.packbits(sequence[i-order:i].reshape(-1, 2, 8)[:, ::-1]).view(np.uint16)).astype(int)
        state_count[state] += 1
    # Calculate probabilities and entropy
    prob = state_count / np.sum(state_count)
    prob=prob[np.nonzero(prob)]
    entropy = -np.sum(prob * np.log2(prob)/order)
    return entropy

Ath=poly_VA_norm_2000[1]
Afl=poly_VA_norm_2000[0]

N_accumulation=np.logspace(2,5,34).astype(int)


H_Baudet_all=np.zeros(np.size(N_accumulation))
H8_all=np.zeros(np.size(N_accumulation))
H8_noFl_all=np.zeros(np.size(N_accumulation))
H8_noTh_all=np.zeros(np.size(N_accumulation))
H8_10Th_all=np.zeros(np.size(N_accumulation))
H8_10Fl_all=np.zeros(np.size(N_accumulation))
H8_100Th_all=np.zeros(np.size(N_accumulation))
H8_100Fl_all=np.zeros(np.size(N_accumulation))

for i in range(np.size(N_accumulation)):
    print('start for N='+str(N_accumulation[i]))
    H_Baudet_all[i]=H_Baudet(2e-9,2.00e-9, Ath, N_accumulation[i])
    H8_noFl_all[i]=entropy(ERO_bits(2e-9,2.00e-9,Ath,0*Afl,N_accumulation[i],int(1e6)),8)
    H8_noTh_all[i]=entropy(ERO_bits(2e-9,2.00e-9,0*Ath,Afl,N_accumulation[i],int(1e6)),8)
    H8_10Th_all[i]=entropy(ERO_bits(2e-9,2.00e-9,10*Ath,Afl,N_accumulation[i],int(1e6)),8)
    H8_10Fl_all[i]=entropy(ERO_bits(2e-9,2.00e-9,Ath,10*Afl,N_accumulation[i],int(1e6)),8)
    H8_100Th_all[i]=entropy(ERO_bits(2e-9,2.00e-9,100*Ath,Afl,N_accumulation[i],int(1e6)),8)
    H8_100Fl_all[i]=entropy(ERO_bits(2e-9,2.00e-9,Ath,100*Afl,N_accumulation[i],int(1e6)),8)
    
    
    H8_all[i]=entropy(ERO_bits(2e-9,2.00e-9,Ath,Afl,N_accumulation[i],int(1e6)),8)
    
plt.figure()
p1=plt.plot(N_accumulation,H_Baudet_all)
p1=plt.plot(N_accumulation,H8_all,marker='o')
p1=plt.plot(N_accumulation,H8_noFl_all,marker='v', linestyle='--')
p1=plt.plot(N_accumulation,H8_noTh_all,marker='v', linestyle='--')
p1=plt.plot(N_accumulation,H8_10Fl_all,marker='^', linestyle='--')
p1=plt.plot(N_accumulation,H8_10Th_all,marker='^', linestyle='--')
p1=plt.plot(N_accumulation,H8_100Fl_all,marker='s', linestyle='--')
p1=plt.plot(N_accumulation,H8_100Th_all,marker='s', linestyle='--')

plt.legend(['Hmin Baudet','H8 reference','H8 no th','H8 no fl','H8 10*th','H8 10*fl','H8 100*th','H8 100*fl',])
plt.title('Entropy')
plt.yscale('linear')
plt.xscale('log')
plt.ylabel('H')
plt.xlabel("accumulation factor N")
plt.savefig("entropy_h8.pdf")

plt.figure()
p1=plt.plot(N_accumulation,H_Baudet_all)
p1=plt.plot(N_accumulation,H8_all,marker='o')
p1=plt.plot(N_accumulation,H8_noFl_all,marker='v', linestyle='--')
p1=plt.plot(N_accumulation,H8_10Fl_all,marker='^', linestyle='--')
p1=plt.plot(N_accumulation,H8_100Fl_all,marker='s', linestyle='--')
p1=plt.plot([100,100000],[0.997,0.997])
p1=plt.plot([8373],[1.003],color='C1',marker='o', markersize=10)
p1=plt.plot([18329],[1.003],color='C2',marker='v', markersize=10)
p1=plt.plot([3934],[1.003],color='C3',marker='^', markersize=10)
p1=plt.plot([1626],[1.003],color='C4',marker='s', markersize=10)

plt.legend(['Hmin Baudet','H8','H8 no fl','H8 10*fl','H8 100*fl','0.997'])
plt.title('Entropy for different flicker noise amplitudes')
plt.yscale('linear')
plt.xscale('log')
plt.ylim((0.96,1.005))
plt.ylabel('H')
plt.xlabel("accumulation factor N")
plt.savefig("entropy_for_different_flicker_noise_amplitude.pdf")

plt.figure()
p1=plt.plot(N_accumulation,H_Baudet_all)
p1=plt.plot(N_accumulation,H8_all,marker='o')
p1=plt.plot(N_accumulation,H8_noTh_all,marker='v', linestyle='--')
p1=plt.plot(N_accumulation,H8_10Th_all,marker='^', linestyle='--')
p1=plt.plot(N_accumulation,H8_100Th_all,marker='s', linestyle='--')
p1=plt.plot([100,100000],[0.997,0.997])
p1=plt.plot([8373],[1.003],color='C1',marker='o', markersize=10)
p1=plt.plot([12011],[1.003],color='C2',marker='v', markersize=10)
p1=plt.plot([1821],[1.003],color='C3',marker='^', markersize=10)
p1=plt.plot([183],[1.003],color='C4',marker='s', markersize=10)

plt.legend(['Hmin Baudet','H8 reference','H8 no th','H8 10*th','H8 100*th','0.997'])
plt.title('Entropy for different thermal noise amplitudes')
plt.yscale('linear')
plt.xscale('log')
plt.ylim((0.96,1.005))
plt.ylabel('H')
plt.xlabel("accumulation factor N")
plt.savefig("entropy_for_different_thermal_noise_amplitude.pdf")

#entropy rate of an n-bit series

H2_all=np.zeros(np.size(N_accumulation))
H4_all=np.zeros(np.size(N_accumulation))
H6_all=np.zeros(np.size(N_accumulation))
H8_all=np.zeros(np.size(N_accumulation))
H16_all=np.zeros(np.size(N_accumulation))


for i in range(np.size(N_accumulation)):
    print('start for N='+str(N_accumulation[i]))
    H2_all[i]=entropy(ERO_bits(2e-9,2.00e-9,Ath,Afl,N_accumulation[i],int(1e6)),2)
    H4_all[i]=entropy(ERO_bits(2e-9,2.00e-9,Ath,Afl,N_accumulation[i],int(1e6)),4)
    H6_all[i]=entropy(ERO_bits(2e-9,2.00e-9,Ath,Afl,N_accumulation[i],int(1e6)),6)
    H8_all[i]=entropy(ERO_bits(2e-9,2.00e-9,Ath,Afl,N_accumulation[i],int(1e6)),8)
    H16_all[i]=entropy_16bits(ERO_bits(2e-9,2.00e-9,Ath,Afl,N_accumulation[i],int(1e6)),16)
    
plt.figure()
p1=plt.plot(N_accumulation,H2_all,marker='o')
p1=plt.plot(N_accumulation,H4_all,marker='o')
p1=plt.plot(N_accumulation,H6_all,marker='o')
p1=plt.plot(N_accumulation,H8_all,marker='o')
p1=plt.plot(N_accumulation,H16_all,marker='o')

plt.legend(['H2','H4','H6','H8','H16'])
plt.title('Entropy')
plt.yscale('linear')
plt.xscale('log')
plt.ylim((0.6,1.05))
plt.ylabel('H')
plt.xlabel("accumulation factor N")
plt.savefig("entropy_h2_to_h16.pdf")



####################################################################################################
########################################## Autocorrelation #########################################
####################################################################################################


    
def autocorrelation(bit_series,k_all=np.linspace(0,100,101).astype(int)):
    #bits = np.array(bit_series,dtype=int)
    bits = np.array(bit_series)
    autocorr = np.zeros(np.size(k_all))
    for k in k_all:
        print(k)
        autocorr[np.where(k_all==k)]=np.mean(np.multiply(bits[k:]-np.mean(bits[k:])*np.ones(np.size(bits[k:])),bits[:len(bits)-k]-np.mean(bits[:len(bits)-k])*np.ones(np.size(bits[:len(bits)-k]))))/np.sqrt(np.var(bits[:len(bits)-k])*np.var(bits[k:]))
        
    return autocorr


####################### Autocorrelation of the raw time series ######################################
size_ACF=int(1e7)
di_thermal_ACF = cn.powerlaw_psd_gaussian(0, size_ACF)
di_flicker_ACF = cn.powerlaw_psd_gaussian(1, size_ACF)
#kflicker=np.logspace(0,5,6)
kflicker=np.array([0,1,10,100,1000])
Nacc=100
#Nacc=np.logspace(0,6,7)
T_raw=2e-9

plt.figure()
for k in kflicker:
    dti_raw=N*T_raw*np.ones(size_ACF)+di_thermal_ACF*np.sqrt(Ath*T_raw*Nacc/factor_th)+di_flicker_ACF*np.sqrt(k*Afl*((Nacc*T_raw)**2)/factor_fl)

    dti_raw=np.diff(np.cumsum(dti_raw))
    
    ACAll_dti_test=autocorrelation(dti_raw,np.linspace(0,10000,101).astype(int))
    p1=plt.plot(np.linspace(0,10000,101).astype(int),ACAll_dti_test,marker='o')

plt.legend(['0*flicker','1*flicker','10*flicker','100*flicker','1000*flicker'])
plt.title('Autocorrelation N='+str(int(N)))
plt.yscale('linear')
plt.xscale('linear')
plt.ylabel('Autocorrelation')
plt.xlabel("k")
plt.savefig("autocorrelation_n={:d}.pdf".format(N))
    
####################### Autocorrelation of the bit series ######################################

#different flicker noise amplitudes

#AC_N1_fl0=autocorrelation(ERO_bits(2e-9,2e-9,Ath,0*Afl,1,int(1e7)))
#AC_N1_fl1=autocorrelation(ERO_bits(2e-9,2e-9,Ath,Afl,1,int(1e7)))
#AC_N1_fl10=autocorrelation(ERO_bits(2e-9,2e-9,Ath,10*Afl,1,int(1e7)))
#AC_N1_fl100=autocorrelation(ERO_bits(2e-9,2e-9,Ath,100*Afl,1,int(1e7)))
#AC_N1_fl1000=autocorrelation(ERO_bits(2e-9,2e-9,Ath,1000*Afl,1,int(1e7)))

AC_N100_fl0=autocorrelation(ERO_bits(2e-9,2e-9,Ath,0*Afl,100,int(1e7)))
AC_N100_fl1=autocorrelation(ERO_bits(2e-9,2e-9,Ath,Afl,100,int(1e7)))
AC_N100_fl10=autocorrelation(ERO_bits(2e-9,2e-9,Ath,10*Afl,100,int(1e7)))
AC_N100_fl100=autocorrelation(ERO_bits(2e-9,2e-9,Ath,100*Afl,100,int(1e7)))
AC_N100_fl1000=autocorrelation(ERO_bits(2e-9,2e-9,Ath,1000*Afl,100,int(1e7)))

AC_N1000_fl0=autocorrelation(ERO_bits(2e-9,2e-9,Ath,0*Afl,1000,int(1e7)))
AC_N1000_fl1=autocorrelation(ERO_bits(2e-9,2e-9,Ath,Afl,1000,int(1e7)))
AC_N1000_fl10=autocorrelation(ERO_bits(2e-9,2e-9,Ath,10*Afl,1000,int(1e7)))
AC_N1000_fl100=autocorrelation(ERO_bits(2e-9,2e-9,Ath,100*Afl,1000,int(1e7)))
AC_N1000_fl1000=autocorrelation(ERO_bits(2e-9,2e-9,Ath,1000*Afl,1000,int(1e7)))

#plt.figure()
#p1=plt.plot(np.round(np.linspace(0,100,101)).astype(int)[:-1],AC_N1_fl0[1:],marker='o')
#p1=plt.plot(np.round(np.linspace(0,100,101)).astype(int)[:-1],AC_N1_fl1[1:],marker='o')
#p1=plt.plot(np.round(np.linspace(0,100,101)).astype(int)[:-1],AC_N1_fl10[1:],marker='o')
#p1=plt.plot(np.round(np.linspace(0,100,101)).astype(int)[:-1],AC_N1_fl100[1:],marker='o')
#p1=plt.plot(np.round(np.linspace(0,100,101)).astype(int)[:-1],AC_N1_fl1000[1:],marker='o')

#plt.legend(['N=1 0*flicker','N=1 1*flicker','N=1 10*flicker','N=1 100*flicker','N=1 1000*flicker','N=1 10000*flicker','N=1 100000*flicker',])
#plt.title('Autocorrelation ')
#plt.yscale('linear')
#plt.xscale('linear')
#plt.ylabel('Autocorrelation')
#plt.xlabel("k")
#plt.savefig("autocorrelation.pdf")

plt.figure()
p1=plt.plot(np.round(np.linspace(0,100,101)).astype(int)[:-1],AC_N100_fl0[1:],marker='o')
p1=plt.plot(np.round(np.linspace(0,100,101)).astype(int)[:-1],AC_N100_fl1[1:],marker='o')
p1=plt.plot(np.round(np.linspace(0,100,101)).astype(int)[:-1],AC_N100_fl10[1:],marker='o')
p1=plt.plot(np.round(np.linspace(0,100,101)).astype(int)[:-1],AC_N100_fl100[1:],marker='o')
p1=plt.plot(np.round(np.linspace(0,100,101)).astype(int)[:-1],AC_N100_fl1000[1:],marker='o')

plt.legend(['N=100 0*flicker','N=100 1*flicker','N=100 10*flicker','N=100 100*flicker','N=100 1000*flicker','N=100 10000*flicker','N=100 100000*flicker',])
plt.title('Autocorrelation ')
plt.yscale('linear')
plt.xscale('linear')
plt.ylabel('Autocorrelation')
plt.xlabel("k")
plt.savefig("autocorrelation_for_different_n.pdf")

plt.figure()
p1=plt.plot(np.round(np.linspace(0,100,101)).astype(int)[:-1],AC_N1000_fl0[1:],marker='o')
p1=plt.plot(np.round(np.linspace(0,100,101)).astype(int)[:-1],AC_N1000_fl1[1:],marker='o')
p1=plt.plot(np.round(np.linspace(0,100,101)).astype(int)[:-1],AC_N1000_fl10[1:],marker='o')
p1=plt.plot(np.round(np.linspace(0,100,101)).astype(int)[:-1],AC_N1000_fl100[1:],marker='o')
p1=plt.plot(np.round(np.linspace(0,100,101)).astype(int)[:-1],AC_N1000_fl1000[1:],marker='o')

plt.legend(['N=10000 0*flicker','N=10000 1*flicker','N=10000 10*flicker','N=10000 100*flicker','N=10000 1000*flicker','N=10000 10000*flicker','N=10000 100000*flicker',])
plt.title('Autocorrelation ')
plt.yscale('linear')
plt.xscale('linear')
plt.ylabel('Autocorrelation')
plt.xlabel("k")
plt.savefig("autocorrelation_for_different_n_extended.pdf")

#different thermal noise amplitudes

AC_N100_th0=autocorrelation(ERO_bits(2e-9,2e-9,0*Ath,Afl,100,int(1e7)))
AC_N100_th1=autocorrelation(ERO_bits(2e-9,2e-9,Ath,Afl,100,int(1e7)))
AC_N100_th10=autocorrelation(ERO_bits(2e-9,2e-9,10*Ath,Afl,100,int(1e7)))
AC_N100_th100=autocorrelation(ERO_bits(2e-9,2e-9,100*Ath,Afl,100,int(1e7)))
AC_N100_th1000=autocorrelation(ERO_bits(2e-9,2e-9,1000*Ath,Afl,100,int(1e7)))
print('N100')

plt.figure()
p1=plt.plot(np.round(np.linspace(0,100,101)).astype(int)[:-1],AC_N100_th0[1:],marker='o')
p1=plt.plot(np.round(np.linspace(0,100,101)).astype(int)[:-1],AC_N100_th1[1:],marker='o')
p1=plt.plot(np.round(np.linspace(0,100,101)).astype(int)[:-1],AC_N100_th10[1:],marker='o')
p1=plt.plot(np.round(np.linspace(0,100,101)).astype(int)[:-1],AC_N100_th100[1:],marker='o')
p1=plt.plot(np.round(np.linspace(0,100,101)).astype(int)[:-1],AC_N100_th1000[1:],marker='o')

plt.legend(['N=100 0*thermal','N=100 1*thermal','N=100 10*thermal','N=100 100*thermal','N=100 1000*thermal','N=100 10000*thermal','N=100 100000*thermal',])
plt.title('Autocorrelation ')
plt.yscale('linear')
plt.xscale('linear')
plt.ylabel('Autocorrelation')
plt.xlabel("k")
plt.savefig("autocorrelation_for_different_thermal_amplitudes.pdf")

####################################ACF for different TRO1######################################
AC_N100_fl1_T0=autocorrelation(ERO_bits(2e-9,2e-9,Ath,Afl,100,int(1e7)),np.linspace(0,100,101).astype(int))
AC_N100_fl1_2T0=autocorrelation(ERO_bits(2*2e-9,2e-9,Ath,Afl,100,int(1e7)),np.linspace(0,1000,101).astype(int))
AC_N100_fl1_4T0=autocorrelation(ERO_bits(4*2e-9,2e-9,Ath,Afl,100,int(1e7)),np.linspace(0,1000,101).astype(int))
#AC_N100_fl1_7T0=autocorrelation(ERO_bits(7*2e-9,2e-9,Ath,Afl,100,int(1e7)),np.linspace(0,1000,101).astype(int))
AC_N100_fl1_10T0=autocorrelation(ERO_bits(10*2e-9,2e-9,Ath,Afl,100,int(1e7)),np.linspace(0,1000,101).astype(int))
#AC_N100_fl1_20T0=autocorrelation(ERO_bits(20*2e-9,2e-9,Ath,Afl,100,int(1e7)),np.linspace(0,10000,1001).astype(int))
#AC_N100_fl1_40T0=autocorrelation(ERO_bits(40*2e-9,2e-9,Ath,Afl,100,int(1e7)),np.linspace(0,10000,1001).astype(int))
#AC_N100_fl1_70T0=autocorrelation(ERO_bits(70*2e-9,2e-9,Ath,Afl,100,int(1e7)),np.linspace(0,10000,1001).astype(int))
#AC_N100_fl1_100T0=autocorrelation(ERO_bits(100*2e-9,2e-9,Ath,Afl,100,int(1e7)),np.linspace(0,10000,1001).astype(int))
#AC_N100_fl1_200T0=autocorrelation(ERO_bits(200*2e-9,2e-9,Ath,Afl,100,int(1e7)),np.linspace(0,10000,101).astype(int))
#AC_N100_fl1_400T0=autocorrelation(ERO_bits(400*2e-9,2e-9,Ath,Afl,100,int(1e7)),np.linspace(0,10000,101).astype(int))
#AC_N100_fl1_700T0=autocorrelation(ERO_bits(700*2e-9,2e-9,Ath,Afl,100,int(1e7)),np.linspace(0,10000,101).astype(int))
#AC_N100_fl1_1000T0=autocorrelation(ERO_bits(1000*2e-9,2e-9,Ath,Afl,100,int(1e7)),np.linspace(0,10000,101).astype(int))

plt.figure()
p1=plt.plot(np.linspace(0,100,101).astype(int),np.abs(AC_N100_fl1_T0),marker='o')
p1=plt.plot(np.linspace(0,1000,101).astype(int),np.abs(AC_N100_fl1_2T0),marker='o')
p1=plt.plot(np.linspace(0,1000,101).astype(int),np.abs(AC_N100_fl1_4T0),marker='o')
#p1=plt.plot(np.linspace(0,1000,101).astype(int),np.abs(AC_N100_fl1_7T0),marker='o')
p1=plt.plot(np.linspace(0,1000,101).astype(int),np.abs(AC_N100_fl1_10T0),marker='o')
#p1=plt.plot(np.linspace(0,1000,1001).astype(int),np.abs(AC_N100_fl1_20T0),marker='o')
#p1=plt.plot(np.linspace(0,10000,1001).astype(int),np.abs(AC_N100_fl1_40T0),marker='o')
#p1=plt.plot(np.linspace(0,10000,1001).astype(int),np.abs(AC_N100_fl1_70T0),marker='o')
#p1=plt.plot(np.linspace(0,10000,1001).astype(int),np.abs(AC_N100_fl1_100T0),marker='o')
#p1=plt.plot(np.linspace(0,10000,101).astype(int),np.abs(AC_N100_fl1_1000T0),marker='o')

plt.legend(['1*T0RO1','2*T0RO1','4*T0RO1','7*T0RO1','10*T0RO1','20*T0RO1','40*T0RO1','70*T0RO1','100*T0RO1'])
plt.title('Autocorrelation N=100')
plt.yscale('log')
plt.xscale('log')
plt.ylim([1e-2,1.1])
plt.ylabel('Autocorrelation')
plt.xlabel("lag k")
plt.savefig("autocorrelation_for_different_tro1.pdf")


########## absolute jitter for different flicker noise amplitudes ###############################
size_ACF=int(1e5)
T_test=2e-9

di_thermal_ACF = cn.powerlaw_psd_gaussian(0, size_ACF)
di_flicker_ACF = cn.powerlaw_psd_gaussian(1, size_ACF)

dti_test=T_test*np.ones(size_ACF)+di_thermal_ACF*np.sqrt(Ath*T_test/factor_th)+di_flicker_ACF*np.sqrt(Afl*(T_test**2)/factor_fl)
ti_test=np.cumsum(dti_test)-2e-9*np.linspace(1,np.size(dti_test),np.size(dti_test))

dti_test2=T_test*np.ones(size_ACF)+di_thermal_ACF*np.sqrt(Ath*T_test/factor_th)+di_flicker_ACF*np.sqrt(10*Afl*(T_test**2)/factor_fl)
ti_test2=np.cumsum(dti_test2)-2e-9*np.linspace(1,np.size(dti_test),np.size(dti_test))

plt.figure()
p1=plt.plot(ti_test)
p1=plt.plot(ti_test2)

min_abs_jitter=np.min([np.min(ti_test),np.min(ti_test2)])
min_abs_jitter_qnt=sampling(min_abs_jitter,2e-9/2)
max_abs_jitter=np.max([np.max(ti_test),np.max(ti_test2)])
max_abs_jitter_qnt=sampling(max_abs_jitter,2e-9/2)
domains=int((max_abs_jitter_qnt-min_abs_jitter_qnt)/(2e-9/2))
levels=min_abs_jitter_qnt*np.ones(domains+1)+np.linspace(0, domains,domains+1)*(2e-9/2)
for k in range(domains):
    p1=plt.plot([1,size_ACF],[levels[k],levels[k]],color='gray',linestyle='--')

plt.legend(['1*flicker','10*flicker'])
plt.yscale('linear')
plt.xscale('linear')
plt.ylabel('Absolute jitter (s)')
plt.xlabel("Period #")
plt.savefig("absolute_jitter.pdf")
