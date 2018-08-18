import matplotlib.pyplot as plt
import numpy as np

Vm = 100
thetav = 0   # Volage amplitude and phase angle
Z = 1.25
gama = 60   # Impedance magnitude and phase angle
thetai = thetav - gama
theta = (thetav-thetai)*np.pi/180   # Degree to radian
Im = Vm/Z                           # Current amplitude
wt = np.arange(0.0, 2*np.pi, 0.05)  # wt from 0 to 2*pi
v = Vm*np.cos(wt)                   # Instantaneous voltage
i = Im*np.cos(wt+thetai*np.pi/180)  # Instantaneous current
p = v*i                             # Instantaneous power
V = Vm/np.sqrt(2)                   # rms voltage and current
I = Im/np.sqrt(2)
P = V*I*np.cos(theta)               # Average power
Q = V*I*np.sin(theta)               # Reactive power
S = P+1j*Q                          # Complex power
pr = P*(1+np.cos(2*(wt+thetav)))    # Eq. (2.6)
px = Q*np.sin(2*(wt+thetav))        # Eq. (2.8)
PP = P*np.ones(len(wt))             # Average power of length w for plot
xline = np.zeros(len(wt))           # Generates a zero vector
wt = wt*180/np.pi                   # Converting radians to degree

# Creating the figure
fig, axes = plt.subplots(2, 2)
for line in [v, i, xline]:
    axes[0, 0].plot(wt, line)
axes[0, 0].set_title("v(t)=Vm cos(wt), i(t)=Im cos(wt+{0})".format(thetai))
axes[0, 0].set_xlabel("wt, degree")

for line in [p, xline]:
    axes[0, 1].plot(wt, line)
axes[0, 1].set_title("p(t)=v(t)*i(t)")
axes[0, 1].set_xlabel("wt, degree")

for line in [pr, PP, xline]:
    axes[1, 0].plot(wt, line)
axes[1, 0].set_title("pr(t)  -  Eq. 2.6")
axes[1, 0].set_xlabel("wt, degree")

for line in [px, xline]:
    axes[1, 1].plot(wt, line)
axes[1, 1].set_title("px(t)  -  Eq. 2.8")
axes[1, 1].set_xlabel("wt, degree")

fig.subplots_adjust(
    wspace=0.3, hspace=0.5)  # Legger til mellomrom mellom figurer.
#plt.show()


# ######### Example 2.2 - Complex power balance ########
print("Example 2.2")
V = 1200 + 1j * 0  # Input Voltage
Z1 = 60 + 1j * 0   # Input impedance 1
Z2 = 6 + 1j * 12   # Input impedance 2
Z3 = 30 - 1j * 30  # Input impedance 3

I1, I2, I3 = [V / Z for Z in [Z1, Z2, Z3]]           # I = V/Z
print("I: ", I1, I2, I3)
S1, S2, S3 = [V * np.conj(I) for I in [I1, I2, I3]]  # S = U*I'
print("S: ", S1, S2, S3)
print("S sum:", S1+S2+S3)


# ######### Example 2.3 - Power factor correction #######
print("Example 2.3")
V = 200 + 1j * 0
Z1 = 100 + 1j * 0
Z2 = 10 + 1j * 20
# Find total real power, total reactive power, PF, total current
S1, S2 = [V**2 / np.conj(Z) for Z in [Z1, Z2]]
Stot = S1 + S2
print(Stot)
Ptot = np.real(Stot)
Qtot = np.imag(Stot)
Itot = np.conj(Stot) / np.conj(V)
print("Total real P: ", Ptot)
print("Total reactive Q: ", Qtot)
print("Total current I:", Itot)
print("PF: ", np.cos(np.angle(Stot)))

# Find the capacitance of the capacitor connected across the loads to
# improve the overall power factor to 0.8 lagging.

# Total power is equal, both for 0,6 and 0,8 power factor.
newTheta = np.arccos(0.8)
newQtot = Ptot * np.tan(newTheta)
Qc = Qtot - newQtot
Zc = V**2 / (1j * Qc)
C = 10**6 / (2 * np.pi * 50 * Zc)           # myF: 10E6/2pi*f*Z
print("New theta", np.arccos(0.8), np.arccos(0.8) * 180 / np.pi)
print("New total Q: ", newQtot)
print("Zc: ", Zc)
print("Installed C:{0} microF".format(C))

Snew = Ptot + 1j * newQtot
Inew = np.conj(Snew) / np.conj(V)
print("New S:", Snew)
print("New I:{0}, ({1} v{2})".format(Inew, np.abs(Inew),
                                     np.angle(Inew) * 180 / np.pi))


# ###### Example 2.4 #########
print("Example 2.4")
V = 1400 + 1j * 0
freq = 60
# Load 1, inductiv, 125 kva at 0.28 pf
# Load 2, capacitiv, 10 kw and 40 kvar
# Load 3, Resistive load of 15 kwar
kv = 1000
theta = np.arccos(0.28)
S1 = (125*np.cos(theta)+1j*125*np.sin(theta))*kv
S2 = (10-1j*40)*kv
S3 = 15 * kv
print("S1:", S1)
print("S2:", S2)
print("S3:", S3)

# Find the total kW, kvar, kVA and the supply PF
Stot = S1 + S2 + S3
print("Total kW: ", np.real(Stot))
print("Total kvar: ", np.imag(Stot))
print("Total kVA: ", np.abs(Stot))
print("PF: ", np.real(Stot) / np.abs(Stot))

# Determine the kvar ratoing of a capacitor and the capacitance in
# microF to improve the total power factor to be 0.8 lagging.
newTheta = np.arccos(0.8)   # Finding the new angle with PF = 0.8
newQtot = np.real(Stot)*np.tan(newTheta)  # Finding the new total Q with PF 0.8
Qc = np.imag(Stot)-newQtot  # Installed Qc = Qold - Qnew
Xc = V**2/np.conj(0-1j*Qc)  # Convert to Z
C = 10**6/(2*np.pi*freq*np.abs(Xc))  # Convert from Z to C

print("newTheta: ", newTheta, newTheta*180/np.pi)
print("newQtot: ", newQtot)
print("Qc kvar: ", Qc)
print("Xc Ohms: ", Xc)
print("C microF: ", C)


# ########## Example 2.5 #################
# Two voltage sources over a line Z      #
def polarToCart(arg, angle=0, deg=False):
    if deg:
        angle = angle * np.pi / 180
    a = arg * (np.cos(angle) + 1j * np.sin(angle))
    return a


def cartToPolar(z, deg=False):
    scale = 1
    if deg:
        scale = 180 / np.pi
    a = [np.abs(z), np.angle(z) * scale]
    return a


V1 = polarToCart(120, -5, deg=True)
V2 = 100
Z = 1 + 1j * 7

I12 = (V1 - V2) / Z
I21 = (V2 - V1) / Z
print("I12: ", cartToPolar(I12, deg=True))
print("I21: ", cartToPolar(I21, deg=True))
S12 = V1 * np.conj(I12)
S21 = V2 * np.conj(I21)
print("S12: {0} W = {1} W".format(S12, cartToPolar(S12, deg=True)))
print("S21: {0} W = {1} W".format(S21, cartToPolar(S21, deg=True)))

# Line losses
Sloss = S12+S21
print("Sloss: ",Sloss)
